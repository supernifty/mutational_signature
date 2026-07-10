import csv
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import mutational_signature.signature_stability as ss


ROOT = Path(__file__).resolve().parents[1]


def write_fixture(tmp_path):
  signatures = tmp_path / "signatures.tsv"
  counts = tmp_path / "counts.tsv"
  signatures.write_text(
    "Sig\tA>A\tC>C\tG>G\n"
    "SIG_A\t1\t0\t0\n"
    "SIG_B\t0\t1\t0\n"
    "SIG_C\t0.48\t0.48\t0.04\n"
  )
  counts.write_text("Variation\tCount\nA>A\t90\nC>C\t10\nG>G\t0\n")
  return signatures, counts


def read_rows(text):
  return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def base_args(signatures, counts, **overrides):
  values = {
    "signatures": str(signatures),
    "counts": str(counts),
    "sample": "S1",
    "replicates": 20,
    "perturbation": "dirichlet-multinomial",
    "alpha": 0.1,
    "fraction": None,
    "mutation_count": None,
    "seed": 7,
    "metric": "cosine",
    "solver": "basin",
    "max_sigs": None,
    "context_cutoff": 1e6,
    "mode": "direct",
    "refine_minimum_error": 0.01,
    "interval": 0.95,
    "detection_threshold_mutations": 0,
    "detection_threshold_proportion": 0,
    "exceeds_threshold_mutations": 0,
    "summary_output": None,
    "replicates_output": None,
    "count_column": "Count",
  }
  values.update(overrides)
  return SimpleNamespace(**values)


def test_perturbation_sample_sizes():
  counts = np.array([9, 1, 0])
  rng = np.random.default_rng(1)
  assert np.sum(ss.perturb_counts(counts, "multinomial", rng)) == 10
  assert np.sum(ss.perturb_counts(counts, "dirichlet-multinomial", rng, alpha=0.1)) == 10
  assert np.sum(ss.perturb_counts(counts, "downsample", rng, fraction=0.5)) == 5
  assert np.sum(ss.perturb_counts(counts, "downsample", rng, mutation_count=4)) == 4


def test_reproducible_with_fixed_seed(tmp_path):
  signatures, counts = write_fixture(tmp_path)
  args = base_args(signatures, counts)
  rows1, reps1 = ss.run(args, out=None)
  rows2, reps2 = ss.run(args, out=None)
  assert rows1 == rows2
  assert reps1 == reps2


def test_mixed_context_file_ignores_nonmatching_contexts(tmp_path):
  signatures, counts = write_fixture(tmp_path)
  counts.write_text("Variation\tCount\nA>A\t90\nC>C\t10\nDEL_T_1_4\t12\n")
  rows, reps = ss.run(base_args(signatures, counts, replicates=3), out=None)
  assert len(rows) == 3
  assert {row["Signature"] for row in reps} == {"SIG_A", "SIG_B", "SIG_C"}
  assert {row["total_mutations"] for row in reps} == {100}


def test_fails_when_no_positive_contexts_match(tmp_path):
  signatures, counts = write_fixture(tmp_path)
  counts.write_text("Variation\tCount\nDEL_T_1_4\t12\nINS_C_1_3\t5\n")
  with pytest.raises(ValueError, match="total count must be greater than zero"):
    ss.run(base_args(signatures, counts, replicates=3), out=None)


def test_missing_signatures_recorded_as_zero(tmp_path):
  signatures, counts = write_fixture(tmp_path)
  replicate_path = tmp_path / "replicates.tsv"
  args = base_args(signatures, counts, mode="refine", replicates_output=str(replicate_path), replicates=5)
  ss.run(args, out=None)
  rows = list(csv.DictReader(replicate_path.open(), delimiter="\t"))
  sigs_by_rep = {}
  for row in rows:
    sigs_by_rep.setdefault(row["Replicate"], set()).add(row["Signature"])
  assert all(sigs == {"SIG_A", "SIG_B", "SIG_C"} for sigs in sigs_by_rep.values())
  assert any(float(row["exposure_mutations"]) == 0 for row in rows)


def test_detection_frequency_and_interval_known_values():
  point = {
    "signature_mutations": np.array([5.0]),
    "signature_proportions": np.array([0.5]),
    "error": (0.1, None),
  }
  rows = [
    {"Signature": "SIG_A", "exposure_mutations": 0.0, "exposure_proportion": 0.0, "detected": "false", "reconstruction_error": 0.1},
    {"Signature": "SIG_A", "exposure_mutations": 10.0, "exposure_proportion": 1.0, "detected": "true", "reconstruction_error": 0.2},
    {"Signature": "SIG_A", "exposure_mutations": 20.0, "exposure_proportion": 2.0, "detected": "true", "reconstruction_error": 0.3},
  ]
  summary = ss.summarise("S1", np.array(["SIG_A"]), point, rows, 0.8, "multinomial", 1, "direct", "cosine", "basin", 0)[0]
  assert summary["detection_frequency"] == 2 / 3
  assert summary["detected_replicates"] == 2
  assert summary["bootstrap_percentile_10_mutations"] == pytest.approx(2.0)
  assert summary["bootstrap_percentile_90_mutations"] == pytest.approx(18.0)


def test_pure_signature_recovered_stably_high_count(tmp_path):
  signatures = tmp_path / "pure_signatures.tsv"
  counts = tmp_path / "pure_counts.tsv"
  signatures.write_text("Sig\tA>A\tC>C\nSIG_A\t1\t0\nSIG_B\t0\t1\n")
  counts.write_text("Variation\tCount\nA>A\t1000\nC>C\t0\n")
  rows, _ = ss.run(base_args(signatures, counts, replicates=20, perturbation="multinomial"), out=None)
  sig_a = next(r for r in rows if r["Signature"] == "SIG_A")
  assert sig_a["point_exposure_proportion"] > 0.99
  assert sig_a["bootstrap_percentile_2_5_proportion"] > 0.99


def test_weak_signature_less_stable_low_count(tmp_path):
  signatures = tmp_path / "weak_signatures.tsv"
  high = tmp_path / "high.count"
  low = tmp_path / "low.count"
  signatures.write_text("Sig\tA>A\tC>C\nSIG_A\t1\t0\nSIG_B\t0\t1\n")
  high.write_text("Variation\tCount\nA>A\t900\nC>C\t100\n")
  low.write_text("Variation\tCount\nA>A\t9\nC>C\t1\n")
  high_rows, _ = ss.run(base_args(signatures, high, replicates=50, seed=2), out=None)
  low_rows, _ = ss.run(base_args(signatures, low, replicates=50, seed=2), out=None)
  high_b = next(r for r in high_rows if r["Signature"] == "SIG_B")
  low_b = next(r for r in low_rows if r["Signature"] == "SIG_B")
  assert low_b["bootstrap_sd_exposure_proportion"] > high_b["bootstrap_sd_exposure_proportion"]


def test_similar_signatures_show_instability(tmp_path):
  signatures = tmp_path / "similar_signatures.tsv"
  counts = tmp_path / "similar.count"
  signatures.write_text("Sig\tA>A\tC>C\nSIG_A\t0.51\t0.49\nSIG_B\t0.49\t0.51\n")
  counts.write_text("Variation\tCount\nA>A\t10\nC>C\t10\n")
  rows, _ = ss.run(base_args(signatures, counts, replicates=40, seed=3), out=None)
  freqs = [r["detection_frequency"] for r in rows]
  assert any(f < 1.0 for f in freqs)


def test_direct_and_refine_cli_run(tmp_path):
  signatures, counts = write_fixture(tmp_path)
  for mode in ("direct", "refine"):
    cmd = [
      sys.executable,
      "-m",
      "mutational_signature.signature_stability",
      "--signatures",
      str(signatures),
      "--counts",
      str(counts),
      "--sample",
      "S1",
      "--replicates",
      "3",
      "--mode",
      mode,
      "--seed",
      "9",
    ]
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=True)
    rows = read_rows(result.stdout)
    assert {row["Signature"] for row in rows} == {"SIG_A", "SIG_B", "SIG_C"}
