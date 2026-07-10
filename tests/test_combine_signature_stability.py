import csv
import subprocess
import sys
from pathlib import Path

import pytest

import mutational_signature.combine_signature_stability as css


ROOT = Path(__file__).resolve().parents[1]


HEADER = [
  "Sample",
  "Signature",
  "point_exposure_mutations",
  "point_exposure_proportion",
  "bootstrap_percentile_2_5_mutations",
  "bootstrap_percentile_97_5_mutations",
  "bootstrap_percentile_2_5_proportion",
  "bootstrap_percentile_97_5_proportion",
  "detection_frequency",
  "total_mutations",
  "perturbation_method",
  "seed",
  "fitting_mode",
  "metric",
  "solver",
  "extra_column",
]


def write_summary(path, sample, rows, header=None):
  header = header or HEADER
  with open(path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=header, lineterminator="\n")
    writer.writeheader()
    for row in rows:
      base = {column: "" for column in header}
      base.update({
        "Sample": sample,
        "Signature": row["Signature"],
        "point_exposure_mutations": row.get("point_exposure_mutations", "10"),
        "point_exposure_proportion": row.get("point_exposure_proportion", "0.1"),
        "bootstrap_percentile_2_5_mutations": row.get("bootstrap_percentile_2_5_mutations", "1"),
        "bootstrap_percentile_97_5_mutations": row.get("bootstrap_percentile_97_5_mutations", "20"),
        "bootstrap_percentile_2_5_proportion": row.get("bootstrap_percentile_2_5_proportion", "0.01"),
        "bootstrap_percentile_97_5_proportion": row.get("bootstrap_percentile_97_5_proportion", "0.2"),
        "detection_frequency": row.get("detection_frequency", "1"),
        "total_mutations": row.get("total_mutations", "100"),
        "perturbation_method": "dirichlet-multinomial",
        "seed": "123",
        "fitting_mode": "direct",
        "metric": "cosine",
        "solver": "basin",
        "extra_column": row.get("extra_column", "extra"),
      })
      writer.writerow({column: base.get(column, "") for column in header})


def read_rows(text):
  return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def test_default_long_combines_key_columns(tmp_path, capsys):
  a = tmp_path / "a.tsv"
  b = tmp_path / "b.tsv"
  write_summary(a, "S1", [{"Signature": "SBS1"}])
  write_summary(b, "S2", [{"Signature": "SBS1"}])
  args = css.build_parser().parse_args(["--files", str(a), str(b)])
  css.run(args)
  rows = read_rows(capsys.readouterr().out)
  columns = list(rows[0].keys())
  assert [row["Sample"] for row in rows] == ["S1", "S2"]
  assert "point_exposure_proportion" in rows[0]
  assert "bootstrap_percentile_2_5_proportion" in rows[0]
  assert "bootstrap_percentile_97_5_proportion" in rows[0]
  point_index = columns.index("point_exposure_proportion")
  assert columns[point_index + 1:point_index + 3] == [
    "bootstrap_percentile_2_5_proportion",
    "bootstrap_percentile_97_5_proportion",
  ]
  assert "detection_frequency" in rows[0]
  assert "extra_column" not in rows[0]
  assert rows[0]["Source"] == str(a)


def test_all_columns_keeps_extra_columns(tmp_path, capsys):
  path = tmp_path / "a.tsv"
  write_summary(path, "S1", [{"Signature": "SBS1", "extra_column": "kept"}])
  args = css.build_parser().parse_args(["--files", str(path), "--all-columns"])
  css.run(args)
  rows = read_rows(capsys.readouterr().out)
  assert rows[0]["extra_column"] == "kept"


def test_wide_output(tmp_path, capsys):
  path = tmp_path / "a.tsv"
  write_summary(path, "S1", [{"Signature": "SBS1"}, {"Signature": "SBS2", "point_exposure_proportion": "0.2"}])
  args = css.build_parser().parse_args(["--files", str(path), "--wide"])
  css.run(args)
  rows = read_rows(capsys.readouterr().out)
  assert rows[0]["Sample"] == "S1"
  assert rows[0]["SBS1_point_exposure_proportion"] == "0.1"
  assert rows[0]["SBS2_point_exposure_proportion"] == "0.2"
  assert rows[0]["SBS1_detection_frequency"] == "1"


def test_include_and_exclude_filters(tmp_path, capsys):
  path = tmp_path / "a.tsv"
  write_summary(path, "S1", [{"Signature": "SBS1"}, {"Signature": "SBS2"}, {"Signature": "SBS3"}])
  args = css.build_parser().parse_args(["--files", str(path), "--include", "SBS1", "SBS2", "--exclude", "SBS2"])
  css.run(args)
  rows = read_rows(capsys.readouterr().out)
  assert [row["Signature"] for row in rows] == ["SBS1"]


def test_duplicate_sample_signature_errors(tmp_path):
  a = tmp_path / "a.tsv"
  b = tmp_path / "b.tsv"
  write_summary(a, "S1", [{"Signature": "SBS1"}])
  write_summary(b, "S1", [{"Signature": "SBS1"}])
  args = css.build_parser().parse_args(["--files", str(a), str(b)])
  with pytest.raises(ValueError, match="duplicate Sample\\+Signature"):
    css.run(args)


def test_non_95_interval_columns_are_detected(tmp_path, capsys):
  header = [
    "Sample",
    "Signature",
    "point_exposure_proportion",
    "bootstrap_percentile_10_proportion",
    "bootstrap_percentile_90_proportion",
    "detection_frequency",
  ]
  path = tmp_path / "a.tsv"
  write_summary(path, "S1", [{"Signature": "SBS1"}], header=header)
  args = css.build_parser().parse_args(["--files", str(path)])
  css.run(args)
  rows = read_rows(capsys.readouterr().out)
  assert "bootstrap_percentile_10_proportion" in rows[0]
  assert "bootstrap_percentile_90_proportion" in rows[0]


def test_cli_help_runs():
  cmd = [sys.executable, "-m", "mutational_signature.combine_signature_stability", "--help"]
  result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=True)
  assert "combine signature-stability summary" in result.stdout.lower()
