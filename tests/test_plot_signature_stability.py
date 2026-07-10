import csv
import subprocess
import sys
from pathlib import Path

import mutational_signature.plot_signature_stability as pss


ROOT = Path(__file__).resolve().parents[1]


def write_summary(path):
  fieldnames = [
    "Sample",
    "Signature",
    "point_exposure_proportion",
    "bootstrap_percentile_2_5_proportion",
    "bootstrap_percentile_97_5_proportion",
    "detection_frequency",
  ]
  rows = [
    ("S1", "SBS1", "0.4", "0.3", "0.5", "1.0"),
    ("S2", "SBS1", "0.2", "0.1", "0.3", "0.97"),
    ("S3", "SBS1", "0.1", "0.0", "0.2", "0.80"),
    ("S1", "SBS2", "0.3", "0.2", "0.4", "0.60"),
    ("S2", "SBS2", "0.0", "0.0", "0.1", "0.20"),
  ]
  with open(path, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(fieldnames)
    writer.writerows(rows)


def test_detection_category_boundaries():
  assert pss.detection_category(0.49)[0] == "<50%"
  assert pss.detection_category(0.50)[0] == "50-75%"
  assert pss.detection_category(0.75)[0] == "75-95%"
  assert pss.detection_category(0.95)[0] == "95-99%"
  assert pss.detection_category(0.99)[0] == ">99%"
  assert pss.detection_category(1.0)[0] == ">99%"


def test_plot_signature_stability_cli_smoke(tmp_path):
  summary = tmp_path / "summary.tsv"
  out_dir = tmp_path / "plots"
  write_summary(summary)
  cmd = [
    sys.executable,
    "-m",
    "mutational_signature.plot_signature_stability",
    "--input",
    str(summary),
    "--output-dir",
    str(out_dir),
    "--include",
    "SBS1",
    "--dpi",
    "80",
  ]
  subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=True)
  plot = out_dir / "SBS1.signature_stability.png"
  index = out_dir / "index.tsv"
  assert plot.exists()
  assert plot.stat().st_size > 0
  assert index.exists()
  rows = list(csv.DictReader(index.open(), delimiter="\t"))
  assert rows[0]["Signature"] == "SBS1"
  assert rows[0]["Samples"] == "3"
