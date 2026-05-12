import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "example"


def run_command(args, stdin_text=None):
  return subprocess.run(
    args,
    input=stdin_text,
    text=True,
    capture_output=True,
    check=True,
    cwd=ROOT,
  )


def test_count_decompose_and_plot_smoke(tmp_path):
  count_file = tmp_path / "mini.count"
  exposures_file = tmp_path / "mini.exposures"
  plot_file = tmp_path / "mini.png"

  count_cmd = [
    sys.executable,
    "-m",
    "mutational_signature.count",
    "--genome",
    str(EXAMPLE / "mini_reference.fa"),
    "--vcf",
    str(EXAMPLE / "mini_sample.vcf"),
  ]
  count_result = run_command(count_cmd)
  count_file.write_text(count_result.stdout)

  rows = list(csv.DictReader(count_result.stdout.splitlines(), delimiter="\t"))
  aca_row = next(row for row in rows if row["Variation"] == "ACA>A")
  assert aca_row["Count"] == "1"

  decompose_cmd = [
    sys.executable,
    "-m",
    "mutational_signature.decompose",
    "--signatures",
    str(EXAMPLE / "mini_signatures.tsv"),
    "--counts",
    str(count_file),
    "--seed",
    "1",
  ]
  decompose_result = run_command(decompose_cmd)
  exposures_file.write_text(decompose_result.stdout)
  assert "SBS_Test_A" in decompose_result.stdout
  assert "Mutations\t1.0" in decompose_result.stdout

  plot_input = count_file.read_text()
  plot_cmd = [
    sys.executable,
    "-m",
    "mutational_signature.plot_counts",
    "--target",
    str(plot_file),
  ]
  run_command(plot_cmd, stdin_text=plot_input)
  assert plot_file.exists()
  assert plot_file.stat().st_size > 0


def test_console_script_help_smoke():
  help_cmd = [sys.executable, "-m", "mutational_signature.decompose", "--help"]
  result = run_command(help_cmd)
  assert "mutational signature finder" in result.stdout.lower()
