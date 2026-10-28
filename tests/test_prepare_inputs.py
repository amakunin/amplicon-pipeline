"""Tests for ampliseq input preparation on the small fixture manifest."""

import gzip
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/small"
SCRIPT = ROOT / "scripts/prepare_inputs.py"
RECORD = b"@r\nACGT\n+\nIIII\n"


def run_prepare(demux_dir: Path, output_dir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--manifest", str(FIXTURE / "manifest.tsv"),
            "--config-dir", str(FIXTURE / "config"),
            "--demux-dir", str(demux_dir),
            "--target", "COI",
            "--output-dir", str(output_dir),
        ],
        capture_output=True,
        text=True,
        check=False,
    )


class PrepareInputsTests(unittest.TestCase):
    def test_samplesheet_excluded_samples_and_params(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            demux = Path(temporary) / "demux"
            demux.mkdir()
            # ANLCx_002_A12 has reads, BFHK_023_H8 is empty, the rest are absent
            with gzip.open(demux / "ont_barcode10_barcode89.fastq.gz", "wb") as handle:
                handle.write(RECORD)
            gzip.open(demux / "ont_barcode16_barcode64.fastq.gz", "wb").close()

            result = run_prepare(demux, Path(temporary) / "ampliseq")

            self.assertEqual(result.returncode, 0, result.stderr)
            target_dir = Path(temporary) / "ampliseq/COI"
            self.assertEqual(
                (target_dir / "samplesheet.tsv").read_text(),
                "sample\tfastq_1\n"
                f"ANLCx_002_A12\t{demux.resolve()}/ont_barcode10_barcode89.fastq.gz\n",
            )
            self.assertEqual(
                (target_dir / "excluded_samples.tsv").read_text(),
                "sample_id\treason\n"
                "BFHK_023_H8\tempty_fastq\n"
                "BFHK_023_A2\tno_demux_file\n"
                "ANLCx_001_A1\tno_demux_file\n"
                "ANLCx_001_B1\tno_demux_file\n",
            )
            parameters = yaml.safe_load((target_dir / "params.yaml").read_text())
            self.assertEqual(parameters["input"], str(target_dir / "samplesheet.tsv"))
            self.assertEqual(parameters["outdir"], str(target_dir / "results"))
            self.assertEqual(parameters["asv_calling"], "savont")


if __name__ == "__main__":
    unittest.main()
