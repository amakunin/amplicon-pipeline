"""Tests for ampliseq input preparation on the small fixture manifest."""

import gzip
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import yaml

from run_pipeline import load_config


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/small"
SCRIPT = ROOT / "scripts/prepare_inputs.py"
RECORD = b"@r\nACGT\n+\nIIII\n"


def run_prepare(preprocessed_dir: Path, output_dir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--manifest", str(FIXTURE / "manifest.tsv"),
            "--config-dir", str(FIXTURE / "config"),
            "--preprocessed-dir", str(preprocessed_dir),
            "--target", "CytB",
            "--output-dir", str(output_dir),
        ],
        capture_output=True,
        text=True,
        check=False,
    )


class PrepareInputsTests(unittest.TestCase):
    def test_bundled_coi_reference_in_generated_params(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "excluded_samples.tsv").write_text(
                "sample_id\treason\n"
                + "".join(
                    f"{sample}\tno_demux_file\n" for sample in (
                        "ANLCx_002_A12", "BFHK_023_H8", "BFHK_023_A2",
                        "ANLCx_001_A1", "ANLCx_001_B1",
                    )
                )
            )
            result = subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--manifest", str(FIXTURE / "manifest.tsv"),
                    "--config-dir", str(ROOT / "config/example"),
                    "--preprocessed-dir", str(directory),
                    "--target", "COI",
                    "--output-dir", str(directory / "ampliseq"),
                ],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            parameters = yaml.safe_load(
                (directory / "ampliseq/COI/params.yaml").read_text()
            )
            self.assertEqual(parameters["sintax_ref_taxonomy"], "coidb")
            self.assertNotIn("sintax_ref_tax_custom", parameters)
            self.assertNotIn("sintax_assign_taxlevels", parameters)

    def test_samplesheet_excluded_samples_and_params(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            preprocessed = Path(temporary) / "preprocessed"
            preprocessed.mkdir()
            with gzip.open(preprocessed / "barcode02_barcode89.fastq.gz", "wb") as handle:
                handle.write(RECORD)
            gzip.open(preprocessed / "barcode08_barcode64.fastq.gz", "wb").close()
            (preprocessed / "excluded_samples.tsv").write_text(
                "sample_id\treason\n"
                "BFHK_023_A2\tempty_fastq\n"
                "ANLCx_001_A1\tno_demux_file\n"
                "ANLCx_001_B1\tno_demux_file\n"
            )

            result = run_prepare(preprocessed, Path(temporary) / "ampliseq")

            self.assertEqual(result.returncode, 0, result.stderr)
            rerun = run_prepare(preprocessed, Path(temporary) / "ampliseq")
            self.assertEqual(rerun.returncode, 0, rerun.stderr)
            target_dir = Path(temporary) / "ampliseq/CytB"
            self.assertEqual(
                (target_dir / "samplesheet.tsv").read_text(),
                "sample\tfastq_1\n"
                f"ANLCx_002_A12\t{preprocessed.resolve()}/barcode02_barcode89.fastq.gz\n"
                f"BFHK_023_H8\t{preprocessed.resolve()}/barcode08_barcode64.fastq.gz\n",
            )
            self.assertEqual(
                (target_dir / "excluded_samples.tsv").read_text(),
                "sample_id\treason\n"
                "BFHK_023_A2\tempty_fastq\n"
                "ANLCx_001_A1\tno_demux_file\n"
                "ANLCx_001_B1\tno_demux_file\n",
            )
            parameters = yaml.safe_load((target_dir / "params.yaml").read_text())
            _, targets = load_config(FIXTURE / "config")
            self.assertEqual(parameters, {
                **targets["CytB"],
                "input": str(target_dir / "samplesheet.tsv"),
                "outdir": str(target_dir / "results"),
            })
            self.assertEqual(parameters["input"], str(target_dir / "samplesheet.tsv"))
            self.assertEqual(parameters["outdir"], str(target_dir / "results"))
            self.assertEqual(parameters["asv_calling"], "savont")
            self.assertEqual(
                parameters["sintax_ref_tax_custom"],
                str(
                    ROOT.parent / "culicidae_ref"
                    / "MIDORI2_UNIQ_NUC_GB272_Cytb_SINTAX.fasta"
                ),
            )
            before = (target_dir / "samplesheet.tsv").read_bytes()
            (preprocessed / "barcode02_barcode89.fastq.gz").unlink()
            missing = run_prepare(preprocessed, Path(temporary) / "ampliseq")
            self.assertNotEqual(missing.returncode, 0)
            self.assertIn("FileNotFoundError", missing.stderr)
            self.assertEqual((target_dir / "samplesheet.tsv").read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
