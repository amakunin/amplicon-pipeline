"""Hierarchical read accounting with small synthetic observations."""

import csv
import gzip
from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts.count_reads import count_reads


class CountReadsTests(unittest.TestCase):
    def test_hierarchy_counts_missing_values_and_qc(self) -> None:
        for method, prefix in (("illumina", "barcode"), ("hybrid", "ont_barcode")):
            with self.subTest(method=method), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                raw = root / "run.fastq.gz"
                with gzip.open(raw, "wb") as handle:
                    handle.write(b"@r\nACGT\n+\nIIII\n" * 14)
                samples = [
                    {"sample_id": name, "barcode_indices": {"first": first, "rear": rear}}
                    for name, first, rear in (("S", "01", "02"), ("T", "02", "03"), ("M", "03", "04"))
                ]
                config = {"demux": {"method": method}}

                def write(relative, text):
                    path = root / relative
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(text)

                write("demux/demux_stats.txt",
                      f"{prefix}01_barcode02\t10\n{prefix}01_barcode03\t1\n"
                      f"{prefix}01_unclassified\t1\n{prefix}09_barcode02\t1\nunclassified\t1\n")
                write("preprocessed/read_counts.tsv",
                      "sample_id\tinput_reads\tporechop_reads\tchopper_reads\n"
                      "S\t10\t9\t8\nT\t0\t\t\nM\t\t\t\n")
                write("summary/asvs.tsv",
                      "target\tasv_id\tlength_status\tcoding_status\ttotal_reads\n"
                      "A\tx\tpass\tone_stop_free_frame\t4\n"
                      "A\ty\ttoo_short\tnot_checked\t1\n"
                      "A\tz\tnot_checked\tuncertain\t1\n")
                write("summary/assignments.tsv",
                      "sample_id\ttarget\tasv_id\tread_count\tcontam_possible\tlocov\n"
                      "S\tA\tx\t3\tfalse\tfalse\nS\tA\ty\t1\tfalse\tfalse\n"
                      "S\tA\tz\t1\tfalse\tfalse\nT\tA\tx\t1\ttrue\ttrue\n")
                write("ampliseq/A/results/overall_summary.tsv",
                      "sample\tcutadapt_total_processed\tcutadapt_passing_filters\tsavont_output\n"
                      "S\t8\t6\t5\nT\t0\t0\tNA\n")
                write("ampliseq/A/results/savont/pooled_feature-table.tsv",
                      "ASV_ID\tS\tT\nx\t3\t1\ny\t1\t0\nz\t1\t0\n")
                warnings = StringIO()
                with (
                    patch("scripts.count_reads.load_config", return_value=(config, {"A": {}})),
                    patch("scripts.count_reads.read_manifest", return_value=samples),
                    redirect_stderr(warnings),
                ):
                    count_reads(raw, root / "manifest.tsv", root, root)
                self.assertEqual(warnings.getvalue(), "")
                with (root / "summary/read_flow.tsv").open() as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    self.assertEqual(reader.fieldnames, [
                        "node_id", "parent_id", "stage", "first_barcode",
                        "second_barcode", "sample_id", "target", "label", "read_count",
                    ])
                    rows = list(reader)
                ids = set()
                for row in rows:
                    if row["parent_id"]:
                        self.assertIn(row["parent_id"], ids)
                    ids.add(row["node_id"])
                self.assertEqual(rows[0]["read_count"], "14")
                self.assertTrue(any(r["first_barcode"] == prefix + "09" and r["label"] == "unexpected" for r in rows))
                self.assertTrue(any(r["second_barcode"] == "barcode03" and r["first_barcode"] == prefix + "01" and r["label"] == "unexpected" for r in rows))
                self.assertEqual(sum(r["label"] == "unmatched" for r in rows), 2)
                observed = {(r["sample_id"], r["stage"]): r["read_count"] for r in rows if r["sample_id"]}
                self.assertEqual(observed[("S", "sequence_qc_pass")], "3")
                self.assertEqual(observed[("S", "both_qc_pass")], "3")
                self.assertEqual(observed[("T", "sequence_qc_pass")], "1")
                self.assertEqual(observed[("T", "both_qc_pass")], "0")
                self.assertEqual(observed[("T", "savont_output")], "NA")
                self.assertEqual(observed[("M", "both_qc_pass")], "")
                self.assertEqual(observed[("T", "porechop")], "")
                self.assertEqual(observed[("T", "sample")], "0")
                write("preprocessed/read_counts.tsv",
                      "sample_id\tinput_reads\tporechop_reads\tchopper_reads\n"
                      "S\t11\t12\t8\nT\t0\t\t\nM\t\t\t\n")
                with gzip.open(raw, "wb") as handle:
                    handle.write(b"@r\nACGT\n+\nIIII\n" * 15)
                warnings = StringIO()
                with (
                    patch("scripts.count_reads.load_config", return_value=(config, {"A": {}})),
                    patch("scripts.count_reads.read_manifest", return_value=samples),
                    redirect_stderr(warnings),
                ):
                    count_reads(raw, root / "manifest.tsv", root, root)
                self.assertIn("S demux/input: 10 -> 11", warnings.getvalue())
                self.assertIn("S porechop: 11 -> 12", warnings.getvalue())
                self.assertIn("raw/demux total: 15 -> 14", warnings.getvalue())
                write("ampliseq/B/results/overall_summary.tsv",
                      "sample\tcutadapt_total_processed\tcutadapt_passing_filters\tsavont_output\n"
                      "S\t7\t6\t4\n")
                write("ampliseq/B/results/savont/pooled_feature-table.tsv", "ASV_ID\tS\n")
                warnings = StringIO()
                with (
                    patch("scripts.count_reads.load_config", return_value=(config, {"A": {}, "B": {}})),
                    patch("scripts.count_reads.read_manifest", return_value=samples),
                    redirect_stderr(warnings),
                ):
                    count_reads(raw, root / "manifest.tsv", root, root)
                self.assertIn("primer inputs A/B: 8 -> 7", warnings.getvalue())
                self.assertIn("S/B Savont/features: 4 -> 0", warnings.getvalue())
                before = (root / "summary/read_flow.tsv").read_bytes()
                path = root / "summary/assignments.tsv"
                path.write_text(path.read_text().replace("false", "invalid", 1))
                with (
                    patch("scripts.count_reads.load_config", return_value=(config, {"A": {}})),
                    patch("scripts.count_reads.read_manifest", return_value=samples),
                ):
                    with self.assertRaisesRegex(ValueError, "Invalid assignment QC flag"):
                        count_reads(raw, root / "manifest.tsv", root, root)
                self.assertEqual((root / "summary/read_flow.tsv").read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
