"""Tests for parse_mirna_targets.py (end to end, from ssearch36 reports). Run from the
repository root with:

    python -m unittest discover tests
"""

import os
import re
import tempfile
import unittest

from helpers import DATA, parse_rows, read_fasta, read_report, run_script


def pipeline(report, *args):
    """ssearch36 report | parse_ssearch.py | parse_mirna_targets.py, as in the README."""
    rev = ["--rev"] if "--rev" in args else []
    parsed = run_script("parse_ssearch.py", rev, stdin=read_report(report))
    return run_script("parse_mirna_targets.py", args, stdin=parsed)


def expected(name):
    with open(os.path.join(DATA, name)) as fh:
        return fh.read()


class TestRegression(unittest.TestCase):
    """Full pipeline output must match the stored reference output.

    If a change intentionally alters the results, regenerate the references with
    tests/make_fixtures.py and review the diff before committing.
    """

    def test_forward_workflow(self):
        self.assertEqual(pipeline("ssearch_fwd.txt.gz"), expected("expected_targets_fwd.tsv"))

    def test_reverse_workflow(self):
        self.assertEqual(pipeline("ssearch_rev.txt.gz", "--rev"),
                         expected("expected_targets_rev.tsv"))

    def test_ssearch_i_gives_same_targets_as_forward(self):
        # Using ssearch36 -i instead of rev-complementing the miRNAs beforehand must give
        # the same targets and scores. q_start/q_end are numbered on the original instead
        # of the rev-complemented miRNA, so ignore those.
        def comparable(tsv):
            return [r[:5] + r[7:] for r in parse_rows(tsv)]

        self.assertEqual(comparable(pipeline("ssearch_i.txt.gz")),
                         comparable(expected("expected_targets_fwd.tsv")))


class TestScores(unittest.TestCase):
    """Penalty scores of the planted target sites."""

    def score(self, rows, query, s_start):
        matches = [float(r[9]) for r in rows if r[0] == query and r[7] == s_start]
        self.assertEqual(len(matches), 1)
        return matches[0]

    def test_planted_sites(self):
        fwd = parse_rows(pipeline("ssearch_fwd.txt.gz"))
        # Perfect site.
        self.assertEqual(self.score(fwd, "ath-miR158a-3p_MIMAT0000176", "501"), 0.0)
        # 2 mismatches at miRNA positions 20-21, outside the 19 scored positions.
        self.assertEqual(self.score(fwd, "ath-miR159a_MIMAT0000177", "1221"), 0.0)
        # 1 mismatch outside the seed region.
        self.assertEqual(self.score(fwd, "ath-miR161.1_MIMAT0000181", "2042"), 1.0)

        rev = parse_rows(pipeline("ssearch_rev.txt.gz", "--rev"))
        # Perfect minus-strand site.
        self.assertEqual(self.score(rev, "ath-miR165b_MIMAT0000188", "2663"), 0.0)
        # 1 mismatch in the seed region: 1 x penalty multiplier (1.5).
        self.assertEqual(self.score(rev, "ath-miR158a-3p_MIMAT0000176", "3284"), 1.5)


class TestECutoff(unittest.TestCase):
    """Pairs with an expectation greater than --E_cutoff are discarded (issue #3)."""

    def test_cutoff_is_an_upper_bound(self):
        counts = {}
        for cutoff in (5.0, 4.0, 1.0):
            rows = parse_rows(pipeline("ssearch_fwd.txt.gz", "--E_cutoff", str(cutoff)))
            with self.subTest(cutoff=cutoff):
                self.assertTrue(rows)
                self.assertTrue(all(float(r[9]) <= cutoff for r in rows))
            counts[cutoff] = len(rows)
        self.assertGreater(counts[5.0], counts[4.0])
        self.assertGreater(counts[4.0], counts[1.0])

    def test_cutoff_is_inclusive(self):
        rows = parse_rows(pipeline("ssearch_fwd.txt.gz", "--E_cutoff", "1.0"))
        self.assertIn("1.0", [r[9] for r in rows])


class TestOptions(unittest.TestCase):

    def test_target_suffix_removed_by_default(self):
        contigs = set(read_fasta("contigs.fna"))
        rows = parse_rows(pipeline("ssearch_fwd.txt.gz"))
        self.assertTrue(rows)
        self.assertTrue(all(r[1] in contigs for r in rows))

    def test_keep_target_suffix(self):
        default = parse_rows(pipeline("ssearch_fwd.txt.gz"))
        kept = parse_rows(pipeline("ssearch_fwd.txt.gz", "--keep_target_suffix"))
        self.assertEqual(len(kept), len(default))
        for k, d in zip(kept, default):
            self.assertRegex(k[1], r"^contig_\d+_\d+$")
            self.assertEqual([k[0], re.sub(r"_\d+$", "", k[1])] + k[2:], d)

    def test_outfile_failed(self):
        with tempfile.TemporaryDirectory() as tmp:
            failed_path = os.path.join(tmp, "failed.tsv")
            passed = parse_rows(pipeline("ssearch_fwd.txt.gz", "-o", failed_path))
            with open(failed_path) as fh:
                failed = parse_rows(fh.read())
        self.assertEqual(passed, parse_rows(expected("expected_targets_fwd.tsv")))
        self.assertTrue(failed)
        self.assertTrue(all(len(r) == 11 for r in failed))
        self.assertTrue(any(float(r[9]) > 5.0 for r in failed))


if __name__ == "__main__":
    unittest.main()
