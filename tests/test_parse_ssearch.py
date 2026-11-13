"""Tests for parse_ssearch.py. Run from the repository root with:

    python -m unittest discover tests
"""

import unittest

from helpers import contig_of, parse_rows, read_fasta, read_report, run_script, ungapped

CONTIGS = read_fasta("contigs.fna")
MIRNAS = read_fasta("mirnas.fna")


def parse(report, *args):
    return parse_rows(run_script("parse_ssearch.py", args, stdin=read_report(report)))


class TestSubjectCoordinates(unittest.TestCase):
    """s_start/s_end must span exactly the subject region shown in subject_aln (issue #4)."""

    def check_coordinates(self, rows):
        checked = 0
        for row in rows:
            query, subject, query_aln, subject_aln = row[0], row[1], row[3], row[4]
            s_start, s_end = int(row[7]), int(row[8])
            contig = CONTIGS[contig_of(subject)]
            with self.subTest(row="\t".join(row)):
                self.assertTrue(1 <= s_start <= s_end <= len(contig))
                # With -W 10, ssearch36 truncates the displayed alignment when more than
                # 10 miRNA nt are unaligned; subject_aln is then incomplete. Skip those.
                if len(ungapped(query_aln)) < len(MIRNAS[query]):
                    continue
                expected = ungapped(subject_aln).replace("U", "T")
                self.assertEqual(contig[s_start - 1:s_end], expected)
                checked += 1
        self.assertGreater(checked, 100)

    def test_forward_workflow(self):
        self.check_coordinates(parse("ssearch_fwd.txt.gz"))

    def test_ssearch_i_workflow(self):
        # ssearch36 -i reports descending query coordinates, e.g. (18-1:9262-9279).
        self.check_coordinates(parse("ssearch_i.txt.gz"))

    def test_reverse_workflow(self):
        self.check_coordinates(parse("ssearch_rev.txt.gz", "--rev"))

    def test_planted_site(self):
        rows = parse("ssearch_fwd.txt.gz")
        coords = {(r[0], contig_of(r[1]), r[7], r[8]) for r in rows}
        self.assertIn(("ath-miR158a-3p_MIMAT0000176", "contig_1", "501", "520"), coords)

    def test_clamped_at_contig_start(self):
        # contig_2 starts with the last 12 nt of a miR165b site.
        rows = [r for r in parse("ssearch_fwd.txt.gz")
                if r[0] == "ath-miR165b_MIMAT0000188" and contig_of(r[1]) == "contig_2"
                and r[7] == "1"]
        self.assertTrue(rows)
        self.assertEqual(rows[0][8], "12")

    def test_clamped_at_contig_end(self):
        # contig_3 ends with the first 12 nt of a miR159c site.
        length = str(len(CONTIGS["contig_3"]))
        rows = [r for r in parse("ssearch_fwd.txt.gz")
                if r[0] == "ath-miR159c_MIMAT0001015" and contig_of(r[1]) == "contig_3"
                and r[8] == length]
        self.assertTrue(rows)
        self.assertEqual(int(rows[0][7]), len(CONTIGS["contig_3"]) - 11)


class TestManyQueries(unittest.TestCase):
    """From the 100th query on, ssearch36 headers have no leading whitespace (issue #2)."""

    def test_all_queries_reported_with_their_own_alignment(self):
        queries = read_fasta("mirnas_105_rc.fna")
        rows = parse("ssearch_105.txt.gz")
        self.assertEqual({r[0] for r in rows}, set(queries))
        for row in rows:
            with self.subTest(query=row[0]):
                self.assertIn(ungapped(row[3]), queries[row[0]])


if __name__ == "__main__":
    unittest.main()
