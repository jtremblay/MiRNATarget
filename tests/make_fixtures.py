#!/usr/bin/env python
"""Regenerate the test fixtures in tests/data/.

Only needed if the fixtures must be rebuilt; the tests themselves do not need
ssearch36. Requires ssearch36 (FASTA36 package) in PATH. Run from anywhere:

    python tests/make_fixtures.py

Creates small synthetic contigs (fixed random seed) with planted miRNA target
sites, including partial sites at a contig start and end to exercise
coordinate clamping, then runs ssearch36 with the parameters from the README
for the three documented workflows:
  fwd: rev-complemented miRNAs                 -> parse without --rev
  i:   original miRNAs with ssearch36 -i       -> parse without --rev (same targets as fwd)
  rev: original miRNAs (minus-strand targets)  -> parse with --rev
"""

import gzip
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
REPO = os.path.dirname(HERE)

SSEARCH_ARGS = ["-f", "-8", "-g", "-3", "-E", "10000", "-T", "1", "-r", "+4/-3",
                "-n", "-U", "-W", "10", "-N", "20000"]

COMP = {"A": "U", "C": "G", "G": "C", "U": "A"}


def read_fasta(path):
    recs, name = [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
            elif line:
                recs.append((name, line))
    return recs


def revcomp_rna(seq):
    return "".join(COMP[c] for c in reversed(seq))


def write_fasta(path, recs, width=60):
    with open(path, "w") as fh:
        for name, seq in recs:
            fh.write(">" + name + "\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def mutate(seq, positions, rng):
    s = list(seq)
    for p in positions:
        s[p] = rng.choice([c for c in "ACGT" if c != s[p]])
    return "".join(s)


def main():
    os.makedirs(DATA, exist_ok=True)
    rng = random.Random(42)
    mirnas = read_fasta(os.path.join(REPO, "data", "test.fna"))
    # Target site (DNA, + strand of the contig) = reverse complement of the miRNA.
    sites = [revcomp_rna(s).replace("U", "T") for _, s in mirnas]

    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    # Minus-strand site (found by the --rev workflow) = the miRNA sequence itself, as DNA.
    minus = [s.replace("U", "T") for _, s in mirnas]

    contig1 = (rand(500) + sites[0]                      # perfect site
               + rand(700) + mutate(sites[1], [0, 1], rng)  # 2 mismatches at the site's 5' end (miRNA 3' end)
               + rand(800) + mutate(sites[2], [3], rng)     # 1 mismatch
               + rand(600) + minus[3]                       # perfect minus-strand site
               + rand(600) + mutate(minus[0], [5], rng)     # minus-strand site, 1 mismatch in the seed
               + rand(600))
    contig2 = sites[3][-12:] + rand(2000)    # partial site at the very start
    contig3 = rand(2000) + sites[5][:12]     # partial site at the very end
    write_fasta(os.path.join(DATA, "contigs.fna"),
                [("contig_1", contig1), ("contig_2", contig2), ("contig_3", contig3)])

    write_fasta(os.path.join(DATA, "mirnas.fna"), mirnas)
    write_fasta(os.path.join(DATA, "mirnas_rc.fna"), [(n, revcomp_rna(s)) for n, s in mirnas])
    # 105 queries (> 99) to exercise query headers without leading whitespace (issue #2).
    many = [("q%d_%s" % (k + 1, mirnas[k % len(mirnas)][0]), revcomp_rna(mirnas[k % len(mirnas)][1]))
            for k in range(105)]
    write_fasta(os.path.join(DATA, "mirnas_105_rc.fna"), many)

    def run(extra, query, out):
        # Reports are stored gzipped to keep the repository small.
        cmd = ["ssearch36"] + extra + [query, "contigs.fna"]
        report = subprocess.run(cmd, cwd=DATA, capture_output=True, check=True).stdout
        with gzip.open(os.path.join(DATA, out), "wb") as fh:
            fh.write(report)
        return report.decode()

    # -b/-d/-H and the 2nd -E value only limit the size of the report
    # (fewer alignments, no histogram); they do not change the alignments.
    small = ["-b", "30", "-d", "30", "-H"]
    reports = {
        "fwd": run(SSEARCH_ARGS + small, "mirnas_rc.fna", "ssearch_fwd.txt.gz"),
        "i": run(SSEARCH_ARGS + small + ["-i"], "mirnas.fna", "ssearch_i.txt.gz"),
        "rev": run(SSEARCH_ARGS + small, "mirnas.fna", "ssearch_rev.txt.gz"),
    }
    args_105 = [a if a != "10000" else "10000 0.0001" for a in SSEARCH_ARGS]
    run(args_105 + ["-b", "1", "-d", "1", "-H"], "mirnas_105_rc.fna", "ssearch_105.txt.gz")

    # Expected end-to-end output (regression reference).
    for direction, extra in (("fwd", []), ("rev", ["--rev"])):
        parsed = subprocess.run([sys.executable, os.path.join(REPO, "parse_ssearch.py")] + extra,
                                input=reports[direction], capture_output=True, text=True,
                                check=True).stdout
        targets = subprocess.run([sys.executable, os.path.join(REPO, "parse_mirna_targets.py")] + extra,
                                 input=parsed, capture_output=True, text=True, check=True).stdout
        with open(os.path.join(DATA, "expected_targets_%s.tsv" % direction), "w") as fh:
            fh.write(targets)

if __name__ == "__main__":
    main()
