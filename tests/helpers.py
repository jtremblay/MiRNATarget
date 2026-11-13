"""Shared helpers for the tests. Standard library only."""

import gzip
import os
import subprocess
import sys

TESTS = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(TESTS, "data")
REPO = os.path.dirname(TESTS)

COLUMNS = ["query", "subject", "match_aln", "query_aln", "subject_aln",
           "q_start", "q_end", "s_start", "s_end"]


def read_fasta(name):
    seqs, current = {}, None
    with open(os.path.join(DATA, name)) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif line:
                seqs[current].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def read_report(name):
    """Text of a gzipped ssearch36 report from tests/data."""
    with gzip.open(os.path.join(DATA, name), "rt") as fh:
        return fh.read()


def run_script(script, args=(), stdin=""):
    """Run one of the repository scripts like a user would; return stdout."""
    result = subprocess.run([sys.executable, os.path.join(REPO, script)] + list(args),
                            input=stdin, capture_output=True, text=True)
    if result.returncode != 0:
        raise AssertionError("%s %s failed:\n%s" % (script, " ".join(args), result.stderr))
    return result.stdout


def parse_rows(tsv):
    """Rows of a parse_ssearch.py / parse_mirna_targets.py output as lists of fields."""
    return [line.split("\t") for line in tsv.splitlines()
            if line and not line.startswith("#")]


def contig_of(subject):
    """parse_ssearch.py appends _<hit number> to subject ids."""
    return subject.rsplit("_", 1)[0]


def ungapped(aln):
    return aln.replace("-", "").replace(" ", "")
