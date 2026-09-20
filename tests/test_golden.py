"""Golden-output regression: a fixed seed must reproduce exact bytes.

The hashes are real output, measured on this checkout. They only move when the
order or arithmetic of the random draws changes, which is exactly what this
test exists to catch.
"""

import hashlib
import sys

import pytest

import seq_sim_v2.seq_sim as S

TREE = "((A:0.2,B:0.2):0.3,(C:0.2,D:0.2):0.3);"

# (argv, output file relative to <prefix>/, expected sha256)
CASES = [
    (["-m", "HKY", "-l", "100", "-z", "12345", "-t", "4", "-o", "f", "-q"],
     "o_1.fasta",
     "a2cd4ce1165bf90b699d2c82d5053f3df350516df147cc7b3673af590e1840c6"),
    (["-m", "GTR", "-l", "100", "-z", "12345", "-r", "1,2,1,1,2,1",
      "-o", "p", "-q"],
     "o.phy",
     "5a7de3829cc47948eda098e7d54455d8d642a2a6f209f5bd930eb7c3da5908a7"),
    (["-m", "WAG", "-l", "80", "-z", "12345", "-o", "n", "-q"],
     "o.nex",
     "c067b9f1e8088f606ffb8f1c44f91ce6d1198a17171838f5ad74c15e3c7a47f7"),
]


@pytest.mark.parametrize("argv,filename,expected", CASES)
def test_golden_output(tmp_path, monkeypatch, argv, filename, expected):
    tree = tmp_path / "t.nwk"
    tree.write_text(TREE, encoding="utf-8")
    prefix = tmp_path / "o"
    monkeypatch.setattr(
        sys, "argv", ["seq_sim", *argv, "-y", str(prefix), str(tree)])
    S.main()
    data = (prefix / filename).read_bytes()
    assert hashlib.sha256(data).hexdigest() == expected
