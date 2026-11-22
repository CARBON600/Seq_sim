"""Shared test scaffolding: repository paths, a fake tree, and a CLI runner."""
import os
import subprocess
import sys
import types

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The CLI writes into the directory it runs in, so runs happen one level above
# the checkout and never litter the tree.
CWD = os.path.dirname(REPO)


def fake_tree(names, seqs):
    """Tree with list-shaped tips, as writers.write_sequences expects."""
    tree = types.SimpleNamespace()
    tree.num_tips = len(names)
    tree.names = {i: n for i, n in enumerate(names)}
    tree.tips = [types.SimpleNamespace(sequence=np.array(s, dtype=np.int64))
                 for s in seqs]
    return tree


def fake_tree_dict(names, seqs):
    """Tree with dict-shaped tips, as writers.write_sequences expects."""
    tree = types.SimpleNamespace()
    tree.num_tips = len(names)
    tree.names = {i: n for i, n in enumerate(names)}
    tree.tips = {i: types.SimpleNamespace(sequence=np.array(s, dtype=np.int64))
                 for i, s in enumerate(seqs)}
    return tree


def run_cli(args, tmp_path):
    """Run python -m seq_sim from CWD, writing under tmp_path/out."""
    return subprocess.run(
        [sys.executable, "-m", "seq_sim"] + args + ["-y", str(tmp_path / "out")],
        cwd=CWD, capture_output=True, text=True)
