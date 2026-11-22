"""Newick and NEXUS parsing: flags, TRANSLATE scope, taxa, concatenation."""

import io

import pytest

import seq_sim.seq_sim as S
from seq_sim.io.nexus_reader import parse_nexus_tree_entries, parse_nexus_trees
from seq_sim.io.writers import write_ancestral_multi, write_sequences
from seq_sim.tree.tree_parser import allocate_sequences, parse_trees
from helpers import fake_tree as _fake_tree, run_cli as _run

def _tree(text):
    return parse_trees(text)[0][2]


# ---- NEXUS [&R] / [&U] must agree with the Newick structure

def test_rooted_flag_rejects_trifurcating_root():
    tree = _tree("(A:0.0,B:0.0,C:0.0);")
    assert tree.root.branch0 is not None
    with pytest.raises(SystemExit):
        S._check_rooted_flags([(0, 1.0, tree)], [True])


def test_unrooted_flag_rejects_bifurcating_root():
    tree = _tree("((A:0.0,B:0.0):0.0,C:0.0);")
    assert tree.root.branch0 is None
    with pytest.raises(SystemExit):
        S._check_rooted_flags([(0, 1.0, tree)], [False])


def test_matching_rooted_flags_are_applied():
    tree = _tree("((A:0.0,B:0.0):0.0,C:0.0);")
    S._check_rooted_flags([(0, 1.0, tree)], [True])
    assert tree.rooted is True


# ---- TRANSLATE is scoped to its TREES block

NEXUS_TWO_BLOCKS = (
    "#NEXUS\nBEGIN TREES;\nTRANSLATE 1 A, 2 B;\nTREE t1 = (1:0.0,2:0.0);\nEND;\n"
    "BEGIN TREES;\nTRANSLATE 1 W, 2 X;\nTREE t2 = (1:0.0,2:0.0);\nEND;\n")


def test_translate_is_block_scoped():
    entries = parse_nexus_tree_entries(NEXUS_TWO_BLOCKS)
    assert len(entries) == 2
    names = []
    for newick, _rooted, block_translate in entries:
        tree = parse_trees(newick, block_translate)[0][2]
        names.append([tree.names[i] for i in range(tree.num_tips)])
    assert names == [["A", "B"], ["W", "X"]]


def test_merged_translate_still_reported_for_compatibility():
    translate, trees = parse_nexus_trees(NEXUS_TWO_BLOCKS)
    assert translate.get(1) == "W" and translate.get(2) == "X"
    assert len(trees) == 2


# ---- a tree without a prefix after a pending partition is refused

def test_unprefixed_tree_after_pending_partition_rejected():
    raw = parse_trees("[60,1](A:0.0,B:0.0);\n(W:0.0,X:0.0);\n")
    with pytest.raises(SystemExit):
        S._group_partitions(raw, 100)


def test_zero_length_partition_rejected():
    raw = parse_trees("[0,1](A:0.0,B:0.0);\n")
    with pytest.raises(SystemExit):
        S._group_partitions(raw, 100)


# ---- partitions must carry exactly the same taxa

def test_extra_taxon_in_later_partition_rejected():
    t1 = _fake_tree(["A", "B"], [[0, 0], [1, 1]])
    t2 = _fake_tree(["A", "B", "C"], [[0, 0], [1, 1], [2, 2]])
    with pytest.raises(ValueError):
        write_sequences(io.StringIO(), "r", [t1, t2], [2, 2], 2, 4, "ACGT")


def test_duplicate_tip_names_rejected():
    t1 = _fake_tree(["A", "A"], [[0, 0], [1, 1]])
    with pytest.raises(ValueError):
        write_sequences(io.StringIO(), "r", [t1], [2], 2, 2, "ACGT")


# ---- ancestral concatenation compares topology, not labels

def test_ancestral_rejects_numeric_label_collision():
    trees = []
    for text in ("[3,1]((5:0.1,A:0.1):0.1,B:0.1);", "[3,1](5:0.1,(A:0.1,B:0.1):0.1);"):
        tree = parse_trees(text)[0][2]
        allocate_sequences(tree, 3)
        trees.append(tree)
    # Both trees render the same flat label list, so a label check accepts this.
    with pytest.raises(ValueError):
        write_ancestral_multi(io.StringIO(), "r", trees, "ACGT")


def test_ancestral_accepts_identical_topologies():
    trees = []
    for _ in range(2):
        tree = _tree("(A:0.1,B:0.1);")
        allocate_sequences(tree, 3)
        trees.append(tree)
    out = io.StringIO()
    write_ancestral_multi(out, "r", trees, "ACGT")
    rows = [ln.split() for ln in out.getvalue().splitlines()[1:] if ln.strip()]
    assert rows and all(len(r[1]) == 6 for r in rows)


# ---- each branch factor is validated, not just the product

def test_nan_branch_multiplier_rejected(tmp_path):
    tree_file = tmp_path / "t.nwk"
    tree_file.write_text("((A:0.1[nan],B:0.1):0.1,C:0.1);\n")
    proc = _run(["-m", "LG", "-l", "40", "-br", "-o", "p", str(tree_file)], tmp_path)
    assert proc.returncode != 0
    assert "not finite" in proc.stderr


def test_negative_length_with_negative_multiplier_rejected(tmp_path):
    tree_file = tmp_path / "t.nwk"
    tree_file.write_text("((A:-1[-1],B:0.1):0.1,C:0.1);\n")
    proc = _run(["-m", "LG", "-l", "40", "-br", "-o", "p", str(tree_file)], tmp_path)
    assert proc.returncode != 0
    assert "negative" in proc.stderr


NAMES = ["A", "'A[west]'", "'Taxon A'", "A_B", "1"]
ANNOTS = ["", "[2.0]", "[f=1,2.0]", "[f=1][2.0]", "[nan]", "[f=+1]", "[r=-2]"]


def _newick(name, annot):
    return "(%s:0.1%s,B:0.1,C:0.1);" % (name, annot)


def _sig(tree):
    # repr so that nan == nan (nan != nan would fail a tuple compare)
    return sorted((n.name, n.tip_no, repr(float(n.length0)),
                   repr(float(n.branch_rate)), n.comp_idx, n.rmat_idx)
                  for n in tree.node_list)


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("annot", ANNOTS)
def test_nexus_matches_newick(name, annot):
    newick = _newick(name, annot)
    nexus = "#NEXUS\nBEGIN TREES;\nTREE t1 = %s\nEND;\n" % newick
    (nx_nwk, __, ___), = parse_nexus_tree_entries(nexus)
    a = parse_trees(newick)[0][2]
    b = parse_trees(nx_nwk)[0][2]
    assert _sig(a) == _sig(b)


def test_translate_keeps_quoted_labels():
    nexus = ("#NEXUS\nBEGIN TREES;\n"
             "TRANSLATE 1 'A,B', 2 'Taxon B', 3 'C[d]';\n"
             "TREE t1 = (1:0.1,2:0.1,3:0.1);\nEND;\n")
    (nwk, __, tr), = parse_nexus_tree_entries(nexus)
    tips = sorted(n.name for n in parse_trees(nwk, tr)[0][2].node_list
                  if n.tip_no >= 0)
    assert tips == ["A,B", "C[d]", "Taxon B"]


def test_partition_prefix_survives():
    nexus = "#NEXUS\nBEGIN TREES;\nTREE t1 = [100,1](A:0.1,B:0.1);\nEND;\n"
    (nwk, __, ___), = parse_nexus_tree_entries(nexus)
    assert nwk.startswith("[100,1](")
