"""The path from command line to file: data plumbing, input guards, parity."""

import io

import numpy as np
import pytest

import seq_sim.seq_sim as S
from seq_sim.core.aa_models import create_aa_model
from seq_sim.core.evolve import _random_sequence_profiles, evolve_sequences
from seq_sim.core.nuc_models import F84Model, GTRModel, HKYModel, JC69Model
from seq_sim.core.profile_mixture import assign_site_profiles, normalise_weights
from seq_sim.core.rate_heterogeneity import normalise_branch_rates
from seq_sim.io.nexus_reader import _strip_bracket_content, parse_nexus_trees
from seq_sim.io.writers import write_ancestral_multi, write_sequences
from seq_sim.tree.tree_parser import allocate_sequences, parse_trees
from helpers import fake_tree as _fake_tree, fake_tree_dict, run_cli as _run

# ---- tip rows are matched by name, not by leaf order

def test_tip_rows_follow_names_not_leaf_order():
    tree1 = _fake_tree(["A", "B"], [[0, 0, 0, 0], [1, 1, 1, 1]])   # AAAA / CCCC
    tree2 = _fake_tree(["B", "A"], [[2, 2], [3, 3]])               # B=GG / A=TT
    out = io.StringIO()
    write_sequences(out, "r", [tree1, tree2], [4, 2], 2, 6, "ACGT")
    rows = {line.split()[0]: line.split()[1]
            for line in out.getvalue().splitlines()[1:] if line.strip()}
    assert rows["A"] == "AAAATT"
    assert rows["B"] == "CCCCGG"


def test_tip_rows_refuse_missing_taxon():
    tree1 = _fake_tree(["A", "B"], [[0, 0], [1, 1]])
    tree2 = _fake_tree(["A", "C"], [[2, 2], [3, 3]])
    with pytest.raises(ValueError):
        write_sequences(io.StringIO(), "r", [tree1, tree2], [2, 2], 2, 4, "ACGT")


# ---- -r is honoured under GFmix

def test_custom_rates_honoured_under_gfmix():
    rates = [1.0 + i for i in range(190)]
    assert S._exchangeability_for("LG", rates) == rates
    assert S._exchangeability_for("LG", None) != rates


def test_general_still_requires_rates():
    with pytest.raises(SystemExit):
        S._exchangeability_for("GENERAL", None)


# ---- root profile slice and root composition

class _FixedModel:
    """A model whose root draw is deterministic (point mass at one state)."""

    def __init__(self, state):
        self.num_states = 4
        self.cum_freqs = np.zeros(4)
        self.cum_freqs[state] = 1.0

    def set_vector(self, vec, state, length):
        vec[:] = 0.0
        vec[state] = 1.0
        np.cumsum(vec, out=vec)


def _rate_info_none():
    return {"type": "none", "invariable": False, "invariable_prop": 0.0}


def test_root_profile_uses_root_models():
    tree = parse_trees("(A:0.0,B:0.0);")[0][2]
    allocate_sequences(tree, 4)
    base = [_FixedModel(0), _FixedModel(0)]        # base profiles -> A
    root = [_FixedModel(3), _FixedModel(3)]        # root edge gamma -> T
    profile_info = {"base_models": base, "per_node": {tree.root: root},
                    "site_profile": np.zeros(4, dtype=np.int64)}
    evolve_sequences(tree, 0, 4, 1.0, {}, _rate_info_none(), None, profile_info)
    for i in range(tree.num_tips):
        assert list(tree.tips[i].sequence) == [3, 3, 3, 3]


def test_root_profile_honours_from_site():
    seq = [-1, -1]
    _random_sequence_profiles(seq, 1, 1, [_FixedModel(0), _FixedModel(3)], [0, 1])
    assert seq[0] == 3


# ---- unrooted normalisation

def test_normalise_branch_rates_unrooted_mean_is_one():
    _, _, tree = parse_trees("(A:1[1],B:2[1],C:5[4]);")[0]
    normalise_branch_rates([[(-1, 1.0, tree)]])
    edges = [(n.length0, n.branch_rate) for n in tree.node_list if n is not tree.root]
    weight = sum(length for length, _ in edges)
    assert sum(length * branch_rate for length, branch_rate in edges) / weight == pytest.approx(1.0)


# ---- numeric NEXUS annotations survive

def test_nexus_keeps_numeric_branch_rate():
    text = _strip_bracket_content("(A:0.1[2.0],B:0.1);", keep_fr=True)
    assert "[2.0]" in text
    assert "[&R]" not in _strip_bracket_content("[&R] (A:0.1,B:0.1);", keep_fr=True)
    _, trees = parse_nexus_trees(
        "#NEXUS\nBEGIN TREES;\nTREE t1 = [&R] (A:0.1[2.0],B:0.1);\nEND;\n")
    newick = trees[0][0] if isinstance(trees[0], tuple) else trees[0]
    assert "[2.0]" in newick


# ---- partition sum must match -l

def test_partition_sum_must_match_length():
    raw = parse_trees("[60,1](A:0.1,B:0.1);")
    with pytest.raises(SystemExit):
        S._group_partitions(raw, 100)
    assert S._group_partitions(raw, 60)


# ---- negative branch lengths

@pytest.mark.parametrize("model", [
    JC69Model(),
    HKYModel([0.25] * 4, 3.0, False),
    F84Model([0.25] * 4, 3.0, False),
    GTRModel([0.25] * 4, [1, 2, 3, 4, 5, 6]),
    create_aa_model("LG"),
])
def test_negative_branch_length_rejected(model):
    vec = np.zeros(model.num_states)
    with pytest.raises(ValueError):
        model.set_vector(vec, 0, -0.5)


# ---- ancestral output covers every partition

def test_ancestral_multi_partition_concatenates():
    trees = []
    for length in (3, 4):
        tree = parse_trees("(A:0.1,B:0.1);")[0][2]
        allocate_sequences(tree, length)
        trees.append(tree)
    out = io.StringIO()
    write_ancestral_multi(out, "r", trees, "ACGT")
    rows = [line.split() for line in out.getvalue().splitlines()[1:] if line.strip()]
    assert rows and all(len(row[1]) == 7 for row in rows)


# ---- codon rates and CLI level

def test_two_trees_two_replicates_do_not_overwrite(tmp_path):
    tree_file = tmp_path / "t.nwk"
    tree_file.write_text("(A:0.1,B:0.1);\n(W:0.1,X:0.1);\n")
    proc = _run(["-m", "JC69", "-l", "20", "-n", "2", "-o", "f",
                 str(tree_file)], tmp_path)
    assert proc.returncode == 0, proc.stderr
    files = sorted((tmp_path / "out").glob("*.fasta"))
    assert len(files) == 4
    heads = [f.read_text().splitlines()[0] for f in files]
    assert ">A" in heads and ">W" in heads


def test_gfmix_negative_weight_rejected(tmp_path):
    tree_file = tmp_path / "t.nwk"
    tree_file.write_text("(A:0.1,B:0.1);\n")
    weights = "0.55,-0.05,0.1,0.1,0.1,0.05,0.05,0.05,0.02,0.03"
    proc = _run(["-m", "LG", "--gfmix-cseries", "C10", "-l", "40",
                 "--gfmix-weights", weights, str(tree_file)], tmp_path)
    assert proc.returncode != 0
    assert "non-negative" in proc.stderr


def test_negative_codon_rate_rejected(tmp_path):
    tree_file = tmp_path / "t.nwk"
    tree_file.write_text("(A:0.5,B:0.5);\n")
    proc = _run(["-m", "JC69", "-l", "30", "-c", "-1", "2", "2",
                 str(tree_file)], tmp_path)
    assert proc.returncode != 0
    assert "non-negative" in proc.stderr


def _parse(argv):
    old = S.sys.argv
    S.sys.argv = ["seq_sim"] + argv
    try:
        return S._parse_args()
    finally:
        S.sys.argv = old


# partition rates are validated before the weighted mean is taken
def test_negative_partition_rate_rejected_before_rescale():
    parts = [(1, -1.0, None), (1, -3.0, None)]
    with pytest.raises(SystemExit):
        S._rescale_partitions(parts, 2)


def test_zero_sum_partition_rates_rejected():
    parts = [(1, 0.0, None), (1, 0.0, None)]
    with pytest.raises(SystemExit):
        S._rescale_partitions(parts, 2)


# full format names are not truncated into "nexu"
def test_full_format_name_nexus_accepted():
    args = _parse(["-m", "JC69", "-o", "nexus"])
    assert args.format == "n"
    assert args.separate is False


def test_trailing_s_still_means_separate_files():
    args = _parse(["-m", "JC69", "-o", "ps"])
    assert args.format == "p"
    assert args.separate is True


# strict PHYLIP must not print two taxa under the same 10-char label
def test_phylip_truncation_collision_rejected():
    tree = fake_tree_dict(["ABCDEFGHIJ1", "ABCDEFGHIJ2"],
                     [np.zeros(4, dtype=np.int64), np.zeros(4, dtype=np.int64)])
    f = io.StringIO()
    with pytest.raises(ValueError):
        write_sequences(f, "p", [tree], [4], 2, 4, "ACGT")


# ancestral labels must stay unique when a tip has a numeric name
def test_ancestral_label_collision_rejected():
    t = parse_trees("(3:0.1,4:0.1);")[0][2]
    allocate_sequences(t, 4)
    f = io.StringIO()
    with pytest.raises(ValueError):
        write_ancestral_multi(f, "p", [t], "ACGT")


# -t on a model without ts/tv is reported, not ignored
def test_tstv_warns_on_gtr(capsys):
    args = _parse(["-m", "GTR", "-t", "4"])
    S._warn_ignored_model_args(args)
    err = capsys.readouterr().err
    assert "has no ts/tv parameter" in err


# non-finite scales are rejected
@pytest.mark.parametrize("flag", [["-s", "nan"], ["-s", "inf"], ["-d", "nan"]])
def test_nonfinite_scale_rejected(flag):
    with pytest.raises(SystemExit):
        _parse(["-m", "JC69"] + flag)


def test_nonfinite_gamma_rejected():
    with pytest.raises(SystemExit):
        _parse(["-m", "JC69", "-a", "nan"])


# accepted weights are renormalised so the cumulative sum ends at 1
def test_weights_renormalised_to_exact_one():
    w = normalise_weights([0.9999999995] + [0.0] * 9)
    assert float(np.cumsum(w)[-1]) == 1.0
    assert w[-1] == 0.0
    # with an exact sum, a draw strictly below 1 can never reach index 10
    idx = int(np.searchsorted(np.cumsum(w), 0.9999999997, side="right"))
    assert idx < 10


# negative control for the normaliser: plain searchsorted on the raw table
# returns 10, one past the last class; the normalised table returns below 10
def test_unnormalised_weights_can_reach_the_last_class():
    raw = np.array([0.9999999995] + [0.0] * 9)
    idx = int(np.searchsorted(np.cumsum(raw), 0.9999999997, side="right"))
    assert idx == 10


# the AA matrix path validates its own length, like set_vector
def test_aa_set_matrix_rejects_negative_length():
    m = create_aa_model("LG")
    buf = np.zeros((20, 20))
    with pytest.raises(ValueError):
        m.set_matrix(buf, -0.5)


# a mixed annotation must survive NEXUS exactly as it does in Newick
def test_nexus_keeps_mixed_model_and_rate_annotation():
    nx = "#NEXUS\nBEGIN TREES;\nTREE t = (A:0.1[f=1,2.0],B:0.1);\nEND;\n"
    _, trees = parse_nexus_trees(nx)
    assert "[f=1,2.0]" in trees[0][0]
    tree = parse_trees(trees[0][0])[0][2]
    node = tree.node_list[1]
    assert node.comp_idx == 1
    assert node.branch_rate == 2.0


def test_nexus_keeps_partition_prefix():
    nx = "#NEXUS\nBEGIN TREES;\nTREE t = [100,1](A:0.1,B:0.1);\nEND;\n"
    _, trees = parse_nexus_trees(nx)
    assert trees[0][0].startswith("[100,1](")


# ancestral strict PHYLIP shares the 10-character label with the tip writer
def test_ancestral_strict_phylip_truncates_names():
    t = parse_trees("(ABCDEFGHIJK:0.05,BBBBBBBBBBK:0.05);")[0][2]
    allocate_sequences(t, 4)
    f = io.StringIO()
    write_ancestral_multi(f, "p", [t], "ACGT")
    rows = f.getvalue().splitlines()[1:]
    assert rows and all(len(r) == 14 for r in rows), rows   # 10 + 4 sites


def test_ancestral_strict_phylip_rejects_collision():
    t = parse_trees("(ABCDEFGHIJ1:0.05,ABCDEFGHIJ2:0.05);")[0][2]
    allocate_sequences(t, 4)
    f = io.StringIO()
    with pytest.raises(ValueError):
        write_ancestral_multi(f, "p", [t], "ACGT")


# JC69 cannot express alternate compositions
def test_jc69_rejects_alt_freqs():
    with pytest.raises(SystemExit):
        _parse(["-m", "JC69", "--alt-freqs", "0.4,0.3,0.2,0.1"])


# one name must not address two edges of the same tree
def test_gfmix_duplicate_branch_name_rejected():
    tree = parse_trees("((A:0.1,B:0.1)H:0.2,(C:0.1,D:0.1)H:0.2);")[0][2]
    with pytest.raises(SystemExit):
        S._apply_branch_gammas([[(2, 1.0, tree)]], {"H": (2.0, 0.5)})


# an overflowing weighted length must not zero the multipliers
def test_normalise_branch_rates_rejects_overflow():
    tree = parse_trees("(A:1[1e308],B:1[1e308]);")[0][2]
    with pytest.raises(ValueError):
        normalise_branch_rates([[(2, 1.0, tree)]])
    assert all(n.branch_rate != 0.0 for n in tree.node_list)


# a zero-weight profile can never be drawn
def test_zero_weight_profile_never_drawn():
    w = np.array([0.5, 0.4999999995, 0.0])
    got = assign_site_profiles(50, w, lambda: 0.9999999997)
    assert set(int(x) for x in got) == {1}


def test_positive_weights_still_span_every_class():
    rng = np.random.default_rng(0).random
    got = assign_site_profiles(20000, np.array([0.1, 0.2, 0.3, 0.4]), rng)
    freq = np.bincount(got, minlength=4) / 20000.0
    assert np.all(np.abs(freq - np.array([0.1, 0.2, 0.3, 0.4])) < 0.02)


# ---- bad input in the tree file is reported, not raised at the user

def test_bad_partition_rate_reports_instead_of_raising(tmp_path):
    tree_file = tmp_path / "t.nwk"
    tree_file.write_text("[100,-1.0] ((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);\n")
    proc = _run(["-m", "JC69", "-l", "100", str(tree_file)], tmp_path)
    assert proc.returncode != 0
    assert "partition rate must be finite and positive" in proc.stderr
    assert "Traceback" not in proc.stderr


def test_impossible_branch_scale_reports_instead_of_raising(tmp_path):
    tree_file = tmp_path / "t.nwk"
    tree_file.write_text("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);\n")
    proc = _run(["-m", "JC69", "-l", "40", "-d", "1e308", str(tree_file)], tmp_path)
    assert proc.returncode != 0
    assert "effective branch scale is not a positive finite number" in proc.stderr
    assert "Traceback" not in proc.stderr
