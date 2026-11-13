"""Regression tests for seq_sim_v2 (run with: ``pytest``)."""

import io
import sys

import numpy as np
import pytest
from scipy.linalg import expm

import seq_sim_v2.seq_sim as S
from seq_sim_v2.core.aa_models import create_aa_model, _EMPIRICAL_MODELS, normalize_freqs
from seq_sim_v2.core.models import build_model_registry, setup_model, AA_MODELS
from seq_sim_v2.core.evolve import _lookup_model, _set_site_rates
from seq_sim_v2.core.rate_heterogeneity import assign_categories, setup_rates
from seq_sim_v2.tree.tree_parser import parse_trees, allocate_sequences, ParseError
from seq_sim_v2.output.writers import write_ancestral
from seq_sim_v2.output.nexus_reader import (
    is_nexus, parse_nexus_trees, _is_fr_annotation, _strip_bracket_content)
from seq_sim_v2.utils import set_seed

TREE_ROOTED = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"


def _parse(argv):
    old = sys.argv
    sys.argv = ["seq_sim"] + argv
    try:
        return S._parse_args()
    finally:
        sys.argv = old


def _cli_equal_freqs(args):
    name = args.model.upper()
    if name == "REV":
        name = "GTR"
    return args.equal_freqs or (args.freqs is None and name not in AA_MODELS)


def _registry(args):
    return build_model_registry(
        args.model, freqs=args.freqs, tstv=args.tstv,
        equal_exchangeability=(args.tstv is None),
        equal_freqs=_cli_equal_freqs(args), relative_rates=args.rates,
        alt_freqs=args.alt_freqs_list, alt_tstv=args.alt_tstv,
        alt_rates=args.alt_rates_list)


# --- Empirical matrices ---

def test_empirical_matrices_shape_and_normalisation():
    for name in ["JTT", "WAG", "PAM", "BLOSUM", "LG"]:
        rates = _EMPIRICAL_MODELS[name]["rates"]
        freqs = list(_EMPIRICAL_MODELS[name]["freqs"])
        normalize_freqs(freqs)
        assert len(rates) == 190, name
        assert len(freqs) == 20, name
        assert all(r >= 0 for r in rates), name
        assert abs(sum(freqs) - 1.0) < 1e-9, name


def test_lg_matches_authoritative_values():
    lg = _EMPIRICAL_MODELS["LG"]["rates"]
    assert lg[0] == pytest.approx(0.425093)
    assert lg[3] == pytest.approx(2.489084)
    assert lg[18] == pytest.approx(2.547870)
    assert lg[56] == pytest.approx(5.24387)
    assert lg[144] == pytest.approx(10.649107)
    assert lg[-1] == pytest.approx(0.249313)


def test_aa_default_uses_empirical_frequencies():
    args = _parse(["-m", "WAG"])
    reg = _registry(args)
    freqs = list(_EMPIRICAL_MODELS["WAG"]["freqs"])
    normalize_freqs(freqs)
    assert np.allclose(reg[(0, 0)].freq, freqs)


# --- Transition-probability correctness ---

def _reconstruct_Q(model):
    Ssym = (model._eigvecs * model._roots[None, :]) @ model._eigvecs.T
    D = np.diag(model._sqrt_pi)
    return np.linalg.inv(D) @ Ssym @ D


@pytest.mark.parametrize("name", ["JTT", "WAG", "PAM", "BLOSUM", "LG"])
def test_aa_Pt_matches_matrix_exponential(name):
    model = create_aa_model(name)
    Q = _reconstruct_Q(model)
    n = model.num_states
    for t in (0.01, 0.1, 0.5, 1.0):
        cum = np.zeros(n * n)
        model.set_matrix(cum, t)
        raw = np.diff(cum.reshape(n, n), axis=1, prepend=0.0)
        assert np.allclose(raw, expm(Q * t), atol=1e-8), (name, t)


def test_gtr_Pt_matches_matrix_exponential():
    from seq_sim_v2.core.nuc_models import GTRModel
    model = GTRModel([0.1, 0.2, 0.3, 0.4], [1.0, 2.0, 0.5, 1.5, 0.3, 2.5])
    Q = _reconstruct_Q(model)
    cum = np.zeros(16)
    model.set_matrix(cum, 0.3)
    raw = np.diff(cum.reshape(4, 4), axis=1, prepend=0.0)
    assert np.allclose(raw, expm(Q * 0.3), atol=1e-8)


def test_discrete_gamma_mean_is_one():
    from seq_sim_v2.core.gamma import discrete_gamma
    for alpha in (0.1, 0.5, 1.0, 2.0):
        for k in (4, 8):
            freqK, rK = discrete_gamma(alpha, k)
            assert abs(sum(p * r for p, r in zip(freqK, rK)) - 1.0) < 1e-9


# --- -i (invariable sites) regression ---

def test_invariant_sites_assigned_when_alone():
    set_seed(1)
    info = setup_rates(100000, invariable=0.5)
    assert info["type"] == "none"
    assign_categories(info, 100000, None, None)
    frac = float(np.mean(info["invariable_sites"]))
    assert abs(frac - 0.5) < 0.02


def test_invariant_sites_zero_rate_in_wr():
    set_seed(2)
    info = setup_rates(100000, invariable=0.5, write_rates=True)
    assign_categories(info, 100000, None, None)
    _set_site_rates(info, 0, 100000, 1.0)
    zeros = int(np.sum(info["site_rates"] == 0.0))
    assert abs(zeros - 50000) < 2000


# --- -wa (ancestral output) regression ---

@pytest.mark.parametrize("fmt", ["p", "f"])
def test_ancestral_writes_each_node_once(fmt):
    _, _, tree = parse_trees(TREE_ROOTED)[0]
    allocate_sequences(tree, 8)
    buf = io.StringIO()
    write_ancestral(buf, fmt, tree, "ACGT")
    text = buf.getvalue()
    if fmt == "p":
        lines = [ln for ln in text.splitlines() if ln.strip()]
        assert int(lines[0].split()[0]) == tree.num_nodes
        assert len(lines) - 1 == tree.num_nodes
    else:
        assert text.count(">") == tree.num_nodes


def test_ancestral_nexus_block():
    _, _, tree = parse_trees(TREE_ROOTED)[0]
    allocate_sequences(tree, 8)
    buf = io.StringIO()
    write_ancestral(buf, "n", tree, "ACGT")
    text = buf.getvalue()
    assert "Begin DATA" in text
    assert f"NTAX={tree.num_nodes} NCHAR=8" in text
    assert "DATATYPE=DNA" in text
    assert "End;" in text


def test_ancestral_relaxed_phylip_full_names():
    _, _, tree = parse_trees(TREE_ROOTED)[0]
    allocate_sequences(tree, 4)
    buf = io.StringIO()
    write_ancestral(buf, "r", tree, "ACGT")
    lines = [ln for ln in buf.getvalue().splitlines() if ln.strip()]
    assert int(lines[0].split()[0]) == tree.num_nodes
    assert "A AAAA" in lines          # tip: full name, single space
    assert "5 AAAA" in lines          # internal node numbered from num_tips+1


# --- Parser and registry diagnostics ---

def test_annotation_before_branch_length():
    _, _, tree = parse_trees("(A[f=1]:0.1,B:0.1);")[0]
    assert tree.tips[0].comp_idx == 1
    assert tree.tips[0].length0 == pytest.approx(0.1)


def test_annotation_after_branch_length():
    _, _, tree = parse_trees("(A:0.1[f=2],B:0.1);")[0]
    assert tree.tips[0].comp_idx == 2


def test_registry_lookup_diagnostic():
    class _Node:
        comp_idx = 1
        rmat_idx = 0

    with pytest.raises(ValueError, match="comp_idx=1"):
        _lookup_model({(0, 0): setup_model("JC69")}, _Node())


# --- NEXUS annotations ---

def test_nexus_preserves_fr_annotations():
    text = (
        "#NEXUS\nBEGIN TREES;\n"
        "TREE t1 = ((A:0.1[f=1],B:0.1):0.1,C:0.2);\n"
        "END;\n"
    )
    assert is_nexus(text)
    _, trees = parse_nexus_trees(text)
    newick = trees[0][0]
    assert "[f=1]" in newick
    _, _, tree = parse_trees(newick)[0]
    assert tree.tips[0].comp_idx == 1


# --- CLI vector separators ---

def test_rate_separators():
    assert _parse(["-m", "GTR", "-r", "1,2,3,4,5,6"]).rates == [1, 2, 3, 4, 5, 6]
    assert _parse(["-m", "GTR", "-r", "1", "2", "3", "4", "5", "6"]).rates == [1, 2, 3, 4, 5, 6]


def test_freqs_separators():
    assert _parse(["-m", "HKY", "-f", "0.2,0.3,0.3,0.2"]).freqs == pytest.approx([0.2, 0.3, 0.3, 0.2])
    args = _parse(["-m", "WAG", "-f", "e"])
    assert args.equal_freqs is True and args.freqs is None


# --- main() validation and output-directory hygiene ---

def _run_main(argv, tree_text, tmp_path, monkeypatch):
    treefile = tmp_path / "t.nwk"
    treefile.write_text(tree_text, encoding="utf-8")
    out = tmp_path / "out"
    monkeypatch.setattr(
        sys, "argv", ["seq_sim"] + argv + ["-y", str(out), str(treefile)])
    with pytest.raises(SystemExit):
        S.main()
    return out


def test_tree_scale_zero_height_errors_without_leftover_dir(tmp_path, monkeypatch):
    out = _run_main(
        ["-m", "HKY", "-d", "0.5"],
        "((A:0,B:0):0,(C:0,D:0):0);", tmp_path, monkeypatch)
    assert not out.exists()


def test_tree_scale_unrooted_errors_without_leftover_dir(tmp_path, monkeypatch):
    out = _run_main(
        ["-m", "HKY", "-d", "0.5"],
        "((A:0.1,B:0.1):0.1,C:0.1,D:0.1);", tmp_path, monkeypatch)
    assert not out.exists()


def test_ancestor_length_validation(tmp_path, monkeypatch):
    _run_main(
        ["-m", "HKY", "-l", "10", "--ancestor", "ACGT"],
        TREE_ROOTED, tmp_path, monkeypatch)


def test_ancestor_invalid_state_validation(tmp_path, monkeypatch):
    _run_main(
        ["-m", "HKY", "-l", "4", "--ancestor", "ACGX"],
        TREE_ROOTED, tmp_path, monkeypatch)


def test_aa_rejects_codon_rates():
    with pytest.raises(SystemExit):
        _parse(["-m", "WAG", "-c", "1", "1", "1"])


def test_treefile_dash_reads_stdin(tmp_path, monkeypatch):
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", ["seq_sim", "-m", "JC69", "-l", "8",
                                      "-z", "1", "-y", str(out), "-"])
    monkeypatch.setattr(sys, "stdin", io.StringIO(TREE_ROOTED))
    S.main()
    assert (out / "out_1.fasta").exists()


# --- NEXUS f/r annotation filtering ---

@pytest.mark.parametrize("content,expected", [
    ("f=1", True), ("r=0", True), ("f=1,r=2", True), ("F = 3", True),
    ("R", False), ("F", False), ("&R", False), ("height=0.5", False),
    ("f=", False), ("f=1,R", False),
])
def test_fr_annotation_detection(content, expected):
    assert _is_fr_annotation(content) is expected


def test_nexus_strips_plain_comments_but_keeps_fr():
    text = "((A:0.1[R],B:0.1[f=1]):0.1,C:0.2);"
    out = _strip_bracket_content(text, keep_fr=True)
    assert "[R]" not in out
    assert "[f=1]" in out


# --- Duplicate branch metadata validation ---

def test_duplicate_branch_length_rejected():
    with pytest.raises(ParseError):
        parse_trees("(A:0.1:0.2,B:0.1);")


def test_duplicate_f_annotation_rejected():
    with pytest.raises(ParseError):
        parse_trees("(A:0.1[f=1][f=2],B:0.1);")


def test_duplicate_rate_annotation_rejected():
    with pytest.raises(ParseError):
        parse_trees("(A:0.1[1.5][2.5],B:0.1);")


def test_multiple_distinct_annotations_allowed():
    _, _, tree = parse_trees("(A:0.1[f=1][r=2],B:0.1);")[0]
    assert tree.tips[0].comp_idx == 1
    assert tree.tips[0].rmat_idx == 2


# --- Input file encoding fallback ---

def test_read_text_utf8_bom(tmp_path):
    p = tmp_path / "bom.nwk"
    p.write_bytes("\ufeff((A:0.1,B:0.1):0.1);".encode("utf-8"))
    text = S._read_text(str(p))
    assert not text.startswith("\ufeff")
    assert text.startswith("(")


def test_read_text_non_utf8_fallback(tmp_path, capsys):
    p = tmp_path / "gbk.nwk"
    p.write_bytes("((\u4e2d\u6587:0.1,B:0.1):0.1);".encode("gbk"))
    text = S._read_text(str(p))
    assert "0.1" in text
    assert "not valid UTF-8" in capsys.readouterr().err
