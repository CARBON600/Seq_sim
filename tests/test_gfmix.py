"""GFmix profile mixture: matrices, branch transform, site-profile assignment."""

import numpy as np
import pytest

from seq_sim.core.aa_models import AA_ORDER, empirical_rates
from seq_sim.core.cseries import CSERIES, load_cseries
from seq_sim.core.profile_mixture import (
    FYMINK,
    GARP,
    assign_site_profiles,
    branch_transform,
    build_profile_models,
)


def _reference_transform(pi, gamma_g, gamma_f):
    """Independent re-statement of Eq. 1.4 for cross-checking."""
    pi = np.asarray(pi, dtype=np.float64)
    g = [AA_ORDER.index(a) for a in GARP]
    f = [AA_ORDER.index(a) for a in FYMINK]
    sum_g = pi[g].sum()
    sum_f = pi[f].sum()
    sum_o = 1.0 - sum_g - sum_f
    mu = 1.0 / (gamma_g * sum_g + gamma_f * sum_f + sum_o)
    out = pi * mu
    out[g] *= gamma_g
    out[f] *= gamma_f
    return out


class TestCSeriesData:
    def test_all_models_load_and_are_consistent(self):
        for name in CSERIES:
            freqs, weights = load_cseries(name)
            c = int(name[1:])
            assert freqs.shape == (c, 20)
            assert weights.shape == (c,)
            assert np.all(freqs >= 0)
            assert np.all(weights >= 0)
            assert np.allclose(freqs.sum(axis=1), 1.0, atol=1e-6)
            assert np.isclose(weights.sum(), 1.0, atol=1e-6)

    def test_unknown_model_raises(self):
        with pytest.raises(ValueError):
            load_cseries("C99")


class TestBranchTransform:
    def test_renormalises_to_one(self):
        freqs, _ = load_cseries("C20")
        for c in range(freqs.shape[0]):
            out = branch_transform(freqs[c], 2.0, 0.5)
            assert np.isclose(out.sum(), 1.0, atol=1e-12)

    def test_sums_to_one_over_the_grid(self):
        """The bound quoted in REFERENCES.md, over the same grid."""
        grid = [(g, f) for g in (0.5, 1.0, 1.5, 2.0, 2.5)
                for f in (0.5, 1.0, 1.5, 2.0, 2.5)]
        assert len(grid) == 25
        worst = 0.0
        for name in CSERIES:
            freqs, _ = load_cseries(name)
            for row in freqs:
                p = row / row.sum()
                for g, f in grid:
                    worst = max(worst,
                                abs(float(branch_transform(p, g, f).sum()) - 1.0))
        assert worst <= 5e-15, worst

    def test_shipped_rows_keep_the_storage_precision(self):
        """The rows as shipped are 10-digit, so the transform cannot do better."""
        worst = 0.0
        for name in CSERIES:
            freqs, _ = load_cseries(name)
            for row in freqs:
                worst = max(worst,
                            abs(float(branch_transform(row, 0.5, 0.5).sum()) - 1.0))
        assert worst <= 1e-9, worst

    def test_identity_reproduces_the_input(self):
        worst = 0.0
        for name in CSERIES:
            freqs, _ = load_cseries(name)
            for row in freqs:
                p = row / row.sum()
                worst = max(worst, float(np.max(np.abs(branch_transform(p, 1.0, 1.0) - p))))
        assert worst <= 1e-15, worst

    def test_matches_reference(self):
        freqs, _ = load_cseries("C20")
        rng = np.random.default_rng(0)
        for c in range(freqs.shape[0]):
            gg, gf = rng.uniform(0.1, 5.0, 2)
            got = branch_transform(freqs[c], gg, gf)
            want = _reference_transform(freqs[c], gg, gf)
            assert np.allclose(got, want, atol=1e-14)

    def test_identity_when_gamma_one(self):
        freqs, _ = load_cseries("C20")
        out = branch_transform(freqs[0], 1.0, 1.0)
        assert np.allclose(out, freqs[0], atol=1e-14)

    def test_class_ratio_scales_correctly(self):
        freqs, _ = load_cseries("C20")
        pi = freqs[0]
        g = [AA_ORDER.index(a) for a in GARP]
        f = [AA_ORDER.index(a) for a in FYMINK]
        gg, gf = 3.0, 0.7
        out = branch_transform(pi, gg, gf)
        ratio_out = out[g].sum() / out[f].sum()
        ratio_in = (gg * pi[g].sum()) / (gf * pi[f].sum())
        assert np.isclose(ratio_out, ratio_in, atol=1e-12)


class TestAssignSiteProfiles:
    def test_assignment_counts_follow_weights(self):
        weights = np.array([0.1, 0.3, 0.6])
        rng = np.random.default_rng(42)
        sites = assign_site_profiles(200000, weights, rng.random)
        counts = np.bincount(sites, minlength=3)
        assert np.allclose(counts / len(sites), weights, atol=0.005)

    def test_indices_in_range(self):
        weights = np.array([0.1, 0.2, 0.7])
        rng = np.random.default_rng(1)
        sites = assign_site_profiles(1000, weights, rng.random)
        assert sites.min() >= 0 and sites.max() < 3


class TestInternalNodeLabels:
    def test_internal_labels_are_parsed_and_gamma_applied(self):
        from seq_sim.tree.tree_parser import parse_trees
        # inner node labelled 'H', tip A also heterogeneous
        text = "((A:0.1,B:0.1)H:0.2,(C:0.1,D:0.1):0.2);"
        (_, _, tree), = parse_trees(text)
        names = {n.name for n in tree.node_list if n.name}
        assert "A" in names and "H" in names
        # locate the labelled internal node and check gamma can be set by name
        het = {n.name: (2.0, 0.5) for n in tree.node_list if n.name in ("A", "H")}
        for node in tree.node_list:
            if node.name in het:
                node.gamma_g, node.gamma_f = het[node.name]
        labelled = [n for n in tree.node_list if n.name == "H"][0]
        assert labelled.gamma_g == 2.0 and labelled.gamma_f == 0.5


class TestBuildProfileModels:
    def test_base_models_use_profile_freqs(self):
        freqs, _ = load_cseries("C10")
        rates = empirical_rates("LG")
        base, per_node, _ = build_profile_models("C10", rates, {})
        assert len(base) == 10
        assert per_node == {}
        assert np.allclose(base[0].freq, freqs[0], atol=1e-12)

    def test_transformed_models_match_branch_transform(self):
        freqs, _ = load_cseries("C10")
        rates = empirical_rates("LG")
        gg, gf = 2.0, 0.4
        base, per_node, _ = build_profile_models("C10", rates, {"n1": (gg, gf)})
        assert len(per_node["n1"]) == 10
        for c in range(10):
            want = branch_transform(freqs[c], gg, gf)
            assert np.allclose(per_node["n1"][c].freq, want, atol=1e-12)
