"""Root sequence sampling, pre-order tree traversal, and branch mutation loop.
Sequences are evolved from the root out to the tips by sampling child states
from each branch's cumulative transition probability matrix, or from per-site
profile models when --gfmix-cseries is active.
input: a parsed tree, a model registry, rate_info, optional ancestor and profiles
output: sequences stored on the tips and nodes of that tree
"""

import sys

import numpy as np

from ..utils import random
from ._validators import validate_cdf_matrix, validate_cdf_row

###function definition

def evolve_sequences(tree, from_site, num_sites, scale, model_registry,
                     rate_info, ancestor=None, profile_info=None):
    """Evolve one partition from the root out to the tips.

    ancestor seeds the root sequence; without it the root is drawn from the
    root model's equilibrium frequencies. With profile_info the root is
    instead drawn per-site from each site's base profile frequencies.
    """
    if rate_info.get("invariable"):
        scale /= (1.0 - rate_info["invariable_prop"])

    root_seq = tree.root.sequence
    if profile_info is None:
        root_model = _lookup_model(model_registry, tree.root)
        if ancestor is not None:
            _set_sequence(root_seq, ancestor, from_site, num_sites, root_model)
        else:
            _random_sequence(root_seq, num_sites, root_model.cum_freqs)
    else:
        base_models = profile_info["base_models"]
        # The root has its own composition (GFmix Eq. 9: the root edge carries a
        # gamma vector), so a root entry in per_node must be used here too.
        root_models = profile_info["per_node"].get(tree.root, base_models)
        if ancestor is not None:
            _set_sequence(root_seq, ancestor, from_site, num_sites,
                          root_models[0])
        else:
            _random_sequence_profiles(root_seq, from_site, num_sites,
                                      root_models,
                                      profile_info["site_profile"])

    _set_site_rates(rate_info, from_site, num_sites, scale)

    _evolve_node(tree.root, tree.root.branch1, from_site, num_sites,
                 scale, model_registry, rate_info, profile_info)
    _evolve_node(tree.root, tree.root.branch2, from_site, num_sites,
                 scale, model_registry, rate_info, profile_info)
    if not tree.rooted:
        _evolve_node(tree.root, tree.root.branch0, from_site, num_sites,
                     scale, model_registry, rate_info, profile_info)


def _raise_recursion_limit(num_nodes):
    """Deep trees need more stack than the default limit allows."""
    min_limit = num_nodes + 100
    if sys.getrecursionlimit() < min_limit:
        sys.setrecursionlimit(max(min_limit, 2000))


def _lookup_model(model_registry, node):
    key = (node.comp_idx, node.rmat_idx)
    try:
        return model_registry[key]
    except KeyError:
        raise ValueError(
            f"No substitution model registered for comp_idx={key[0]}, "
            f"rmat_idx={key[1]}. Available keys: {sorted(model_registry)}. "
            f"Check the [f=N]/[r=N] annotations against --alt-freqs/"
            f"--alt-rates/--alt-tstv.") from None


def _random_sequence(seq, num_sites, cum_freqs):
    """Draw the root sequence from cumulative equilibrium frequencies."""
    last = len(cum_freqs) - 1
    for i in range(num_sites):
        r = random()
        for s, af in enumerate(cum_freqs):
            if r < af:
                seq[i] = s
                break
        else:
            # cum_freqs can land one ULP below 1.0; treat those draws as the
            # final state rather than leaving the site unset.
            seq[i] = last


def _set_sequence(seq, source, from_site, num_sites, model):
    """Copy an ancestral string into the node's integer buffer."""
    char_map = {c: i for i, c in enumerate(model.state_chars)}
    for i in range(num_sites):
        ch = source[from_site + i]
        idx = char_map.get(ch.upper())
        if idx is None:
            raise ValueError(f"Invalid state '{ch}' in ancestral sequence")
        seq[i] = idx


def _evolve_node(anc, des, from_site, num_sites, scale, model_registry,
                 rate_info, profile_info=None):
    """Copy the parent sequence onto des, mutate it, then descend."""
    if profile_info is None:
        branch_model = _lookup_model(model_registry, des)
    else:
        branch_model = profile_info["per_node"].get(
            des, profile_info["base_models"])
    site_profile = None if profile_info is None else profile_info["site_profile"]
    for _what, _val in (("branch length", des.length0),
                        ("branch rate multiplier", des.branch_rate),
                        ("global branch scale", scale)):
        if not np.isfinite(_val):
            raise ValueError(f"{_what} is not finite: {_val!r}")
        if _val < 0.0:
            raise ValueError(f"{_what} is negative: {_val!r}")
    length = des.length0 * scale * des.branch_rate
    des.sequence[:] = anc.sequence[:]
    _mutate_sequence(des.sequence, from_site, num_sites, length,
                     branch_model, rate_info, site_profile)
    if des.tip_no == -1:
        _evolve_node(des, des.branch1, from_site, num_sites,
                     scale, model_registry, rate_info, profile_info)
        _evolve_node(des, des.branch2, from_site, num_sites,
                     scale, model_registry, rate_info, profile_info)


def _mutate_sequence(seq, from_site, num_sites, length, model, rate_info,
                     site_profile=None):
    """Sample a new state for every site of one branch.

    With site_profile, model is a list of per-profile models and every
    site evolves under its own profile; otherwise model is a single model
    shared by all sites.
    """
    if site_profile is not None:
        _mutate_sequence_profiles(seq, from_site, num_sites, length, model,
                                  rate_info, site_profile)
        return

    htype = rate_info["type"]
    n_states = model.num_states
    inv = _inv_mask(rate_info, from_site, num_sites)

    if htype == "discrete_gamma":
        cat_rates = rate_info["cat_rates"]
        num_cats = rate_info["num_cats"]
        matrices = [np.zeros(n_states * n_states) for _ in range(num_cats)]
        for c in range(num_cats):
            model.set_matrix(matrices[c], cat_rates[c] * length)
        for m in matrices:
            validate_cdf_matrix(m, n_states)
        cats = rate_info["categories"][from_site:from_site + num_sites]
        for i in range(num_sites):
            if inv is not None and inv[i]:
                continue
            row = int(seq[i]) * n_states
            seq[i] = _pick_state_from_matrix(matrices[cats[i]], row, n_states)
        return

    if htype == "codon":
        cats = np.array(rate_info["cat_rates"], dtype=np.float64)
        site_lengths = cats[np.arange(num_sites) % 3] * length
    elif htype == "continuous_gamma":
        site_lengths = (rate_info["gamma_rates"][from_site:from_site + num_sites]
                        * length)
    elif htype == "none":
        site_lengths = None          # uniform: every site gets length
    else:
        raise ValueError(f"unknown rate heterogeneity type: {htype}")

    vec = np.zeros(n_states)
    for i in range(num_sites):
        if inv is not None and inv[i]:
            continue
        site_length = length if site_lengths is None else float(site_lengths[i])
        model.set_vector(vec, int(seq[i]), site_length)
        seq[i] = _pick_state(vec)


def _random_sequence_profiles(seq, from_site, num_sites, models, site_profile):
    """Draw the root sequence per-site from each site's profile frequencies."""
    profs = site_profile[from_site:from_site + num_sites]
    for i in range(num_sites):
        af = models[profs[i]].cum_freqs
        last = len(af) - 1
        r = random()
        for st, v in enumerate(af):
            if r < v:
                seq[i] = st
                break
        else:
            seq[i] = last


def _mutate_sequence_profiles(seq, from_site, num_sites, length, models,
                              rate_info, site_profile):
    """Branch mutation under a profile mixture.

    models is a length-C list sharing one exchangeability matrix but differing
    in equilibrium frequencies; site_profile maps every site to its profile
    index. Discrete gamma precomputes C*K matrices, the other rate models build
    one cumulative row per site.
    """
    htype = rate_info["type"]
    n_states = models[0].num_states
    inv = _inv_mask(rate_info, from_site, num_sites)
    profs = site_profile[from_site:from_site + num_sites]

    if htype == "discrete_gamma":
        cat_rates = rate_info["cat_rates"]
        num_cats = rate_info["num_cats"]
        n_prof = len(models)
        matrices = [np.zeros(n_states * n_states)
                    for _ in range(n_prof * num_cats)]
        for p in range(n_prof):
            for c in range(num_cats):
                models[p].set_matrix(matrices[p * num_cats + c],
                                     cat_rates[c] * length)
        for m in matrices:
            validate_cdf_matrix(m, n_states)
        cats = rate_info["categories"][from_site:from_site + num_sites]
        for i in range(num_sites):
            if inv is not None and inv[i]:
                continue
            row = int(seq[i]) * n_states
            m = matrices[profs[i] * num_cats + cats[i]]
            seq[i] = _pick_state_from_matrix(m, row, n_states)
        return

    if htype == "continuous_gamma":
        site_lengths = (rate_info["gamma_rates"][from_site:from_site + num_sites]
                        * length)
    elif htype == "none":
        site_lengths = None
    else:
        # "codon" is nucleotide-only and cannot reach a profile-mixture path.
        raise ValueError(f"unexpected rate heterogeneity type: {htype}")

    vec = np.zeros(n_states)
    for i in range(num_sites):
        if inv is not None and inv[i]:
            continue
        site_length = length if site_lengths is None else float(site_lengths[i])
        models[profs[i]].set_vector(vec, int(seq[i]), site_length)
        seq[i] = _pick_state(vec)


def _inv_mask(rate_info, from_site, num_sites):
    if rate_info.get("invariable_sites") is None:
        return None
    return rate_info["invariable_sites"][from_site:from_site + num_sites]


def _pick_state(cumulative_vec):
    """Draw a state from a cumulative probability vector."""
    validate_cdf_row(cumulative_vec)
    r = random()
    for s, v in enumerate(cumulative_vec):
        if r < v:
            return s
    # A valid CDF ends at 1.0, so the loop always returns.  Falling through
    # means the row held NaN/inf (every comparison false): refuse to guess.
    raise ValueError("transition row is not a valid CDF: " + repr(list(cumulative_vec)))


def _pick_state_from_matrix(matrix, row_start, n_states=4):
    """Draw a state from one cumulative row of a flat probability matrix."""
    r = random()
    for j in range(n_states):
        if r < matrix[row_start + j]:
            return j
    raise ValueError("transition row is not a valid CDF at offset "
                     + repr(row_start))


def _set_site_rates(rate_info, from_site, num_sites, scale):
    """Fill the per-site rate table that -wr reports."""
    site_rates = rate_info.get("site_rates")
    if site_rates is None:
        return

    htype = rate_info["type"]
    cat_rates = rate_info.get("cat_rates")
    if htype == "continuous_gamma":
        rates = (rate_info["gamma_rates"][from_site:from_site + num_sites]
                 * scale)
    elif htype == "discrete_gamma":
        cats = rate_info["categories"][from_site:from_site + num_sites]
        rates = [cat_rates[c] * scale for c in cats]
    elif htype == "codon":
        rates = [cat_rates[i % 3] * scale for i in range(num_sites)]
    else:
        rates = [scale] * num_sites

    out = site_rates[from_site:from_site + num_sites]
    inv = _inv_mask(rate_info, from_site, num_sites)
    if inv is None:
        out[:] = rates
    else:
        for i in range(num_sites):
            out[i] = 0.0 if inv[i] else rates[i]
