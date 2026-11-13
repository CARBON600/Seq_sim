"""Per-site rate variation: codon positions, gamma categories, invariable sites.

Configures the ``rate_info`` state dictionary consumed by the evolution loop
and rate writer, and redraws per-site category assignments and invariable-site
masks at the start of each replicate.
"""

from typing import TypedDict

import numpy as np

from .gamma import discrete_gamma, rndgamma
from ..utils import random as _rand


class RateInfo(TypedDict, total=False):
    """Description of the per-site rate model for one run."""

    type: str
    invariable: bool
    invariable_prop: float
    invariable_sites: np.ndarray
    cat_rates: list
    num_cats: int
    gamma_shape: float
    categories: np.ndarray
    gamma_rates: np.ndarray
    site_rates: np.ndarray


def setup_rates(num_sites, *, codon_rates=None, gamma_shape=None,
                gamma_cats=None, invariable=None, write_rates=False):
    """Describe the per-site rate model for one run.

    Precedence is codon positions > discrete gamma > continuous gamma; the
    invariable-site mask is independent of all three. Returns the ``rate_info``
    mapping that the evolution loop and the writers read.
    """
    info = {
        "type": "none",
        "invariable": invariable is not None,
        "invariable_prop": invariable or 0.0,
    }
    if invariable is not None:
        info["invariable_sites"] = np.zeros(num_sites, dtype=np.int32)

    if codon_rates:
        info.update(type="codon", cat_rates=list(codon_rates), num_cats=3)
    elif gamma_cats:
        info.update(type="discrete_gamma", num_cats=gamma_cats,
                    gamma_shape=gamma_shape)
        _, info["cat_rates"] = discrete_gamma(gamma_shape, gamma_cats)
        info["categories"] = np.zeros(num_sites, dtype=np.int32)
    elif gamma_shape is not None:
        info.update(type="continuous_gamma", gamma_shape=gamma_shape)
        info["gamma_rates"] = np.zeros(num_sites)

    if write_rates:
        info["site_rates"] = np.zeros(num_sites)
    return info


def assign_categories(rate_info, num_sites, gamma_shape, num_cats):
    """Redraw the invariable mask and the gamma draws for one replicate."""
    if rate_info.get("invariable_sites") is not None:
        prop = rate_info["invariable_prop"]
        inv = rate_info["invariable_sites"]
        for i in range(num_sites):
            inv[i] = 1 if _rand() < prop else 0

    htype = rate_info["type"]
    if htype == "continuous_gamma":
        gr = rate_info["gamma_rates"]
        for i in range(num_sites):
            gr[i] = rndgamma(gamma_shape) / gamma_shape
    elif htype == "discrete_gamma":
        cats = rate_info["categories"]
        for i in range(num_sites):
            cats[i] = int(_rand() * num_cats)


def normalize_branch_rates(datasets):
    """Rescale per-branch multipliers so the length-weighted mean is 1.0."""
    weight = 0.0
    weighted = 0.0
    for parts in datasets:
        for _, _, tree in parts:
            for node in tree.node_list:
                if node.branch0 is not None:
                    weight += node.length0
                    weighted += node.length0 * node.param
    if weighted <= 0:
        return
    scale = weight / weighted
    for parts in datasets:
        for _, _, tree in parts:
            for node in tree.node_list:
                node.param *= scale
