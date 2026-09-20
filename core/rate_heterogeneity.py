"""Per-site rate models: codon positions, gamma categories, invariable sites.

Model       : the per-site rate model for one run, held in one ``rate_info``
              mapping that both the evolution loop and the writers read.
Equation    : precedence is codon rates > discrete gamma > continuous gamma; the
              invariable mask is independent of all three. Variable sites are
              scaled by 1 / (1 - p) so the mean rate over all sites stays 1.
Source      : Yang, Z. (1994) for the gamma axis; Seq-Gen (Rambaut & Grassly
              1997) for the invariable-site convention, which scales the
              remaining sites rather than dropping the invariable ones.
Convention  : codon rates are relative. The command line normalises them to mean
              1 (sum 3) before they reach this module, so ``setup_rates`` uses
              them as given. The invariable proportion p is in [0, 1) and its
              mask is redrawn once per replicate. The 1/(1 - p) rescaling of the
              variable sites is applied in core/evolve.py, not here.
Invariants  : with no codon rates the mean relative rate over all sites is 1;
              discrete-gamma category rates are normalised to mean 1; the
              invariable mask has expected density p.
"""

from typing import Optional, TypedDict

import numpy as np

from .gamma import discrete_gamma, rndgamma
from ..tree.tree_parser import TTree
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


Partition = tuple
Dataset = list


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
