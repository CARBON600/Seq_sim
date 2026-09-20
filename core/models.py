"""Substitution-model factory and the per-node model registry.

Model       : one substitution model per (composition, rate-matrix) pair. A node
              selects a pair through the [f=N] and [r=N] annotations on the
              branch that leads into it.
Equation    : not applicable. This module assembles models; the transition
              probabilities live in the model classes, and the per-site rates in
              core/rate_heterogeneity.py.
Source      : the idea is Foster, P. G. (2004) Modeling compositional
              heterogeneity, Systematic Biology 53(3): 485-495, and its reference
              implementation P4. See docs/REFERENCES.md. The
              (comp_idx, rmat_idx) keying and the annotation syntax are this
              project's own; P4 is GPL v2 and no P4 code is used here.
Convention  : index 0 on each axis is the model built from the command line; the
              remaining entries come from --alt-freqs, --alt-rates and
              --alt-tstv. equal_freqs applies to axis entry 0 only; a later
              composition is equal-frequency solely when none was supplied. JC69
              ignores both axes, and HKY/F84 fall back to equal exchangeability
              whenever no ts/tv ratio is given for that matrix.
Invariants  : (0, 0) is always present, and the registry holds exactly
              len(compositions) * max(len(rate matrices), len(ts/tv ratios))
              entries. A key outside it raises a ValueError naming the available
              keys, so an annotation typo cannot silently fall back to another
              model.
"""

from typing import Protocol

from numpy.typing import NDArray

from ._validators import normalize_freqs
from .aa_models import create_aa_model
from .nuc_models import F84Model, GTRModel, HKYModel, JC69Model

NUC_MODELS = {"HKY", "F84", "GTR", "REV", "JC69"}
AA_MODELS = {"JTT", "WAG", "PAM", "BLOSUM", "LG", "GENERAL"}


class SubstitutionModel(Protocol):
    """Interface shared by nucleotide and amino-acid models."""

    num_states: int
    state_chars: str
    freq: NDArray
    add_freq: NDArray

    def set_matrix(self, matrix: NDArray, length: float) -> None:
        ...

    def set_vector(self, vector: NDArray, state: int, length: float) -> None:
        ...


ModelRegistry = dict[tuple[int, int], SubstitutionModel]


def setup_model(name, tstv=None, equal_exchangeability=True, freqs=None,
                equal_freqs=False, relative_rates=None):
    """Instantiate one substitution model by name."""
    name = "GTR" if name.upper() == "REV" else name.upper()

    if name in NUC_MODELS:
        f = [0.25, 0.25, 0.25, 0.25] if freqs is None else list(freqs)
        if freqs is not None:
            normalize_freqs(f)
        if name == "JC69":
            return JC69Model()
        if name == "HKY":
            return HKYModel(f, tstv, equal_exchangeability)
        if name == "F84":
            return F84Model(f, tstv, equal_exchangeability)
        if name == "GTR":
            if relative_rates is None:
                relative_rates = [1.0] * 6
            return GTRModel(f, list(relative_rates))
    elif name in AA_MODELS:
        return create_aa_model(name, user_freqs=freqs,
                               user_rates=relative_rates,
                               equal_freqs=equal_freqs)
    raise ValueError(f"Unknown model: {name}")


def build_model_registry(name, freqs, tstv, equal_exchangeability, equal_freqs,
                         relative_rates, alt_freqs, alt_tstv, alt_rates):
    """One model per (composition, rate-matrix) pair.

    Trees pick a pair per branch with ``[f=N]`` / ``[r=N]`` on the branch that
    leads into a node. Index 0 on each axis is the model built from the command
    line's own arguments; the remaining entries come from --alt-freqs,
    --alt-rates and --alt-tstv.
    """
    name = "GTR" if name.upper() == "REV" else name.upper()

    comps = [freqs, *(alt_freqs or [])]
    rmats = [relative_rates, *(alt_rates or [])]
    tstvs = [tstv, *(alt_tstv or [])]
    n_rmats = max(len(rmats), len(tstvs))

    registry = {}
    for c, comp in enumerate(comps):
        # The command line's ``equal_freqs`` only applies to axis entry 0;
        # extra entries are equal-frequency solely when none was supplied.
        eq_freqs = equal_freqs if c == 0 else comp is None
        for r in range(n_rmats):
            if name == "JC69":
                registry[(c, r)] = JC69Model()
            elif name in ("HKY", "F84"):
                tv = tstvs[r] if r < len(tstvs) else None
                eq_exch = equal_exchangeability if (c == 0 and r == 0) else tv is None
                registry[(c, r)] = setup_model(
                    name, tstv=tv, equal_exchangeability=eq_exch, freqs=comp,
                    equal_freqs=eq_freqs)
            else:
                rate = rmats[r] if r < len(rmats) else None
                registry[(c, r)] = setup_model(
                    name, freqs=comp, equal_freqs=eq_freqs,
                    relative_rates=rate)
    return registry
