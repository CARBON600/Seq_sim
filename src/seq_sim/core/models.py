"""Substitution-model factory and per-node model registry.
Builds one model per (comp_idx, rmat_idx) pair referenced by [f=N] / [r=N]
branch annotations (Foster 2004). Index 0 on each axis is the primary CLI
model; further indices come from --alt-freqs, --alt-rates and --alt-tstv.
input: model name and its options, or the annotation pairs of a parsed tree
output: model instances and a registry keyed by (comp_idx, rmat_idx)
"""

from typing import Protocol

from numpy.typing import NDArray

from ._validators import normalise_freqs
from .aa_models import create_aa_model
from .nuc_models import F84Model, GTRModel, HKYModel, JC69Model

NUC_MODELS = {"HKY", "F84", "GTR", "REV", "JC69"}
AA_MODELS = {"JTT", "WAG", "PAM", "BLOSUM", "LG", "GENERAL"}


class SubstitutionModel(Protocol):
    """Interface shared by nucleotide and amino-acid models."""

    num_states: int
    state_chars: str
    freq: NDArray
    cum_freqs: NDArray

    def set_matrix(self, matrix, length):
        ...

    def set_vector(self, vector, state, length):
        ...


ModelRegistry = dict[tuple[int, int], SubstitutionModel]


def setup_model(name, tstv=None, equal_exchangeability=None, freqs=None,
                equal_freqs=False, relative_rates=None):
    """Instantiate one substitution model by name.

    equal_exchangeability is derived from tstv when it is left as None.
    """
    if equal_exchangeability is None:
        equal_exchangeability = tstv is None
    name = "GTR" if name.upper() == "REV" else name.upper()

    if name in NUC_MODELS:
        # equal_freqs and an absent -f mean the same thing on this path, so the
        # two routes agree; reading the flag stops a caller who passes both from
        # getting the supplied vector without an error.
        if equal_freqs or freqs is None:
            f = [0.25, 0.25, 0.25, 0.25]
        else:
            f = list(freqs)
            normalise_freqs(f)
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
    """Build one model per (comp_idx, rmat_idx) pair.

    A branch selects its pair with [f=N] / [r=N] on the branch leading into a
    node. Index 0 is the model built from the primary command line; further
    indices come from --alt-freqs, --alt-rates and --alt-tstv.
    """
    name = "GTR" if name.upper() == "REV" else name.upper()

    comps = [freqs, *(alt_freqs or [])]
    rmats = [relative_rates, *(alt_rates or [])]
    tstvs = [tstv, *(alt_tstv or [])]
    n_rmats = max(len(rmats), len(tstvs))

    registry = {}
    for c, comp in enumerate(comps):
        # The command line's equal_freqs only applies to axis entry 0;
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
