"""Argument checks shared by the substitution models."""

import math

import numpy as np

###function definition

def _check_positive(freqs, kind, order=None):
    for i, pi in enumerate(freqs):
        if not np.isfinite(pi) or pi <= 0.0:
            label = f" ({order[i]})" if order else ""
            raise ValueError(
                f"all {kind} frequencies must be finite and positive, "
                f"got freq[{i}]={pi}{label}")


def _check_sums_to_one(freqs, kind):
    total = sum(freqs)
    if not np.isclose(total, 1.0, rtol=1e-10, atol=1e-12):
        raise ValueError(f"{kind} frequencies must sum to 1.0, got {total}")


def validate_nuc_freqs(freqs):
    """Four A,C,G,T frequencies, positive and summing to 1."""
    if len(freqs) != 4:
        raise ValueError(
            f"freqs must contain exactly four values in A,C,G,T order, "
            f"got {len(freqs)}")
    _check_positive(freqs, "base")
    _check_sums_to_one(freqs, "base")


def validate_aa_freqs(freqs, aa_order=None):
    """Amino-acid frequencies, positive and summing to 1."""
    _check_positive(freqs, "AA", aa_order)
    _check_sums_to_one(freqs, "AA")


def check_relative_rates(rates, expected, order=None):
    """Exchangeability values: right count, finite, non-negative, not all zero."""
    if len(rates) != expected:
        tail = (f" values in order {order}" if order else " values")
        raise ValueError(f"relative_rates must contain exactly {expected}"
                         + tail + f", got {len(rates)}")
    for i, r in enumerate(rates):
        if not np.isfinite(r):
            raise ValueError(f"relative_rates must be finite, got rate[{i}]={r}")
        if r < 0.0:
            raise ValueError("relative_rates must be non-negative, "
                             f"got rate[{i}]={r}")
    if all(r == 0.0 for r in rates):
        raise ValueError("at least one relative rate must be positive")


def validate_state(state, num_states=4):
    """Return state as an int, checked against [0, num_states)."""
    if isinstance(state, float) and not float(state).is_integer():
        raise ValueError(f"state must be an integer, got {state}")
    state = int(state)
    if state < 0 or state >= num_states:
        raise ValueError(f"state must be in [0, {num_states}), got {state}")
    return state


def validate_length(length):
    """Return length as a finite, non-negative float."""
    try:
        length = float(length)
    except (TypeError, ValueError):
        raise ValueError(f"length must be a real number, got {length!r}")
    if not np.isfinite(length):
        raise ValueError(f"length must be finite, got {length}")
    if length < 0.0:
        raise ValueError(f"length must be non-negative, got {length}")
    return length


def validate_output_array(arr, expected_size):
    """Check arr is a writable float array of expected_size elements."""
    if arr.size != expected_size:
        raise ValueError(f"output array must contain exactly {expected_size} "
                         f"elements, got {arr.size}")
    if not np.issubdtype(arr.dtype, np.floating):
        raise TypeError(f"output array must have floating dtype, got {arr.dtype}")
    if not arr.flags.writeable:
        raise ValueError("output array must be writable")


def normalise_freqs(arr):
    """Rescale arr in place to sum to 1."""
    for x in arr:
        if not np.isfinite(x):
            raise ValueError(f"frequencies must be finite, got {x!r}")
        if x < 0.0:
            raise ValueError(f"frequencies must be non-negative, got {x!r}")
    total = sum(arr)
    if total == 0.0:
        raise ValueError("cannot normalise frequencies: sum is zero")
    for i in range(len(arr)):
        arr[i] /= total


def validate_cdf_row(row):
    for v in row:
        if not math.isfinite(v):
            raise ValueError("transition row is not finite")
    prev = row[0]
    for v in row[1:]:
        if v - prev < -1e-12:
            raise ValueError("transition row is not monotone")
        prev = v
    if abs(float(row[-1]) - 1.0) > 1e-9:
        raise ValueError("transition row does not end at 1")


def validate_cdf_matrix(matrix, n_states):
    """Check every row of a flat n_states*n_states cumulative matrix."""
    for start in range(0, len(matrix), n_states):
        validate_cdf_row(matrix[start:start + n_states])
