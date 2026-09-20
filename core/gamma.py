"""Discrete and continuous gamma rate variation (Yang 1994).

Model       : among-site rate variation with a Gamma(alpha, alpha) distribution,
              either discretised into equally probable bins or drawn per site.
Equation    : bin rate i is the conditional mean
              k * [P(alpha + 1, q_{i+1} alpha) - P(alpha + 1, q_i alpha)],
              then all rates are normalised to mean 1.
Source      : Yang, Z. (1994) J. Mol. Evol. 39: 306-314. Seq-Gen 1.3.5
              implements the same conditional-mean convention and additionally
              exposes a median variant (source/gamma.c, DiscreteGamma).
Convention  : rates are conditional means by default. If the quantiles misbehave
              the code falls back to bin medians and then to an equal-rate
              repair, so a run always gets k positive rates. Reported rates are
              relative, and their mean over sites is 1.
Invariants  : the returned rates are positive and normalised to mean 1 for every
              (alpha, k) the caller can pass.
"""

import numpy as np
from scipy.special import gammainc
from scipy.stats import gamma as _gamma_dist

from ..utils import gamma as _rng_gamma


def rndgamma(shape):
    """Gamma(shape, 1) deviate from the seeded generator."""
    return _rng_gamma(shape)


def _quantile(prob, alpha):
    """Quantile of Gamma(alpha, rate=alpha); scipy takes the scale."""
    return _gamma_dist.ppf(prob, alpha, scale=1.0 / alpha)


def _mean_rates(alpha, k):
    """Conditional mean of each equiprobable bin (Yang 1994, the usual method).

    The bin mean is k * [P(alpha+1, q_i+1*alpha) - P(alpha+1, q_i*alpha)], where
    q is the bin boundary of Gamma(alpha, rate=alpha) and P is the regularised
    lower incomplete gamma.
    """
    factor = float(k)
    tail = []
    for i in range(k - 1):
        v = gammainc(alpha + 1.0, _quantile((i + 1.0) / k, alpha) * alpha)
        if not np.isfinite(v) or v < 0:
            raise ValueError("incomplete gamma out of range")
        tail.append(v)

    rates = [0.0] * k
    rates[0] = tail[0] * factor
    rates[k - 1] = (1.0 - tail[k - 2]) * factor
    for i in range(1, k - 1):
        rates[i] = (tail[i] - tail[i - 1]) * factor

    total = sum(rates)
    return [r * (factor / total) for r in rates]


def _median_rates(alpha, k):
    """Median of each equiprobable bin."""
    step = 1.0 / (2.0 * k)
    rates = [_quantile((2 * i + 1) * step, alpha) for i in range(k)]
    total = sum(rates)
    return [r * (k / total) for r in rates]


def _equal_rates(alpha, k):
    """Last resort when the quantiles misbehave: bin medians, non-positive
    entries lifted off zero, then renormalised."""
    step = 1.0 / (2.0 * k)
    rates = []
    for i in range(k):
        v = _quantile((2 * i + 1) * step, alpha)
        rates.append(v if (v >= 0 and np.isfinite(v)) else 0.0)

    positive = [r for r in rates if r > 0]
    if positive:
        floor = min(positive) * 0.5
        rates = [r if r > 0 else floor for r in rates]

    total = sum(rates)
    return [r * (k / total) for r in rates]


def _usable(rates):
    return (all(r >= 0 and np.isfinite(r) for r in rates)
            and any(r > 1e-12 for r in rates))


def discrete_gamma(alpha, k, median=False):
    """k equiprobable rates for Gamma(alpha, alpha), normalised to mean 1.

    Falls back to bin medians, then to the equal-rate repair, when the
    incomplete gamma evaluation is not usable. Returns (frequencies, rates).
    """
    first = _median_rates if median else _mean_rates
    for make in (first, _median_rates):
        try:
            rates = make(alpha, k)
        except (ValueError, ZeroDivisionError, OverflowError):
            continue
        if _usable(rates):
            return [1.0 / k] * k, rates
    return [1.0 / k] * k, _equal_rates(alpha, k)
