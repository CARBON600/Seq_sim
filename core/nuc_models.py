"""Nucleotide substitution models: JC69, HKY, F84, and GTR.

Implements four-state reversible Markov models in ``A, C, G, T`` order. JC69,
HKY, and F84 evaluate transition probabilities in closed form; GTR builds Q,
scales it to mean rate 1 and symmetrises it as S = sqrt(pi) Q / sqrt(pi) before
eigendecomposition. Both ``set_matrix`` and ``set_vector`` write row-wise
cumulative transition probabilities so the evolution loop can sample states
without renormalising on the fly.
"""

import numpy as np

from ._validators import (validate_length, validate_nuc_freqs,
                          validate_output_array, validate_state)

NUM_NUC = 4

A, C, G, T = 0, 1, 2, 3


class JC69Model:
    """Jukes-Cantor 1969: equal frequencies, equal exchangeabilities.

    P_ii(t) = 1/4 + 3/4 exp(-4t/3);  P_ij(t) = 1/4 - 1/4 exp(-4t/3).
    """

    def __init__(self):
        self.num_states = NUM_NUC
        self.state_chars = "ACGT"
        pi = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)
        self.freq = pi
        self.add_freq = np.cumsum(pi)

    def set_matrix(self, matrix, length):
        validate_output_array(matrix, self.num_states * self.num_states)
        length = validate_length(length)
        m = matrix.reshape(self.num_states, self.num_states)
        a = np.exp(-4.0 * length / 3.0)
        p_same = 0.25 + 0.75 * a
        p_diff = 0.25 - 0.25 * a

        m[:] = p_diff
        for i in range(self.num_states):
            m[i, i] = p_same
        np.cumsum(m, axis=1, out=m)
        m[:, -1] = 1.0

    def set_vector(self, vector, state, length):
        # Runs once per site, so the array checks set_matrix does are skipped.
        state = validate_state(state, self.num_states)
        a = np.exp(-4.0 * length / 3.0)
        p_same = 0.25 + 0.75 * a
        p_diff = 0.25 - 0.25 * a

        v = vector[:self.num_states]
        v[:] = p_diff
        v[state] = p_same
        np.cumsum(v, out=v)
        v[-1] = 1.0


class HKYModel:
    """Hasegawa-Kishino-Yano 1985, closed form.

    ``kappa`` is back-solved from the requested transition/transversion ratio
    so that the *expected* ts/tv ratio equals ``tstv``.
    """

    def __init__(self, freqs, tstv, equal_exchangeability):
        validate_nuc_freqs(freqs)
        self.num_states = NUM_NUC
        self.state_chars = "ACGT"
        pi = np.array(freqs, dtype=np.float64)
        self.freq = pi
        self.add_freq = np.cumsum(pi)

        pi_A, pi_C, pi_G, pi_T = pi
        pi_R = pi_A + pi_G
        pi_Y = pi_C + pi_T
        pi_AG = pi_A * pi_G
        pi_CT = pi_C * pi_T

        if equal_exchangeability:
            kappa = 1.0
        else:
            if tstv is None or not np.isfinite(tstv) or tstv <= 0.0:
                raise ValueError(f"tstv must be finite and positive, got {tstv}")
            kappa = (tstv * pi_R * pi_Y) / (pi_AG + pi_CT)

        self._pi = (pi_A, pi_C, pi_G, pi_T)
        self._tabs = _hky_tabs(pi_A, pi_C, pi_G, pi_T, pi_R, pi_Y)
        self._beta = -1.0 / (2.0 * (pi_R * pi_Y + kappa * (pi_AG + pi_CT)))
        self._beta_A_R = self._beta * (1.0 + pi_R * (kappa - 1.0))
        self._beta_A_Y = self._beta * (1.0 + pi_Y * (kappa - 1.0))

    def set_matrix(self, matrix, length):
        validate_output_array(matrix, self.num_states * self.num_states)
        length = validate_length(length)
        aa = np.exp(self._beta * length)
        bbR = np.exp(self._beta_A_R * length)
        bbY = np.exp(self._beta_A_Y * length)
        _fill_matrix(matrix.reshape(-1), self._pi, self._tabs, aa, bbR, bbY)

    def set_vector(self, vector, state, length):
        state = validate_state(state, self.num_states)
        aa = np.exp(self._beta * length)
        bbR = np.exp(self._beta_A_R * length)
        bbY = np.exp(self._beta_A_Y * length)
        _fill_vector(vector, state, self._pi, self._tabs, aa, bbR, bbY)


class F84Model:
    """Felsenstein 1984, closed form.

    ``kappa`` follows the xi/xv parameterisation; a value below -min(pi_R,
    pi_Y) would make an instantaneous rate negative, so it is rejected.
    """

    def __init__(self, freqs, tstv, equal_exchangeability):
        validate_nuc_freqs(freqs)
        self.num_states = NUM_NUC
        self.state_chars = "ACGT"
        pi = np.array(freqs, dtype=np.float64)
        self.freq = pi
        self.add_freq = np.cumsum(pi)

        pi_A, pi_C, pi_G, pi_T = pi
        pi_R = pi_A + pi_G
        pi_Y = pi_C + pi_T
        pi_AG = pi_A * pi_G
        pi_CT = pi_C * pi_T
        sum_sq = pi_A**2 + pi_C**2 + pi_G**2 + pi_T**2
        sum_sq_grouped = (pi_A**2 / pi_R + pi_C**2 / pi_Y
                          + pi_G**2 / pi_R + pi_T**2 / pi_Y)

        if equal_exchangeability:
            kappa = 0.0
        else:
            if tstv is None or not np.isfinite(tstv) or tstv <= 0.0:
                raise ValueError(f"tstv must be finite and positive, got {tstv}")
            xi = pi_R * pi_Y * (pi_R * pi_Y * tstv - pi_AG - pi_CT)
            xv = pi_CT * pi_R + pi_AG * pi_Y
            kappa = xi / xv

        min_pi_ry = min(pi_R, pi_Y)
        if kappa < -min_pi_ry:
            raise ValueError(
                f"tstv={tstv} produces kappa={kappa:.6f} < "
                f"-min(pi_R, pi_Y)={-min_pi_ry:.6f}, resulting in negative "
                f"instantaneous substitution rates. Use a larger tstv ratio.")

        denom = (1.0 - sum_sq) + (kappa * (1.0 - sum_sq_grouped))
        if not np.isfinite(denom) or denom <= 0.0:
            raise ValueError(
                f"invalid F84 rate normalization denominator: {denom}. "
                f"Check tstv and base frequencies.")

        self._pi = (pi_A, pi_C, pi_G, pi_T)
        self._tabs = _hky_tabs(pi_A, pi_C, pi_G, pi_T, pi_R, pi_Y)
        self._mu = -1.0 / denom
        self._mu_kappa_1 = self._mu * (kappa + 1.0)

    def set_matrix(self, matrix, length):
        validate_output_array(matrix, self.num_states * self.num_states)
        length = validate_length(length)
        aa = np.exp(self._mu * length)
        bb = np.exp(self._mu_kappa_1 * length)
        _fill_matrix(matrix.reshape(-1), self._pi, self._tabs, aa, bb, bb)

    def set_vector(self, vector, state, length):
        state = validate_state(state, self.num_states)
        aa = np.exp(self._mu * length)
        bb = np.exp(self._mu_kappa_1 * length)
        _fill_vector(vector, state, self._pi, self._tabs, aa, bb, bb)


def _hky_tabs(pi_A, pi_C, pi_G, pi_T, pi_R, pi_Y):
    """The twelve closed-form coefficients HKY and F84 share."""
    return (
        pi_A * pi_Y / pi_R, pi_G / pi_R, pi_A / pi_R,
        pi_C * pi_R / pi_Y, pi_T / pi_Y, pi_C / pi_Y,
        pi_G * pi_Y / pi_R, pi_A / pi_R, pi_G / pi_R,
        pi_T * pi_R / pi_Y, pi_C / pi_Y, pi_T / pi_Y,
    )


class GTRModel:
    """General time reversible, 4x4 symmetric eigendecomposition.

    ``relative_rates`` are the exchangeabilities r_AC, r_AG, r_AT, r_CG, r_CT,
    r_GT; frequencies are in A, C, G, T order. Q is scaled to mean rate 1 and
    symmetrised as S = sqrt(pi) Q / sqrt(pi), so one eigh call covers all t.
    """

    def __init__(self, freqs, relative_rates):
        validate_nuc_freqs(freqs)
        _check_relative_rates(relative_rates, 6, "AC, AG, AT, CG, CT, GT")

        self.num_states = NUM_NUC
        self.state_chars = "ACGT"
        pi = np.array(freqs, dtype=np.float64)
        self.freq = pi
        self.add_freq = np.cumsum(pi)

        R = np.zeros((NUM_NUC, NUM_NUC))
        k = 0
        for i in range(NUM_NUC - 1):
            for j in range(i + 1, NUM_NUC):
                R[i][j] = R[j][i] = relative_rates[k]
                k += 1

        Q = R * pi[None, :]
        for i in range(NUM_NUC):
            Q[i][i] = -np.sum(Q[i])

        mean_rate = np.sum(pi * (-np.diag(Q)))
        if not np.isfinite(mean_rate) or mean_rate <= 0.0:
            raise ValueError(f"invalid GTR mean substitution rate: {mean_rate}")
        Q /= mean_rate

        sqrt_pi = np.sqrt(pi)
        S = sqrt_pi[:, None] * Q * (1.0 / sqrt_pi)[None, :]
        asym = np.max(np.abs(S - S.T))
        if asym > 1e-12:
            raise FloatingPointError(
                f"GTR reversible transform is unexpectedly asymmetric: {asym}")
        S = 0.5 * (S + S.T)

        self._roots, self._eigvecs = np.linalg.eigh(S)
        self._sqrt_pi = sqrt_pi

    def set_matrix(self, matrix, length):
        validate_output_array(matrix, self.num_states * self.num_states)
        length = validate_length(length)
        m = matrix.reshape(self.num_states, self.num_states)
        if length == 0.0:
            m[:] = np.eye(self.num_states)
        else:
            E = self._eigvecs * np.exp(self._roots * length)[None, :]
            E = E @ self._eigvecs.T
            m[:] = E * (self._sqrt_pi[None, :] / self._sqrt_pi[:, None])
        _clean_rows(m, "GTR transition matrix")

    def set_vector(self, vector, state, length):
        state = validate_state(state, self.num_states)
        v = vector[:self.num_states]
        if length == 0.0:
            v[:] = 0.0
            v[state] = 1.0
        else:
            E_row = self._eigvecs[state] * np.exp(self._roots * length)
            v[:] = E_row @ self._eigvecs.T
            v[:] *= self._sqrt_pi / self._sqrt_pi[state]
        _clean_row(v, "GTR transition vector")


def _check_relative_rates(rates, expected, order):
    if len(rates) != expected:
        raise ValueError(f"relative_rates must contain exactly {expected} "
                         f"values in order {order}, got {len(rates)}")
    for i, r in enumerate(rates):
        if not np.isfinite(r):
            raise ValueError(f"relative_rates must be finite, got rate[{i}]={r}")
        if r < 0.0:
            raise ValueError("relative_rates must be non-negative, "
                             f"got rate[{i}]={r}")
    if all(r == 0.0 for r in rates):
        raise ValueError("at least one relative rate must be positive")


def _clean_rows(m, what):
    """Clip tiny negative probabilities and renormalise each row."""
    min_prob = np.min(m)
    if min_prob < -1e-12:
        raise FloatingPointError(
            f"{what} contains significantly negative probability: {min_prob}")
    m[m < 0.0] = 0.0
    m /= m.sum(axis=1, keepdims=True)
    np.cumsum(m, axis=1, out=m)
    m[:, -1] = 1.0


def _clean_row(v, what):
    min_prob = np.min(v)
    if min_prob < -1e-12:
        raise FloatingPointError(
            f"{what} contains significantly negative probability: {min_prob}")
    v[v < 0.0] = 0.0
    total = v.sum()
    if total > 0.0:
        v /= total
    np.cumsum(v, out=v)
    v[-1] = 1.0


# --- HKY / F84 closed forms -------------------------------------------------
#
# For a state s the raw probabilities are pi_j (1 - exp) for a transversion and
# pi_j (1 + t1_j*exp + t2_j*exp_st) for a transition, where exp is exp(_beta*t)
# and exp_st is exp(_beta_A_R|Y*t) for a transition out of a purine/pyrimidine.
# pi[0..3] is A, C, G, T; tabs are the twelve coefficients from _hky_tabs.


def _fill_matrix(mat, pi, tabs, aa, bbR, bbY):
    """Fill a 4x4 raw HKY/F84 matrix and make it row-wise cumulative."""
    m = mat.reshape(4, 4)
    compl = 1.0 - aa

    m[0, 0] = pi[0] + tabs[0] * aa + tabs[1] * bbR
    m[0, 1] = pi[1] * compl
    m[0, 2] = pi[2] + tabs[6] * aa - tabs[8] * bbR
    m[0, 3] = pi[3] * compl

    m[1, 0] = pi[0] * compl
    m[1, 1] = pi[1] + tabs[3] * aa + tabs[4] * bbY
    m[1, 2] = pi[2] * compl
    m[1, 3] = pi[3] + tabs[9] * aa - tabs[11] * bbY

    m[2, 0] = pi[0] + tabs[0] * aa - tabs[2] * bbR
    m[2, 1] = pi[1] * compl
    m[2, 2] = pi[2] + tabs[6] * aa + tabs[7] * bbR
    m[2, 3] = pi[3] * compl

    m[3, 0] = pi[0] * compl
    m[3, 1] = pi[1] + tabs[3] * aa - tabs[5] * bbY
    m[3, 2] = pi[2] * compl
    m[3, 3] = pi[3] + tabs[9] * aa + tabs[10] * bbY

    np.cumsum(m, axis=1, out=m)
    m[:, 3] = 1.0


def _fill_vector(vec, state, pi, tabs, aa, bbR, bbY):
    """Fill the raw HKY/F84 row for *state*, then accumulate it."""
    v = vec[:4]
    compl = 1.0 - aa
    rows = (
        (pi[0] + tabs[0] * aa + tabs[1] * bbR, pi[1] * compl,
         pi[2] + tabs[6] * aa - tabs[8] * bbR, pi[3] * compl),
        (pi[0] * compl, pi[1] + tabs[3] * aa + tabs[4] * bbY,
         pi[2] * compl, pi[3] + tabs[9] * aa - tabs[11] * bbY),
        (pi[0] + tabs[0] * aa - tabs[2] * bbR, pi[1] * compl,
         pi[2] + tabs[6] * aa + tabs[7] * bbR, pi[3] * compl),
        (pi[0] * compl, pi[1] + tabs[3] * aa - tabs[5] * bbY,
         pi[2] * compl, pi[3] + tabs[9] * aa + tabs[10] * bbY),
    )
    if state < 0 or state > 3:
        raise ValueError(f"invalid nucleotide state: {state}")
    v[:] = rows[state]
    np.cumsum(v, out=v)
    v[3] = 1.0
