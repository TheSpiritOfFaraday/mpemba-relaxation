"""Constrained two-basin Arrhenius model. Discovery calculations, not certificates.

All energies and barriers use k_B T_b = 1; nu_0 = 1.
States 0,3 form basin A and states 1,2 basin B. Two intrabasin
barrier excesses and one shared interbasin excess are chosen before search.
No independently fitted edge rates, prefactors, or spectral modes.
"""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

PAIRS = [(i, j) for i in range(4) for j in range(i + 1, 4)]


def construct(parameters):
    """parameters = three gaps, h_A, h_B, h_X, T_c, T_h."""
    gaps = np.asarray(parameters[:3])
    E = np.r_[0., np.cumsum(gaps)]
    hA, hB, hX, Tc, Th = parameters[3:]
    M = np.zeros((4, 4))
    B = np.zeros((4, 4))
    for i, j in PAIRS:
        h = hA if (i, j) == (0, 3) else hB if (i, j) == (1, 2) else hX
        B[i, j] = B[j, i] = E[j] + h
        M[i, j] = np.exp(-h)
        M[j, i] = np.exp(-h - (E[j] - E[i]))
    M -= np.diag(M.sum(axis=0))
    return E, B, M


def thermal(E, T):
    q = np.exp(-E / T)
    return q / q.sum()


class Evolution:
    def __init__(self, parameters):
        self.parameters = np.asarray(parameters)
        self.E, self.B, self.M = construct(parameters)
        self.pi = thermal(self.E, 1)
        root = np.sqrt(self.pi)
        S = self.M * root[None, :] / root[:, None]
        self.lam, self.U = eigh(S)
        self.tau = -1 / self.lam[-2]
        self.lam[-1] = 0
        self.root = root
        self.initial = np.array([thermal(self.E, t) for t in parameters[-2:]])
        self.coeff = (self.initial / root) @ self.U

    def probabilities(self, s):
        s = np.atleast_1d(s)
        # s is time in slow-relaxation units.
        return np.einsum('ik,tk,jk,j->itj', self.coeff,
                         np.exp(s[:, None] * self.tau * self.lam),
                         self.U, self.root)

    def distances(self, s, metric="kl"):
        p = self.probabilities(s)
        if metric == "kl":
            # Stable KL near equilibrium; sum of deviations is zero analytically.
            z = (p - self.pi) / self.pi
            return np.sum(self.pi * ((1 + z) * np.log1p(z) - z), axis=-1)
        if metric == "tv":
            return .5 * np.sum(np.abs(p - self.pi), axis=-1)
        return np.sum((p - self.pi)**2 / self.pi, axis=-1)

    def hitting(self, delta, preparation):
        return brentq(lambda s: self.distances(s)[preparation, 0] - delta,
                      0, 20, xtol=1e-12)

    def metrics(self):
        deltas = np.geomspace(.001, .01, 31)
        times = np.array([[self.hitting(d, j) for d in deltas] for j in range(2)])
        s = np.linspace(.01, 5, 1000)
        D = self.distances(s)
        valid = (D[1] <= .9 * D[0]) & (D[0] >= .001)
        ids = np.flatnonzero(valid)
        blocks = np.split(ids, np.flatnonzero(np.diff(ids) != 1)+1)
        block = max(blocks, key=len)
        window = None if not len(block) else [float(s[block[0]]), float(s[block[-1]])]
        off = self.M[~np.eye(4, dtype=bool)]
        return dict(parameters=self.parameters.tolist(), energies=self.E.tolist(),
                    initial_KL=self.distances(0)[:, 0].tolist(),
                    eigenvalues=self.lam.tolist(), tau=float(self.tau),
                    rate_ratio=float(max(off)/min(off)),
                    min_time_saving_sampled=float(min(1-times[1]/times[0])),
                    distance_window_sampled=window,
                    thresholds=deltas.tolist(), hitting_times_tau=times.tolist())
