"""1D shallow-water solver: HLL + MUSCL/MC + forward Euler + dry tolerance.

Cartesian form:
    d_t h    + d_x (h u)              = 0
    d_t (h u)+ d_x (h u^2 + g h^2/2)  = -g h d_x h0

Axisymmetric (r-weighted) form:
    d_t (r h)    + d_r (r h u)                    = 0
    d_t (r h u)  + d_r (r h u^2 + (g/2) r h^2)    = (g/2) h^2 - g r h d_r h0

The r-weighting removes the 1/r singularity; u(0)=0 at the symmetry axis.
"""
import numpy as np


class SWE1D:
    def __init__(self, x, h0=None, g=9.81, geometry="cartesian",
                 riemann="hll", limiter="mc", cfl=0.3,
                 dry_tol=1e-8, h_min=1e-10):
        self.x = np.asarray(x, dtype=float)
        self.n = self.x.size
        self.dx = self.x[1] - self.x[0]
        self.g = g
        self.geometry = geometry
        self.riemann = riemann
        self.limiter = limiter
        self.cfl = cfl
        self.dry_tol = dry_tol
        self.h_min = h_min
        self.h0 = np.zeros(self.n) if h0 is None else np.asarray(h0, dtype=float)

        self.xf = np.empty(self.n + 1)
        self.xf[1:-1] = 0.5 * (self.x[:-1] + self.x[1:])
        self.xf[0] = self.x[0] - 0.5 * self.dx
        self.xf[-1] = self.x[-1] + 0.5 * self.dx
        self.r = self.x
        self.rf = self.xf

        self.h = np.zeros(self.n)
        self.u = np.zeros(self.n)
        self.t = 0.0
        self.steps = 0

    # ------------------------------------------------------------------ state
    def set_state(self, h, u):
        self.h = np.maximum(np.asarray(h, dtype=float).copy(), self.h_min)
        self.u = np.asarray(u, dtype=float).copy()
        self.u[self.h < self.dry_tol] = 0.0
        self.t = 0.0
        self.steps = 0

    # ------------------------------------------------------------------ ghost
    def _ghost(self, h, u, bc):
        hp = np.empty(self.n + 2)
        up = np.empty(self.n + 2)
        hp[1:-1] = h
        up[1:-1] = u

        if bc[0] == "transmissive":
            hp[0], up[0] = h[0], u[0]
        elif bc[0] in ("wall", "symmetry"):
            hp[0], up[0] = h[0], -u[0]
        else:
            raise ValueError(f"unknown left BC {bc[0]!r}")

        if bc[1] == "transmissive":
            hp[-1], up[-1] = h[-1], u[-1]
        elif bc[1] == "wall":
            hp[-1], up[-1] = h[-1], -u[-1]
        else:
            raise ValueError(f"unknown right BC {bc[1]!r}")
        return hp, up

    # ---------------------------------------------------------------- limiter
    @staticmethod
    def _minmod(a, b, c):
        return np.where((a > 0) & (b > 0) & (c > 0), np.minimum(np.minimum(a, b), c),
               np.where((a < 0) & (b < 0) & (c < 0), np.maximum(np.maximum(a, b), c), 0.0))

    def _slope(self, q):
        dqL = np.empty_like(q)
        dqR = np.empty_like(q)
        dqL[1:-1] = q[1:-1] - q[:-2]
        dqR[1:-1] = q[2:] - q[1:-1]
        dqL[0] = dqR[0] = q[1] - q[0]
        dqL[-1] = dqR[-1] = q[-1] - q[-2]
        dqC = 0.5 * (dqL + dqR)
        if self.limiter == "mc":
            return self._minmod(2.0 * dqL, dqC, 2.0 * dqR)
        if self.limiter == "minmod":
            return self._minmod(dqL, dqC, dqR)
        raise ValueError(f"unknown limiter {self.limiter!r}")

    # -------------------------------------------------------------- Riemann
    def _hll(self, hL, uL, hR, uR, rf):
        g = self.g
        cL = np.sqrt(g * hL)
        cR = np.sqrt(g * hR)
        SL = np.minimum(uL - cL, uR - cR)
        SR = np.maximum(uL + cL, uR + cR)

        FL0 = hL * uL
        FL1 = hL * uL * uL + 0.5 * g * hL * hL
        FR0 = hR * uR
        FR1 = hR * uR * uR + 0.5 * g * hR * hR
        UL0, UL1 = hL, hL * uL
        UR0, UR1 = hR, hR * uR

        F0 = np.where(SL >= 0, FL0,
             np.where(SR <= 0, FR0,
             (SR * FL0 - SL * FR0 + SL * SR * (UR0 - UL0)) / (SR - SL)))
        F1 = np.where(SL >= 0, FL1,
             np.where(SR <= 0, FR1,
             (SR * FL1 - SL * FR1 + SL * SR * (UR1 - UL1)) / (SR - SL)))
        return F0 * rf, F1 * rf

    def _hllc(self, hL, uL, hR, uR, rf):
        g = self.g
        cL = np.sqrt(g * hL)
        cR = np.sqrt(g * hR)
        SL = np.minimum(uL - cL, uR - cR)
        SR = np.maximum(uL + cL, uR + cR)
        pL = 0.5 * g * hL * hL
        pR = 0.5 * g * hR * hR
        den = hL * (SL - uL) - hR * (SR - uR)
        num = pR - pL + hL * uL * (SL - uL) - hR * uR * (SR - uR)
        SM = np.where(np.abs(den) > 1e-14, num / np.where(np.abs(den) > 1e-14, den, 1.0),
                      0.5 * (uL + uR))

        hLs = np.maximum(hL * (SL - uL) / np.where(np.abs(SL - SM) > 1e-14, SL - SM, 1.0), 1e-12)
        hRs = np.maximum(hR * (SR - uR) / np.where(np.abs(SR - SM) > 1e-14, SR - SM, 1.0), 1e-12)
        qL0, qL1 = hL, hL * uL
        qR0, qR1 = hR, hR * uR
        qLs0, qLs1 = hLs, hLs * SM
        qRs0, qRs1 = hRs, hRs * SM
        FL0, FL1 = hL * uL, hL * uL * uL + pL
        FR0, FR1 = hR * uR, hR * uR * uR + pR
        FLs0 = FL0 + SL * (qLs0 - qL0)
        FLs1 = FL1 + SL * (qLs1 - qL1)
        FRs0 = FR0 + SR * (qRs0 - qR0)
        FRs1 = FR1 + SR * (qRs1 - qR1)

        F0 = np.where(SL >= 0, FL0,
             np.where(SM >= 0, FLs0,
             np.where(SR > 0, FRs0, FR0)))
        F1 = np.where(SL >= 0, FL1,
             np.where(SM >= 0, FLs1,
             np.where(SR > 0, FRs1, FR1)))
        return F0 * rf, F1 * rf

    # -------------------------------------------------------------------- rhs
    def _rhs(self, h, u, bc):
        hp, up = self._ghost(h, u, bc)
        sh = self._slope(hp)
        su = self._slope(up)

        hL = np.maximum(hp[:-1] + 0.5 * sh[:-1], self.h_min)
        hR = np.maximum(hp[1:] - 0.5 * sh[1:], self.h_min)
        uL = up[:-1] + 0.5 * su[:-1]
        uR = up[1:] - 0.5 * su[1:]

        if self.geometry == "cartesian":
            rf = np.ones(self.n + 1)
        else:
            rf = self.rf

        if self.riemann == "hll":
            F0, F1 = self._hll(hL, uL, hR, uR, rf)
        elif self.riemann == "hllc":
            F0, F1 = self._hllc(hL, uL, hR, uR, rf)
        else:
            raise ValueError(f"unknown riemann {self.riemann!r}")

        dh0 = np.gradient(self.h0, self.x) if self.n > 1 else np.zeros(self.n)
        if self.geometry == "cartesian":
            S1 = -self.g * h * dh0
        else:
            S1 = 0.5 * self.g * h * h - self.g * self.r * h * dh0

        dq0 = -(F0[1:] - F0[:-1]) / self.dx
        dq1 = -(F1[1:] - F1[:-1]) / self.dx + S1
        return dq0, dq1

    # ------------------------------------------------------------------- step
    def step(self, dt, bc):
        q0 = self.h.copy() if self.geometry == "cartesian" else self.r * self.h
        q1 = q0 * self.u
        dq0, dq1 = self._rhs(self.h, self.u, bc)
        q0 = q0 + dt * dq0
        q1 = q1 + dt * dq1
        if self.geometry == "cartesian":
            h = q0
        else:
            h = q0 / self.r
        h = np.maximum(h, self.h_min)
        u = q1 / np.maximum(q0, 1e-30)
        if self.geometry == "radial":
            u[0] = 0.0
            q1[0] = 0.0
        u[h < self.dry_tol] = 0.0
        self.h = h
        self.u = u
        self.t += dt
        self.steps += 1

    # -------------------------------------------------------------------- run
    def run(self, t_end, bc):
        while self.t < t_end - 1e-12:
            c = np.sqrt(self.g * self.h)
            dt = self.cfl * self.dx / np.max(np.abs(self.u) + c)
            if self.t + dt > t_end:
                dt = t_end - self.t
            self.step(dt, bc)
        return self.h, self.u
