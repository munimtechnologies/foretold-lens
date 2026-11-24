# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Point-lens detection test, time-slide calibration, injection-recovery and frozen verdicts.

Model for the normalised flux y of a blend (lens + source) in group j (one ZTF field/CCD/quadrant/filter)
with band b:  y = c_j + a_b * (A(u(t)) - 1),  u(t) = sqrt(u0^2 + ((t - t0)/tE)^2),  a_b >= 0,
where a_b is the source's share of the blended flux in band b. For fixed (t0, u0, tE) the best c_j and
a_b follow in closed form, so a detection statistic D = chi2_flat - chi2_lens is cheap.

stdlib only. Made with AI agents working in MV2's name (MV2, Munim, Inc.).
"""
import bisect
import math

TE_GRID = [2.0 * 1.3 ** k for k in range(12)]          # 2 .. 35.8 days
U0_GRID = [0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.3, 1.7, 2.2, 3.0]
T0_STEP = 0.5
U_MAX = 8.0
D_MIN = 25.0          # detection floor on D
D_MARGINAL = 9.0
N_SLIDES_DET = 19     # time-slides needed to call a detection
N_SLIDES_REC = 10     # time-slides needed for injection-recovery verdicts
MIN_POINTS = 8
MIN_SIDE = 2
CLIP_NSIG = 5.0
CLIP_HALF = 2.0       # days, neighbour window for the flare/outlier clip


def magnification(u):
    u2 = u * u
    return (u2 + 2.0) / (u * math.sqrt(u2 + 4.0))


class Series:
    """Sorted light curve. y normalised per group, w = 1/sigma^2, grp index, band per group."""

    def __init__(self, t, y, s, grp, band_of_group):
        order = sorted(range(len(t)), key=lambda i: t[i])
        self.t = [t[i] for i in order]
        self.y = [y[i] for i in order]
        self.s = [s[i] for i in order]
        self.w = [1.0 / (s[i] * s[i]) for i in order]
        self.g = [grp[i] for i in order]
        self.band_of_group = dict(band_of_group)
        self.groups = sorted(set(self.g))
        self._flat()

    def _flat(self):
        self.Sw = {j: 0.0 for j in self.groups}
        self.Swy = {j: 0.0 for j in self.groups}
        for w, y, j in zip(self.w, self.y, self.g):
            self.Sw[j] += w
            self.Swy[j] += w * y

    def with_y(self, y):
        s = Series.__new__(Series)
        s.t, s.y, s.s, s.w, s.g = self.t, list(y), self.s, self.w, self.g
        s.band_of_group, s.groups = self.band_of_group, self.groups
        s._flat()
        return s

    def count(self, lo, hi):
        return bisect.bisect_right(self.t, hi) - bisect.bisect_left(self.t, lo)


def model_m(t, t0, u0, te):
    u = math.sqrt(u0 * u0 + ((t - t0) / te) ** 2)
    return magnification(u) - 1.0


def delta_chi2(ser, t0, u0, te):
    """Return (D, {band: a}) for the point-lens model at (t0, u0, tE)."""
    span = te * math.sqrt(max(U_MAX * U_MAX - u0 * u0, 0.0))
    lo = bisect.bisect_left(ser.t, t0 - span)
    hi = bisect.bisect_right(ser.t, t0 + span)
    Swm, Swmm, Swmy = {}, {}, {}
    for i in range(lo, hi):
        u = math.sqrt(u0 * u0 + ((ser.t[i] - t0) / te) ** 2)
        u2 = u * u
        m = (u2 + 2.0) / (u * math.sqrt(u2 + 4.0)) - 1.0
        j = ser.g[i]
        w = ser.w[i]
        wm = w * m
        Swm[j] = Swm.get(j, 0.0) + wm
        Swmm[j] = Swmm.get(j, 0.0) + wm * m
        Swmy[j] = Swmy.get(j, 0.0) + wm * ser.y[i]
    N, Dn = {}, {}
    for j in Swm:
        b = ser.band_of_group[j]
        N[b] = N.get(b, 0.0) + Swmy[j] - Swm[j] * ser.Swy[j] / ser.Sw[j]
        Dn[b] = Dn.get(b, 0.0) + Swmm[j] - Swm[j] * Swm[j] / ser.Sw[j]
    tot, amp = 0.0, {}
    for b in N:
        if Dn[b] > 1e-12 and N[b] > 0:
            amp[b] = N[b] / Dn[b]
            tot += N[b] * N[b] / Dn[b]
        else:
            amp[b] = 0.0
    return tot, amp


def frange(lo, hi, step):
    n = int(math.floor((hi - lo) / step + 1e-9))
    return [lo + k * step for k in range(n + 1)]


def scan(ser, c, w, te_grid=TE_GRID, u0_grid=U0_GRID, t0_step=T0_STEP):
    """Max D over t0 in [c-w, c+w] and the grids. Returns dict(D, t0, u0, tE, amp)."""
    best = dict(D=0.0, t0=None, u0=None, tE=None, amp={})
    for t0 in frange(c - w, c + w, t0_step):
        for te in te_grid:
            for u0 in u0_grid:
                d, amp = delta_chi2(ser, t0, u0, te)
                if d > best["D"]:
                    best = dict(D=d, t0=t0, u0=u0, tE=te, amp=amp)
    return best


def t0_profile(ser, c, w, te_grid=TE_GRID, u0_grid=U0_GRID, t0_step=T0_STEP):
    out = []
    for t0 in frange(c - w, c + w, t0_step):
        bd = 0.0
        for te in te_grid:
            for u0 in u0_grid:
                d, _ = delta_chi2(ser, t0, u0, te)
                bd = max(bd, d)
        out.append((t0, bd))
    return out


def coverage_ok(ser, c, w, te):
    half = w + te
    n = ser.count(c - half, c + half)
    left = ser.count(c - half, c - 1e-9)
    right = ser.count(c + 1e-9, c + half)
    return n >= MIN_POINTS and left >= MIN_SIDE and right >= MIN_SIDE


def inject(ser, t0, u0, te, amp_by_band):
    y = list(ser.y)
    span = te * U_MAX * 1.5
    for i in range(len(y)):
        if abs(ser.t[i] - t0) <= span:
            y[i] += amp_by_band.get(ser.band_of_group[ser.g[i]], 0.0) * model_m(ser.t[i], t0, u0, te)
    return ser.with_y(y)


def slide_centers(tc, w, exclusion, t_min, t_max):
    step = max(7.0, w)
    out = []
    k = 1
    while True:
        moved = False
        for c in (tc - k * step, tc + k * step):
            if t_min - w <= c <= t_max + w:
                moved = True
                if abs(c - tc) >= exclusion:
                    out.append(c)
        if not moved:
            break
        k += 1
    return sorted(out)


def quantile_95(vals):
    s = sorted(vals)
    if not s:
        return float("nan")
    return s[max(0, int(math.ceil(0.95 * len(s))) - 1)]


def grade(ser, tc, w, exclusion, te_pred, u0_pred, u0_pess, amp_pred, err_tca_days,
          te_grid=TE_GRID, u0_grid=U0_GRID, t0_step=T0_STEP):
    """Frozen verdict for one event. All times in MJD."""
    res = dict(tc=tc, window_half=w, exclusion=exclusion, te_pred=te_pred, u0_pred=u0_pred,
               u0_pess=u0_pess, amp_pred=amp_pred, n_points=len(ser.t))
    if not ser.t or not coverage_ok(ser, tc, w, te_pred):
        res.update(verdict="NO_COVERAGE",
                   n_near=ser.count(tc - w - te_pred, tc + w - 0 + te_pred) if ser.t else 0)
        return res
    obs = scan(ser, tc, w, te_grid, u0_grid, t0_step)
    res["observed"] = obs
    centers = [c for c in slide_centers(tc, w, exclusion, ser.t[0], ser.t[-1])
               if coverage_ok(ser, c, w, te_pred)]
    nulls, inj_med, inj_pess = [], [], []
    for c in centers:
        nulls.append(scan(ser, c, w, te_grid, u0_grid, t0_step)["D"])
    for c in centers:
        inj_med.append(scan(inject(ser, c, u0_pred, te_pred, amp_pred), c, w, te_grid, u0_grid, t0_step)["D"])
        inj_pess.append(scan(inject(ser, c, u0_pess, te_pred, amp_pred), c, w, te_grid, u0_grid, t0_step)["D"])
    n = len(centers)
    thr = max(D_MIN, quantile_95(nulls)) if n else D_MIN
    p = (1 + sum(1 for d in nulls if d >= obs["D"])) / float(1 + n)
    rec_med = (sum(1 for d in inj_med if d >= thr) / float(n)) if n else None
    rec_pess = (sum(1 for d in inj_pess if d >= thr) / float(n)) if n else None
    res.update(slide_centers=centers, null_D=nulls, inj_D_median=inj_med, inj_D_pess=inj_pess,
               n_slides=n, threshold=thr, p_slide=p, recovery_median=rec_med, recovery_pess=rec_pess)
    D = obs["D"]
    if n >= N_SLIDES_DET and D >= thr:
        v = "DETECTED"
    elif D >= D_MIN:
        v = "CANDIDATE"
    elif D >= D_MARGINAL:
        v = "MARGINAL"
    elif n < N_SLIDES_REC:
        v = "NOT_SEEN_UNCALIBRATED"
    elif rec_med >= 0.9:
        v = "NOT_SEEN_PREDICTION_EXCLUDED"
    elif rec_med >= 0.5:
        v = "NOT_SEEN_PREDICTION_DISFAVOURED"
    else:
        v = "NOT_SEEN_INSENSITIVE"
    res["verdict"] = v
    if obs["t0"] is not None:
        res["t0_offset_days"] = obs["t0"] - tc
        res["t0_offset_sigma"] = (obs["t0"] - tc) / err_tca_days if err_tca_days > 0 else None
        prof = t0_profile(ser, tc, w, te_grid, u0_grid, t0_step)
        best = max(d for _, d in prof)
        inside = [t for t, d in prof if d >= best - 1.0]
        res["t0_interval"] = [min(inside), max(inside)] if inside else None
        a = obs["amp"]
        if "g" in a and "r" in a:
            res["chromatic_g_ge_r"] = a["g"] >= a["r"]
    return res


def mass_profile(ser, t0_best, mu_rel, pi_rel, dist, err_dist, kappa=8.144,
                 theta_grid=None, t0_half=3.0, t0_step=0.25):
    """Exploratory: geometry-constrained fit (u0 = d0/thetaE, tE = thetaE/mu_rel) with a Gaussian
    prior on the Gaia impact parameter d0. Returns profile [(thetaE, chi2_min)], 68% interval, mass."""
    if theta_grid is None:
        theta_grid = [0.5 + 0.25 * k for k in range(79)]
    d_lo = max(0.05, dist - 3 * err_dist)
    d_hi = dist + 3 * err_dist
    d_grid = frange(d_lo, d_hi, max((d_hi - d_lo) / 40.0, 0.05))
    prof = []
    for th in theta_grid:
        te = th / mu_rel * 365.25
        best = None
        for d0 in d_grid:
            u0 = d0 / th
            if u0 > 6.0:
                continue
            prior = ((d0 - dist) / err_dist) ** 2
            for t0 in frange(t0_best - t0_half, t0_best + t0_half, t0_step):
                D, _ = delta_chi2(ser, t0, u0, te)
                val = prior - D
                if best is None or val < best:
                    best = val
        if best is not None:
            prof.append((th, best))
    if not prof:
        return dict(profile=[], interval=None)
    m = min(v for _, v in prof)
    inside = [th for th, v in prof if v <= m + 1.0]
    th_best = min(prof, key=lambda p: p[1])[0]

    def mass(th):
        return th * th / (kappa * pi_rel)

    return dict(profile=prof, theta_best=th_best, theta_interval=[min(inside), max(inside)],
                mass_best=mass(th_best), mass_interval=[mass(min(inside)), mass(max(inside))],
                at_grid_edge=(min(inside) <= theta_grid[0] or max(inside) >= theta_grid[-1]))


def normalise_and_clip(points, baseline_ok):
    """points: list of dicts(mjd, rel, rel_err, group, band). Normalise each group to the median of its
    baseline points (all points if none), rescale errors by the robust baseline scatter
    1.4826*median(|y-1|/sigma) (never below 1), then flag 5-sigma outliers relative to the median of
    >= 3 same-band neighbours within +-2 d (single-epoch flares/artefacts; multi-day changes survive)."""
    from phot import median
    by = {}
    for p in points:
        by.setdefault(p["group"], []).append(p)
    out = []
    for gkey, ps in by.items():
        base = [p for p in ps if baseline_ok(p["mjd"])] or ps
        m = median([p["rel"] for p in base])
        if not m > 0:
            continue
        for p in ps:
            p["y"] = p["rel"] / m
            p["s"] = p["rel_err"] / m
        bb = [p for p in ps if baseline_ok(p["mjd"])] or ps
        z = [abs(p["y"] - 1.0) / p["s"] for p in bb if p["s"] > 0]
        scale = max(1.0, 1.4826 * median(z)) if z else 1.0   # robust: flares cannot inflate it
        for p in ps:
            p["s"] *= scale
            p["err_scale"] = scale
        out.extend(ps)
    out.sort(key=lambda p: p["mjd"])
    for p in out:
        nb = [q["y"] for q in out if q is not p and q["band"] == p["band"]
              and abs(q["mjd"] - p["mjd"]) <= CLIP_HALF]
        p["clipped"] = False
        if len(nb) >= 3:
            med = median(nb)
            if abs(p["y"] - med) > CLIP_NSIG * p["s"]:
                p["clipped"] = True
    return out
