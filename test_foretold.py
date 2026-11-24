# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Unit tests for foretold-lens (stdlib only, no network).

Made with AI agents working in MV2's name (MV2, Munim, Inc.).
"""
import math
import random
import unittest

import events as E
import lensfit as L
import phot as P
import ztf_io as Z


def synth(nx, ny, stars, bkg=100.0, noise=0.0, seed=0, fwhm=2.0, hdr=None):
    sig = fwhm / 2.3548
    norm = 1.0 / (2 * math.pi * sig * sig)
    data = [bkg] * (nx * ny)
    for x0, y0, flux in stars:
        for j in range(max(0, int(y0) - 8), min(ny, int(y0) + 9)):
            for i in range(max(0, int(x0) - 8), min(nx, int(x0) + 9)):
                data[j * nx + i] += flux * norm * math.exp(-((i - x0) ** 2 + (j - y0) ** 2) / (2 * sig * sig))
    if noise:
        rng = random.Random(seed)
        data = [v + rng.gauss(0.0, noise) for v in data]
    h = {"SEEING": fwhm, "GAIN": 6.1, "SATURATE": 1.0e6}
    if hdr:
        h.update(hdr)
    return P.Image(h, nx, ny, data)


class TestMagnification(unittest.TestCase):
    def test_values(self):
        self.assertAlmostEqual(L.magnification(1.0), 3.0 / math.sqrt(5.0), places=9)
        self.assertAlmostEqual(L.magnification(0.5), 2.25 / (0.5 * math.sqrt(4.25)), places=9)
        self.assertGreater(L.magnification(0.1), 9.9)


class TestEvents(unittest.TestCase):
    def test_julian_year(self):
        self.assertAlmostEqual(E.jyear_to_mjd(2000.0), 51544.5, places=6)
        self.assertAlmostEqual(E.jyear_to_mjd(2025.5475410761549), 60875.7394, places=3)

    def test_selection_rule(self):
        self.assertEqual(len(E.EVENTS), 6)
        for ev in E.EVENTS:
            self.assertGreater(ev["magnification"], 0.003)
            self.assertGreater(ev["ob_dec"], -28.0)
            self.assertTrue(2018.21 <= ev["tca"] <= 2025.80)

    def test_headline_derived(self):
        ev = [e for e in E.EVENTS if e["event_id"] == E.HEADLINE][0]
        d = E.derived(ev)
        self.assertAlmostEqual(d["u0"], 0.86025, places=4)
        self.assertTrue(250.0 < d["mu_rel"] < 259.0)
        self.assertTrue(9.0 < d["te_days"] < 9.3)
        self.assertTrue(0.095 < d["fs"]["g"] < 0.099)
        self.assertTrue(0.030 < d["fs"]["r"] < 0.034)
        self.assertEqual(d["window_half"], 10.0)
        self.assertTrue(46.0 < d["exclusion"] < 47.0)
        self.assertLess(abs(d["center_ra"] - 286.4309), 3e-4)
        self.assertLess(abs(d["center_dec"] - 31.1761), 3e-4)
        # Kluter's u equals dist / thetaE for this event
        self.assertAlmostEqual(ev["dist"] / ev["theta_e"], ev["u"], places=3)


class TestFits(unittest.TestCase):
    def test_roundtrip(self):
        data = [0.5 * k for k in range(35)]
        blob = P.write_fits([("GAIN", 6.1), ("SEEING", 2.0), ("SATURATE", 50000.0), ("CTYPE1", "RA---TPV")],
                            7, 5, data)
        hdr, nx, ny, arr = P.read_fits(blob)
        self.assertEqual((nx, ny), (7, 5))
        self.assertEqual(arr[3 * 7 + 2], 11.5)
        self.assertEqual(hdr["GAIN"], 6.1)
        self.assertEqual(hdr["CTYPE1"], "RA---TPV")
        self.assertEqual(hdr["BITPIX"], -32)

    def test_not_fits(self):
        with self.assertRaises(ValueError):
            P.read_fits(b"<html>error</html>")


class TestPhot(unittest.TestCase):
    def test_aperture_and_centroid(self):
        img = synth(61, 61, [(30.3, 29.8, 10000.0)])
        m = P.measure_star(img, 30.0, 30.0, 5.0, 6.1, 3.0, 3.0)
        self.assertIsNotNone(m)
        self.assertLess(abs(m["x"] - 30.3), 0.1)
        self.assertLess(abs(m["y"] - 29.8), 0.1)
        self.assertLess(abs(m["flux"] - 10000.0), 50.0)

    def test_find_stars(self):
        stars = [(50, 50, 20000.0), (30, 70, 15000.0), (75, 30, 8000.0), (76, 35, 8000.0), (10, 50, 9000.0)]
        img = synth(101, 101, stars, noise=5.0, seed=1)
        got = [(i, j) for _, i, j in P.find_stars(img)]
        self.assertEqual(got, [(50, 50), (30, 70)])
        got = [(i, j) for _, i, j in P.find_stars(img, sat=6000.0)]
        self.assertEqual(got, [(30, 70)])

    def test_measure_group(self):
        comps = [(-30, 25, 8000.0), (28, -27, 12000.0), (35, 30, 6000.0), (-32, -30, 10000.0)]
        fac = [1.0] * 12
        fac[5] = 1.2
        items = []
        for k in range(12):
            tr = 0.9 + 0.02 * k  # transparency changes cancel in the ratio
            stars = [(60, 60, 5000.0 * fac[k] * tr)] + [(60 + dx, 60 + dy, f * tr) for dx, dy, f in comps]
            img = synth(121, 121, stars, noise=3.0, seed=10 + k)
            items.append(dict(mjd=60000.0 + k, img=img, meta=dict(group="zg_test", band="g")))
        out, diag = P.measure_group(items)
        self.assertEqual(len(out), 12, diag)
        med = P.median([o["rel"] for o in out])
        for o in out:
            k = int(round(o["mjd"] - 60000.0))
            self.assertLess(abs(o["rel"] / med - fac[k]), 0.01)
            self.assertGreater(o["rel_err"], 0.0)


def make_series(n, tc=None, u0=0.5, te=6.0, amp=None, noise=0.005, seed=3):
    rng = random.Random(seed)
    t, y, s, g = [], [], [], []
    for k in range(n):
        for j, b in ((0, "g"), (1, "r")):
            tt = k + 0.3 * j
            v = 1.0 + (rng.gauss(0.0, noise) if noise else 0.0)
            if tc is not None:
                v += amp[b] * L.model_m(tt, tc, u0, te)
            t.append(tt)
            y.append(v)
            s.append(noise if noise else 0.005)
            g.append(j)
    return L.Series(t, y, s, g, {0: "g", 1: "r"})


class TestFit(unittest.TestCase):
    def test_closed_form(self):
        ser = make_series(40, tc=20.2, amp={"g": 0.08, "r": 0.03})
        D, amp = L.delta_chi2(ser, 20.0, 0.6, 5.0)

        def chi2(ag, ar):
            tot = 0.0
            for j, a in ((0, ag), (1, ar)):
                idx = [i for i in range(len(ser.t)) if ser.g[i] == j]
                m = [L.model_m(ser.t[i], 20.0, 0.6, 5.0) for i in idx]
                sw = sum(ser.w[i] for i in idx)
                c = sum(ser.w[i] * (ser.y[i] - a * mm) for i, mm in zip(idx, m)) / sw
                tot += sum(ser.w[i] * (ser.y[i] - c - a * mm) ** 2 for i, mm in zip(idx, m))
            return tot

        flat = chi2(0.0, 0.0)
        best = chi2(amp["g"], amp["r"])
        self.assertAlmostEqual(flat - best, D, delta=1e-6 * max(1.0, D))
        for dg, dr in ((0.01, 0), (-0.01, 0), (0, 0.01), (0, -0.01)):
            self.assertGreater(chi2(amp["g"] + dg, amp["r"] + dr), best)

    def test_scan_finds_event(self):
        ser = make_series(200, tc=100.3, amp={"g": 0.1, "r": 0.05})
        best = L.scan(ser, 100.0, 10.0, te_grid=[3.0, 6.0, 12.0], u0_grid=[0.3, 0.5, 1.0], t0_step=0.5)
        self.assertGreater(best["D"], 100.0)
        self.assertLess(abs(best["t0"] - 100.3), 1.0)

    def test_slide_centers(self):
        c = L.slide_centers(200.0, 10.0, 40.0, 0.0, 399.0)
        self.assertIn(160.0, c)
        self.assertIn(240.0, c)
        self.assertNotIn(170.0, c)
        self.assertEqual(c, sorted(c))

    def _grade(self, ser):
        return L.grade(ser, 200.0, 10.0, 40.0, 6.0, 0.5, 0.8, {"g": 0.1, "r": 0.05}, 2.0,
                       te_grid=[3.0, 6.0, 12.0], u0_grid=[0.3, 0.5, 1.0], t0_step=1.0)

    def test_grade_detected(self):
        ser = make_series(400, tc=200.4, amp={"g": 0.1, "r": 0.05})
        r = self._grade(ser)
        self.assertEqual(r["verdict"], "DETECTED")
        self.assertGreaterEqual(r["n_slides"], 19)
        self.assertLessEqual(abs(r["t0_offset_days"]), 1.5)
        self.assertTrue(r["chromatic_g_ge_r"])

    def test_grade_not_seen(self):
        ser = make_series(400, noise=0.0)
        r = self._grade(ser)
        self.assertEqual(r["verdict"], "NOT_SEEN_PREDICTION_EXCLUDED")
        self.assertEqual(r["recovery_median"], 1.0)

    def test_no_coverage(self):
        t = [float(k) for k in range(100)]
        ser = L.Series(t, [1.0] * 100, [0.01] * 100, [0] * 100, {0: "g"})
        r = L.grade(ser, 300.0, 10.0, 40.0, 6.0, 0.5, 0.8, {"g": 0.1}, 2.0,
                    te_grid=[6.0], u0_grid=[0.5], t0_step=1.0)
        self.assertEqual(r["verdict"], "NO_COVERAGE")

    def test_mass_profile(self):
        mu, th, d0 = 250.0, 6.0, 4.8
        te = th / mu * 365.25
        t, y, s, g = [], [], [], []
        for k in range(121):
            for j, b, a in ((0, "g", 0.1), (1, "r", 0.04)):
                tt = float(k) + 0.25 * j
                t.append(tt)
                y.append(1.0 + a * L.model_m(tt, 60.0, d0 / th, te))
                s.append(0.003)
                g.append(j)
        ser = L.Series(t, y, s, g, {0: "g", 1: "r"})
        r = L.mass_profile(ser, 60.0, mu, 26.7, 4.8, 1.5, theta_grid=[4.0, 5.0, 6.0, 7.0, 8.0],
                           t0_half=1.0, t0_step=0.5)
        self.assertEqual(r["theta_best"], 6.0)
        self.assertLess(abs(r["mass_best"] - 36.0 / (8.144 * 26.7)), 1e-9)

    def test_normalise_and_clip(self):
        pts = []
        for k in range(30):
            rel = 2.0
            if k == 10:
                rel = 2.5
            if 18 <= k <= 22:
                rel = 2.16
            pts.append(dict(mjd=60000.0 + k, rel=rel, rel_err=0.02, group="zg_a", band="g"))
        out = L.normalise_and_clip(pts, lambda m: True)
        flags = {int(round(p["mjd"] - 60000.0)): p["clipped"] for p in out}
        self.assertTrue(flags[10])
        for k in range(18, 23):
            self.assertFalse(flags[k])
        self.assertEqual(sum(flags.values()), 1)
        self.assertEqual(out[0]["err_scale"], 1.0)


class TestZtfIo(unittest.TestCase):
    def test_cutout_url(self):
        row = dict(filefracday="20250720301447", field="684", filtercode="zg", ccdid="1", qid="1",
                   imgtypecode="o")
        u = Z.cutout_url(row, 286.4309, 31.1761)
        self.assertTrue(u.startswith(Z.IBE_DATA + "/2025/0720/301447/ztf_20250720301447_000684_zg_c01_o_q1_sciimg.fits?"))
        self.assertIn("size=180arcsec", u)
        self.assertIn("gzip=false", u)
        self.assertIn("center=286.430900%2C31.176100", u)
        self.assertEqual(Z.group_key(row), "zg_f684_c1_q1")

    def test_quality(self):
        ok = dict(imgtypecode="o", infobits="0", seeing="2.0", maglimit="20.5")
        self.assertTrue(Z.passes_quality(ok))
        self.assertFalse(Z.passes_quality(dict(ok, infobits="8")))
        self.assertFalse(Z.passes_quality(dict(ok, seeing="4.5")))
        self.assertFalse(Z.passes_quality(dict(ok, maglimit="19.0")))
        self.assertFalse(Z.passes_quality(dict(ok, seeing="")))


if __name__ == "__main__":
    unittest.main()
