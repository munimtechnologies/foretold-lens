# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""foretold-lens job: did the Gaia-predicted photometric microlensing events happen on schedule?

For each frozen event (events.py): find every public ZTF science image of the position from
tca-420 d to min(tca+150 d, DR24 cut-off), apply the frozen quality cuts, fetch 180" cut-outs,
measure forced differential aperture photometry (phot.py) and grade with the frozen point-lens
test (lensfit.py). Writes results.json, per-event light curves, model curves and summary.md.

  python job.py --data-dir /data --out-dir /out

Made with AI agents working in MV2's name (MV2, Munim, Inc.).
"""
import argparse
import csv
import json
import os
import sys
import time
import traceback
from collections import Counter

import events as E
import lensfit as L
import phot as P
import ztf_io as Z

VERSION = "0.1.0"
PRE_DAYS, POST_DAYS = 420.0, 150.0
MAX_IMAGES = 1500


def log(*a):
    print("[%s]" % time.strftime("%H:%M:%S"), *a, flush=True)


def process_event(ev, data_dir, out_dir, max_images=MAX_IMAGES):
    d = E.derived(ev)
    tc = d["tca_mjd"]
    res = dict(event_id=ev["event_id"], name=ev["name"], derived=d, catalogue=ev)
    lo = tc - PRE_DAYS
    hi = min(tc + POST_DAYS, E.DR24_END_MJD)
    res["data_window_mjd"] = [lo, hi]
    if hi <= lo:
        res["verdict"] = "NO_COVERAGE"
        return res
    rows = Z.search(d["center_ra"], d["center_dec"], lo + 2400000.5, hi + 2400000.5, log=log)
    res["n_images_found"] = len(rows)
    good = [r for r in rows if Z.passes_quality(r)]
    res["n_images_quality"] = len(good)
    good.sort(key=lambda r: abs(float(r["obsjd"]) - 2400000.5 - tc))
    good = good[:max_images]
    log(ev["event_id"], "images found", len(rows), "after quality cuts", len(good))
    cache = os.path.join(data_dir, "cutouts", ev["event_id"])
    items, failed = [], 0
    for k, r in enumerate(good):
        try:
            blob = Z.get_cutout(r, d["center_ra"], d["center_dec"], cache, log=log)
            hdr, nx, ny, data = P.read_fits(blob)
            items.append(dict(mjd=float(r["obsjd"]) - 2400000.5, img=P.Image(hdr, nx, ny, data),
                              meta=dict(group=Z.group_key(r), band=r["filtercode"].strip()[-1],
                                        filefracday=r["filefracday"], seeing=float(r["seeing"]),
                                        maglimit=float(r["maglimit"]))))
        except Exception as e:
            failed += 1
            log("cut-out failed:", str(e)[:200])
        if (k + 1) % 50 == 0:
            log("  %d/%d cut-outs" % (k + 1, len(good)))
    res["n_cutouts_ok"] = len(items)
    res["n_cutouts_failed"] = failed
    if items:
        shape = Counter((it["img"].nx, it["img"].ny) for it in items).most_common(1)[0][0]
        before = len(items)
        items = [it for it in items if (it["img"].nx, it["img"].ny) == shape]
        res["cutout_shape"] = list(shape)
        res["n_cutouts_edge_dropped"] = before - len(items)
    groups = {}
    for it in items:
        groups.setdefault(it["meta"]["group"], []).append(it)
    points, gdiag = [], {}
    for gk in sorted(groups):
        out, diag = P.measure_group(groups[gk], log=log)
        gdiag[gk] = diag
        for o in out:
            points.append(dict(mjd=o["mjd"], rel=o["rel"], rel_err=o["rel_err"], group=gk,
                               band=o["meta"]["band"], centroid_ok=o["centroid_ok"],
                               dx=o["dx"], dy=o["dy"], seeing=o["meta"]["seeing"]))
    res["groups"] = gdiag
    for it in items:  # free memory
        it["img"] = None
    excl = d["exclusion"]
    pts = L.normalise_and_clip(points, lambda m: abs(m - tc) > excl)
    use = [p for p in pts if not p["clipped"]]
    res["n_points"] = len(pts)
    res["n_clipped"] = len(pts) - len(use)
    lc_path = os.path.join(out_dir, "lc_%s.csv" % ev["event_id"])
    with open(lc_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mjd", "days_from_tca", "band", "group", "rel", "rel_err", "y", "sigma", "err_scale",
                    "clipped", "centroid_ok", "dx_px", "dy_px", "seeing_arcsec"])
        for p in pts:
            w.writerow(["%.6f" % p["mjd"], "%.4f" % (p["mjd"] - tc), p["band"], p["group"],
                        "%.7g" % p["rel"], "%.4g" % p["rel_err"], "%.6f" % p["y"], "%.6f" % p["s"],
                        "%.3f" % p["err_scale"], int(p["clipped"]), int(bool(p["centroid_ok"])),
                        "%.2f" % p["dx"], "%.2f" % p["dy"], "%.2f" % p["seeing"]])
    if not use:
        res["verdict"] = "NO_COVERAGE"
        return res
    gidx = {g: k for k, g in enumerate(sorted(set(p["group"] for p in use)))}
    band_of = {gidx[p["group"]]: p["band"] for p in use}
    ser = L.Series([p["mjd"] for p in use], [p["y"] for p in use], [p["s"] for p in use],
                   [gidx[p["group"]] for p in use], band_of)
    t0 = time.time()
    g = L.grade(ser, tc, d["window_half"], excl, d["te_days"], d["u0"], d["u0_pess"], d["fs"],
                d["err_tca_days"])
    res["grading_seconds"] = round(time.time() - t0, 1)
    res["grading"] = g
    res["verdict"] = g["verdict"]
    mus = ev.get("mustill2018")
    if mus and g.get("observed", {}).get("t0") is not None:
        res["t0_in_mustill_window"] = mus["mjd_lo"] <= g["observed"]["t0"] <= mus["mjd_hi"] + 1.0
    if g["verdict"] in ("DETECTED", "CANDIDATE"):
        res["mass_exploratory"] = L.mass_profile(ser, g["observed"]["t0"], d["mu_rel"], d["pi_rel"],
                                                 ev["dist"], ev["err_dist"])
    # model curves
    obs = g.get("observed") or {}
    with open(os.path.join(out_dir, "model_%s.csv" % ev["event_id"]), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mjd", "days_from_tca", "band", "predicted_median", "predicted_pessimistic", "best_fit"])
        for b in sorted(set(band_of.values())):
            for k in range(-240, 241):
                t = tc + 0.25 * k
                pm = 1.0 + d["fs"].get(b, 0.0) * L.model_m(t, tc, d["u0"], d["te_days"])
                pp = 1.0 + d["fs"].get(b, 0.0) * L.model_m(t, tc, d["u0_pess"], d["te_days"])
                bf = ""
                if obs.get("t0") is not None:
                    bf = "%.6f" % (1.0 + obs["amp"].get(b, 0.0) * L.model_m(t, obs["t0"], obs["u0"], obs["tE"]))
                w.writerow(["%.4f" % t, "%.2f" % (t - tc), b, "%.6f" % pm, "%.6f" % pp, bf])
    return res


def summary_md(results):
    lines = ["# foretold-lens results v%s" % VERSION, "",
             "| event | lens | predicted peak (MJD) | verdict | D | threshold | slides | p | recovery (median/pess.) |",
             "|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        g = r.get("grading", {})
        o = g.get("observed", {})
        lines.append("| %s | %s | %.2f | %s | %s | %s | %s | %s | %s/%s |" % (
            r["event_id"], r["name"], r["derived"]["tca_mjd"], r.get("verdict"),
            "%.1f" % o["D"] if o else "-", "%.1f" % g["threshold"] if "threshold" in g else "-",
            g.get("n_slides", "-"), "%.3f" % g["p_slide"] if "p_slide" in g else "-",
            "%.2f" % g["recovery_median"] if g.get("recovery_median") is not None else "-",
            "%.2f" % g["recovery_pess"] if g.get("recovery_pess") is not None else "-"))
    return "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--events", default="", help="comma-separated event_ids (default: all)")
    ap.add_argument("--max-images", type=int, default=MAX_IMAGES)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(a.data_dir, exist_ok=True)
    want = set(x for x in a.events.split(",") if x)
    evs = sorted(E.EVENTS, key=lambda e: e["event_id"] != E.HEADLINE)
    results = []
    for ev in evs:
        if want and ev["event_id"] not in want:
            continue
        log("=== event", ev["event_id"], ev["name"])
        try:
            r = process_event(ev, a.data_dir, a.out_dir, a.max_images)
        except Exception as e:
            traceback.print_exc()
            r = dict(event_id=ev["event_id"], name=ev["name"], derived=E.derived(ev),
                     verdict="PIPELINE_ERROR", error=str(e)[:500])
        log("verdict", ev["event_id"], r.get("verdict"))
        results.append(r)
        with open(os.path.join(a.out_dir, "results.json"), "w") as f:
            json.dump(dict(version=VERSION, generated=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                           events=results), f, indent=1, default=str)
        with open(os.path.join(a.out_dir, "summary.md"), "w") as f:
            f.write(summary_md(results))
    log("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
