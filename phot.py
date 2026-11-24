# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Minimal FITS reader and forced differential aperture photometry on ZTF cut-outs (stdlib only).

Every cut-out of one group (same ZTF field, CCD, quadrant and filter) is centred by IRSA on the same
sky position (the event position at tca), so the target sits at the cut-out centre and comparison
stars sit at fixed offsets from it. No WCS is needed; positions are re-centroided in every image.

Made with AI agents working in MV2's name (MV2, Munim, Inc.).
"""
import math
import sys
from array import array

# frozen photometry parameters (PREREGISTRATION.md section 4)
R_AP_MIN, R_AP_MAX, R_AP_FWHM = 2.5, 5.0, 1.0     # aperture radius = clamp(1.0 * FWHM, 2.5, 5.0) px
R_IN, R_OUT = 12.0, 20.0                          # sky annulus, px
EDGE = 22                                         # comparison stars at least this far from the edge
TARGET_SEARCH = 3.0                               # px, target re-centroid limit
COMP_SEARCH = 2.0                                 # px, comparison re-centroid limit
N_COMP_MAX = 20
COMP_NSIG, ISO_NSIG, ISO_SEP = 30.0, 10.0, 8.0
SAT_FRAC = 0.6
MIN_TARGET_SEP = 15.0
COMP_KEEP_FRAC = 0.9
MIN_COMPS = 3
MIN_GROUP_IMAGES = 10


# ----------------------------------------------------------------------------- FITS
def _parse_value(s):
    s = s.strip()
    if s.startswith("'"):
        end = s.find("'", 1)
        while end != -1 and end + 1 < len(s) and s[end + 1] == "'":
            end = s.find("'", end + 2)
        return s[1:end if end != -1 else None].replace("''", "'").rstrip()
    v = s.split("/", 1)[0].strip()
    if v == "T":
        return True
    if v == "F":
        return False
    try:
        if any(c in v for c in ".EeDd"):
            return float(v.replace("D", "E").replace("d", "e"))
        return int(v)
    except ValueError:
        return v


def read_fits(blob):
    """Return (header dict, nx, ny, flat list of floats) for the primary HDU of a FITS file."""
    if not blob[:6] == b"SIMPLE":
        raise ValueError("not a FITS file")
    hdr, pos = {}, 0
    while True:
        card = blob[pos:pos + 80].decode("ascii", "replace")
        if len(card) < 80:
            raise ValueError("truncated FITS header")
        pos += 80
        key = card[:8].strip()
        if key == "END":
            break
        if card[8:10] == "= ":
            hdr[key] = _parse_value(card[10:])
    data_start = int(math.ceil(pos / 2880.0)) * 2880
    bitpix = int(hdr["BITPIX"])
    nx, ny = int(hdr.get("NAXIS1", 0)), int(hdr.get("NAXIS2", 0))
    n = nx * ny
    code = {-32: "f", -64: "d", 16: "h", 32: "i", 8: "B"}[bitpix]
    a = array(code)
    size = a.itemsize
    if bitpix == 32 and size != 4:
        a = array("l") if array("l").itemsize == 4 else array("i")
        size = a.itemsize
    raw = blob[data_start:data_start + n * abs(bitpix) // 8]
    if len(raw) < n * abs(bitpix) // 8:
        raise ValueError("truncated FITS data")
    a.frombytes(raw)
    if sys.byteorder == "little" and bitpix != 8:
        a.byteswap()
    bscale, bzero = float(hdr.get("BSCALE", 1.0)), float(hdr.get("BZERO", 0.0))
    if bscale != 1.0 or bzero != 0.0:
        data = array("d", (v * bscale + bzero for v in a))
    elif code in ("f", "d"):
        data = a  # compact storage (a few hundred cut-outs stay in memory at once)
    else:
        data = array("d", (float(v) for v in a))
    return hdr, nx, ny, data


def write_fits(hdr_items, nx, ny, data):
    """Tiny writer (BITPIX -32) used by the tests."""
    cards = ["SIMPLE  = %20s" % "T", "BITPIX  = %20d" % -32, "NAXIS   = %20d" % 2,
             "NAXIS1  = %20d" % nx, "NAXIS2  = %20d" % ny]
    for k, v in hdr_items:
        if isinstance(v, str):
            cards.append("%-8s= '%s'" % (k, v))
        else:
            cards.append("%-8s= %20s" % (k, repr(v)))
    cards.append("END")
    h = "".join(c.ljust(80) for c in cards)
    h += " " * ((2880 - len(h) % 2880) % 2880)
    a = array("f", data)
    if sys.byteorder == "little":
        a.byteswap()
    body = a.tobytes()
    body += b"\0" * ((2880 - len(body) % 2880) % 2880)
    return h.encode("ascii") + body


# ----------------------------------------------------------------------------- statistics
def median(v):
    s = sorted(v)
    n = len(s)
    if n == 0:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def mad_sigma(v, med=None):
    if not v:
        return float("nan")
    m = median(v) if med is None else med
    return 1.4826 * median([abs(x - m) for x in v])


def clipped_stats(vals, nsig=3.0, iters=3):
    v = [x for x in vals if x == x]
    for _ in range(iters):
        if len(v) < 5:
            break
        m = median(v)
        s = mad_sigma(v, m)
        if not s > 0:
            break
        nv = [x for x in v if abs(x - m) <= nsig * s]
        if len(nv) == len(v):
            break
        v = nv
    m = median(v)
    return m, mad_sigma(v, m), len(v)


# ----------------------------------------------------------------------------- image helpers
class Image:
    def __init__(self, hdr, nx, ny, data):
        self.hdr, self.nx, self.ny, self.data = hdr, nx, ny, data
        finite = [x for x in data if x == x and abs(x) != float("inf")]
        self.bkg, self.sig, _ = clipped_stats(finite[::3] if len(finite) > 30000 else finite)

    def px(self, x, y):
        return self.data[y * self.nx + x]

    def center(self):
        return (self.nx - 1) / 2.0, (self.ny - 1) / 2.0


def recentroid(img, x0, y0, radius, sky, iters=3, max_shift=None):
    max_shift = radius if max_shift is None else max_shift
    x, y = x0, y0
    for _ in range(iters):
        sw = sx = sy = 0.0
        for j in range(int(math.floor(y - radius)), int(math.ceil(y + radius)) + 1):
            if j < 0 or j >= img.ny:
                continue
            for i in range(int(math.floor(x - radius)), int(math.ceil(x + radius)) + 1):
                if i < 0 or i >= img.nx or (i - x) ** 2 + (j - y) ** 2 > radius * radius:
                    continue
                v = img.px(i, j)
                if v != v:
                    continue
                v -= sky
                if v > 0:
                    sw += v
                    sx += v * i
                    sy += v * j
        if sw <= 0:
            return x0, y0, False
        nx_, ny_ = sx / sw, sy / sw
        if math.hypot(nx_ - x0, ny_ - y0) > max_shift:
            return x0, y0, False
        x, y = nx_, ny_
    return x, y, True


def annulus_sky(img, cx, cy, r_in=R_IN, r_out=R_OUT):
    vals = []
    for j in range(int(math.floor(cy - r_out)), int(math.ceil(cy + r_out)) + 1):
        if j < 0 or j >= img.ny:
            continue
        for i in range(int(math.floor(cx - r_out)), int(math.ceil(cx + r_out)) + 1):
            if i < 0 or i >= img.nx:
                continue
            d2 = (i - cx) ** 2 + (j - cy) ** 2
            if r_in * r_in <= d2 < r_out * r_out:
                v = img.px(i, j)
                if v == v:
                    vals.append(v)
    if len(vals) < 20:
        return None
    return clipped_stats(vals)


def aperture(img, cx, cy, r, sky, sub=5):
    """Sky-subtracted flux, summed weight (pixels) and positive flux. None if a NaN is inside."""
    flux = npix = pos = 0.0
    for j in range(int(math.floor(cy - r - 1)), int(math.ceil(cy + r + 1)) + 1):
        for i in range(int(math.floor(cx - r - 1)), int(math.ceil(cx + r + 1)) + 1):
            d = math.hypot(i - cx, j - cy)
            if d >= r + 0.71:
                continue
            if i < 0 or j < 0 or i >= img.nx or j >= img.ny:
                return None
            if d <= r - 0.71:
                w = 1.0
            else:
                inside = 0
                for a in range(sub):
                    for b in range(sub):
                        xx = i - 0.5 + (a + 0.5) / sub
                        yy = j - 0.5 + (b + 0.5) / sub
                        if (xx - cx) ** 2 + (yy - cy) ** 2 <= r * r:
                            inside += 1
                w = inside / float(sub * sub)
                if w == 0:
                    continue
            v = img.px(i, j)
            if v != v:
                return None
            flux += w * (v - sky)
            npix += w
            pos += w * max(v - sky, 0.0)
    return flux, npix, pos


def measure_star(img, x, y, r, gain, search, max_shift):
    """Re-centroid then aperture photometry. Returns dict or None."""
    s0 = annulus_sky(img, x, y)
    if s0 is None:
        return None
    cx, cy, ok = recentroid(img, x, y, search, s0[0], max_shift=max_shift)
    sk = annulus_sky(img, cx, cy) if ok else s0
    if sk is None:
        return None
    sky, ssig, nsky = sk
    ap = aperture(img, cx, cy, r, sky)
    if ap is None:
        return None
    flux, npix, pos = ap
    var = pos / gain + npix * ssig * ssig + npix * npix * ssig * ssig / max(nsky, 1)
    return dict(x=cx, y=cy, flux=flux, var=var, centroid_ok=ok)


def find_stars(img, nsig=COMP_NSIG, iso_nsig=ISO_NSIG, iso_sep=ISO_SEP, edge=EDGE,
               sat=None, sat_frac=SAT_FRAC):
    """Isolated, unsaturated local maxima above bkg + nsig*sigma, brightest first: [(peak, x, y)]."""
    lo = img.bkg + iso_nsig * img.sig
    hi = img.bkg + nsig * img.sig
    peaks = []
    for j in range(2, img.ny - 2):
        row = j * img.nx
        for i in range(2, img.nx - 2):
            v = img.data[row + i]
            if not (v >= lo):
                continue
            is_max = True
            for jj in range(j - 2, j + 3):
                for ii in range(i - 2, i + 3):
                    if (ii != i or jj != j) and img.data[jj * img.nx + ii] > v:
                        is_max = False
                        break
                if not is_max:
                    break
            if is_max:
                peaks.append((v, i, j))
    out = []
    for v, i, j in peaks:
        if v < hi or i < edge or j < edge or i >= img.nx - edge or j >= img.ny - edge:
            continue
        if sat is not None and v >= sat_frac * sat:
            continue
        if any((p[1] - i) ** 2 + (p[2] - j) ** 2 < iso_sep * iso_sep for p in peaks if (p[1], p[2]) != (i, j)):
            continue
        out.append((v, i, j))
    out.sort(reverse=True)
    return out


def ap_radius(fwhm_px):
    if fwhm_px is None or not fwhm_px > 0:
        fwhm_px = 2.0
    return min(R_AP_MAX, max(R_AP_MIN, R_AP_FWHM * fwhm_px))


def measure_group(items, log=print):
    """items: list of dicts with keys mjd, img (Image), meta. All of one field/ccd/quadrant/filter.

    Returns list of dicts (mjd, rel, rel_err, meta) with rel = target / comparison-ensemble flux
    (not yet normalised), plus a dict of diagnostics.
    """
    diag = dict(n_in=len(items))
    if len(items) < MIN_GROUP_IMAGES:
        diag["dropped"] = "fewer than %d images" % MIN_GROUP_IMAGES
        return [], diag
    ref = min(items, key=lambda it: float(it["img"].hdr.get("SEEING", 99.0) or 99.0))
    rimg = ref["img"]
    rcx, rcy = rimg.center()
    sat = rimg.hdr.get("SATURATE")
    stars = find_stars(rimg, sat=float(sat) if sat is not None else None)
    comps = []
    for v, i, j in stars:
        if math.hypot(i - rcx, j - rcy) < MIN_TARGET_SEP:
            continue
        comps.append((i - rcx, j - rcy))
        if len(comps) >= N_COMP_MAX:
            break
    diag["n_comp_candidates"] = len(comps)
    if len(comps) < MIN_COMPS:
        diag["dropped"] = "fewer than %d comparison stars" % MIN_COMPS
        return [], diag
    rows = []
    for it in items:
        img = it["img"]
        gain = float(img.hdr.get("GAIN", 6.1) or 6.1)
        r = ap_radius(img.hdr.get("SEEING"))
        cx, cy = img.center()
        t = measure_star(img, cx, cy, r, gain, TARGET_SEARCH, TARGET_SEARCH)
        cf = []
        for dx, dy in comps:
            m = measure_star(img, cx + dx, cy + dy, r, gain, COMP_SEARCH, COMP_SEARCH)
            cf.append(m if (m is not None and m["flux"] > 0) else None)
        rows.append((it, t, cf))
    keep = [k for k in range(len(comps))
            if sum(1 for _, _, cf in rows if cf[k] is not None) >= COMP_KEEP_FRAC * len(rows)]
    rows = [(it, t, cf) for it, t, cf in rows
            if t is not None and t["flux"] > 0 and all(cf[k] is not None for k in keep)]
    if len(keep) < MIN_COMPS or len(rows) < MIN_GROUP_IMAGES:
        diag["dropped"] = "too few usable comparison stars or images"
        return [], diag
    # reject variable comparison stars (one pass)
    scat = {}
    for k in keep:
        q = []
        for _, _, cf in rows:
            rest = sum(cf[j]["flux"] for j in keep if j != k)
            q.append(cf[k]["flux"] / rest if rest > 0 else float("nan"))
        q = [x for x in q if x == x]
        m = median(q)
        scat[k] = mad_sigma([x / m for x in q]) if m > 0 else float("inf")
    ms = median(list(scat.values()))
    keep2 = [k for k in keep if scat[k] <= 3.0 * ms] if ms > 0 else keep
    if len(keep2) < MIN_COMPS:
        keep2 = keep
    diag["n_comps_used"] = len(keep2)
    out = []
    for it, t, cf in rows:
        e = sum(cf[k]["flux"] for k in keep2)
        ve = sum(cf[k]["var"] for k in keep2)
        if e <= 0:
            continue
        rel = t["flux"] / e
        rerr = rel * math.sqrt(t["var"] / t["flux"] ** 2 + ve / e ** 2)
        out.append(dict(mjd=it["mjd"], rel=rel, rel_err=rerr, meta=it["meta"],
                        centroid_ok=t["centroid_ok"], dx=t["x"] - img_center_x(it), dy=t["y"] - img_center_y(it)))
    diag["n_out"] = len(out)
    return out, diag


def img_center_x(it):
    return it["img"].center()[0]


def img_center_y(it):
    return it["img"].center()[1]
