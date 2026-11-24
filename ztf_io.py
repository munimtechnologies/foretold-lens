# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""ZTF image metadata search and science-image cut-outs from IRSA IBE (public data; stdlib only).

Made with AI agents working in MV2's name (MV2, Munim, Inc.).
"""
import csv
import io
import os
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "foretold-lens/0.1 (open research code; MV2, Munim, Inc.)"}
IBE_SEARCH = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci"
IBE_DATA = "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci"
COLUMNS = "obsjd,filtercode,field,ccdid,qid,imgtypecode,filefracday,seeing,maglimit,infobits"
CUTOUT_ARCSEC = 180

# frozen image-quality cuts (PREREGISTRATION.md section 4)
MAX_SEEING = 4.0
MIN_MAGLIMIT = 19.5


def fetch(url, timeout=90, tries=4, log=print):
    last = None
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            last = e
            log("fetch failed (%d/%d) %s: %s" % (a + 1, tries, url[:180], e))
            time.sleep(5 * (a + 1))
    raise RuntimeError("fetch failed: %s" % last)


def search_url(ra, dec, jd_lo, jd_hi):
    q = urllib.parse.urlencode(dict(POS="%.6f,%.6f" % (ra, dec), ct="csv",
                                    WHERE="obsjd BETWEEN %.5f AND %.5f" % (jd_lo, jd_hi),
                                    COLUMNS=COLUMNS))
    return IBE_SEARCH + "?" + q


def search(ra, dec, jd_lo, jd_hi, log=print):
    text = fetch(search_url(ra, dec, jd_lo, jd_hi), log=log).decode("utf-8", "replace")
    if text.lstrip().startswith("<"):
        raise RuntimeError("IBE search returned XML/HTML: %s" % text[:300])
    return list(csv.DictReader(io.StringIO(text)))


def passes_quality(row):
    try:
        return (row.get("imgtypecode", "o") == "o" and int(float(row["infobits"])) == 0
                and float(row["seeing"]) <= MAX_SEEING and float(row["maglimit"]) >= MIN_MAGLIMIT)
    except (KeyError, ValueError, TypeError):
        return False


def cutout_url(row, ra, dec, size_arcsec=CUTOUT_ARCSEC):
    ffd = str(row["filefracday"]).strip()
    fname = "ztf_%s_%06d_%s_c%02d_%s_q%d_sciimg.fits" % (
        ffd, int(row["field"]), row["filtercode"].strip(), int(row["ccdid"]),
        row.get("imgtypecode", "o").strip() or "o", int(row["qid"]))
    path = "/".join([IBE_DATA, ffd[0:4], ffd[4:8], ffd[8:], fname])
    return path + "?" + urllib.parse.urlencode(dict(center="%.6f,%.6f" % (ra, dec),
                                                    size="%darcsec" % size_arcsec, gzip="false"))


def group_key(row):
    return "%s_f%s_c%s_q%s" % (row["filtercode"].strip(), int(row["field"]), int(row["ccdid"]), int(row["qid"]))


def get_cutout(row, ra, dec, cache_dir, log=print):
    os.makedirs(cache_dir, exist_ok=True)
    name = "%s_%s.fits" % (str(row["filefracday"]).strip(), group_key(row))
    p = os.path.join(cache_dir, name)
    if os.path.exists(p) and os.path.getsize(p) > 2880:
        with open(p, "rb") as f:
            return f.read()
    blob = fetch(cutout_url(row, ra, dec), log=log, tries=3)
    if not blob[:6] == b"SIMPLE":
        raise RuntimeError("cut-out is not FITS: %r" % blob[:120])
    with open(p, "wb") as f:
        f.write(blob)
    return blob
