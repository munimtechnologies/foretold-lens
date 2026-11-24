# Pre-registration: did the microlensing events that Gaia foretold happen on schedule?

Frozen 2026-11-24, before any photometry of the event positions was measured. Licence: Apache-2.0 (code),
CC BY 4.0 (text). Made with AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com); no
human expert has checked this design yet.

## 1. The near-miss

Gaia's astrometry makes it possible to predict, years ahead, when a nearby star will pass almost exactly in
front of a distant one and briefly brighten it by gravitational microlensing.
- **Mustill, Davies & Lindegren (2018, A&A 617, A135; [arXiv:1805.11638](https://arxiv.org/abs/1805.11638))**
  searched Gaia DR2 for photometric events up to 2035. Their top event by probability (77% of passing
  within one Einstein radius) was the M dwarf **G 207-14** (Gaia DR2 2042703905028908928, G = 14.7)
  lensing the G = 18.6 star Gaia DR2 2042703900727995008, with a peak in **2025 July 20-27** and a
  peak net brightening of 2.13 (+7.59/-1.43) x 10^-2.
  Their second-best 2025 event was Gaia DR2 2030898139472914688 (peak 2025 September 1-16).
- **Kluter, Bastian, Demleitner & Wambsganss (2022, AJ 163, 176; [arXiv:2112.12152](https://arxiv.org/abs/2112.12152))**
  re-predicted events with Gaia eDR3 and published the catalogue (GAVO `plc3.data`). For G 207-14 it
  gives closest approach tca = J2025.5475 +- 0.0067 (2025-07-19.7 +- 2.4 d), impact 5.48 +- 1.63 mas,
  Einstein radius 6.37 mas, u0 = 0.86, peak net magnification 0.0132 (+0.0123/-0.0056) mag.
- McGill et al. (2020, MNRAS 498, L6; [arXiv:2006.13958](https://arxiv.org/abs/2006.13958)) showed that
  many Gaia DR2 predictions during the mission were spurious, and Kluter et al. (2022) wrote that Gaia
  "is not sufficient to securely predict photometric microlensing events".
  McGill et al. (2023) followed up an astrometric prediction (LAWD 37). As far as our search found
  (notes/scoop-check-2026-11-24.md), **nobody has reported whether any photometric event predicted from
  Gaia happened**, including G 207-14's.
- ZTF DR24 (public 2026-01-22; public-survey data up to 2025-10-22) contains about one g and one r image per
  night of the G 207-14 field from 2025-06-14 to 2025-08-18, across the predicted peak. ZTF's catalogue
  light curves lose G 207-14 after 2025-03 (it has moved about 2" since the reference images), so we
  measure the images ourselves.

## 2. Question and wow

Did G 207-14 (and every other event in the frozen sample) brighten its background star when and as Gaia
predicted? **Wow, if positive:** the first photometric microlensing event foretold years in advance from
Gaia astrometry and then seen on schedule, with a direct mass for a nearby red dwarf. That shows such
predictions work, ahead of the Gaia DR4/DR5 prediction catalogues. A clean non-detection where the
predicted signal would have been seen is also a result: it would show the prediction failed.

## 3. Frozen sample (events.py)

Every event in Kluter et al.'s eDR3 catalogue with magnification > 0.003 mag, source Dec > -28 deg and
2018.21 <= tca <= 2025.80 (ZTF start to the DR24 public cut-off): 6 events. Values were copied verbatim
from GAVO on 2026-11-24. The headline event is G 207-14. The others are graded by the same rules; those
with no ZTF data near tca will be labelled NO_COVERAGE.

## 4. Data and photometry (ztf_io.py, phot.py)

- IRSA IBE search of ZTF science images containing the event position (the brighter star's position at
  tca), with obsjd from tca - 420 d to min(tca + 150 d, 2025-10-22).
- Quality cuts: imgtypecode o, infobits = 0, seeing <= 4.0", maglimit >= 19.5. At most 1500 images,
  closest to tca first.
- 180" cut-outs centred on that position. We keep only cut-outs with the most common shape (edge cut-outs
  are dropped).
- Groups are (filter, field, CCD, quadrant). For each group the reference image is the best-seeing one.
  Comparison stars are up to 20 isolated (no 10-sigma peak within 8 px), unsaturated (peak < 0.6 SATURATE),
  >= 30-sigma stars, at least 15 px from the target and 22 px from the edge.
- Forced aperture photometry, re-centroided (target within 3 px, comparison stars within 2 px). The
  aperture radius is clamp(1.0 x FWHM, 2.5, 5) px, the sky is a 3-sigma-clipped median in a 12-20 px
  annulus, and the error comes from CCD noise.
- Comparison stars must be valid in >= 90% of images; images missing a kept comparison star are dropped.
  One pass rejects comparison stars whose robust scatter is > 3x the median. A group needs >= 3
  comparison stars and >= 10 images.
- Target flux / ensemble flux is normalised to its median over baseline epochs (|t - tca| > E, where
  E = w + max(30 d, 4 tE_pred)). Errors are scaled by max(1, 1.4826 median(|y-1|/sigma)) over the
  baseline. Points > 5 sigma from the median of >= 3 same-band neighbours within +-2 d are flagged and
  dropped (single-night flares and artefacts).

## 5. Test (lensfit.py)

- Model: y = c_group + a_band (A(u) - 1), with u = sqrt(u0^2 + ((t - t0)/tE)^2) and a_band >= 0 (the
  source's share of the blend in that band).
- Detection statistic: D = chi2_flat - chi2_lens, maximised over t0 in [tca - w, tca + w] (step 0.5 d),
  where w = max(10 d, 3 err_tca), tE in 2 x 1.3^k d (k = 0..11) and u0 in
  {0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.3, 1.7, 2.2, 3.0}.
- Time-slides (false-alarm calibration): the same search is centred every max(7 d, w) days at least E from
  tca, wherever there is coverage. A window has coverage if it has >= 8 points within +-(w + tE_pred) and
  >= 2 on each side. Threshold T = max(25, 95th percentile of the slide D values). p = (1 + #slides with
  D >= D_obs) / (1 + N_slides).
- Injection-recovery: the predicted light curve is added at each slide centre and searched the same way.
  The median prediction uses u0 = Kluter's u, tE = thetaE / mu_rel, a_g = Gaia BP source share and
  a_r = Gaia RP source share (G where BP/RP is missing). These shares are conservative for ZTF g and r,
  because the source is bluer than the red-dwarf lens. The pessimistic prediction uses u0 = u + its upper
  error. Recovered if D >= T.

## 6. Verdicts (frozen, per event)

| label | rule |
|---|---|
| NO_COVERAGE | fewer than 8 points within tca +- (w + tE_pred), or fewer than 2 on a side |
| DETECTED | N_slides >= 19 and D_obs >= T |
| CANDIDATE | D_obs >= 25 but not DETECTED (too few slides, or below the slide threshold) |
| MARGINAL | 9 <= D_obs < 25 |
| NOT_SEEN_PREDICTION_EXCLUDED | D_obs < 9, N_slides >= 10, median-prediction recovery >= 90% |
| NOT_SEEN_PREDICTION_DISFAVOURED | same, recovery 50-90% |
| NOT_SEEN_INSENSITIVE | same, recovery < 50% |
| NOT_SEEN_UNCALIBRATED | D_obs < 9 and N_slides < 10 |

Also reported (not gates):
- the best-fit t0 minus tca (in days and in units of err_tca), the Delta-D <= 1 t0 interval, and whether
  t0 falls in Mustill et al.'s window;
- whether a_g >= a_r (expected for G 207-14, whose source is bluer than the lens);
- for DETECTED/CANDIDATE only, an **exploratory** lens mass. It comes from a geometry-constrained fit
  (u0 = d0/thetaE, tE = thetaE/mu_rel with Gaia's relative proper motion, a Gaussian prior on Gaia's
  impact parameter d0), with M = thetaE^2 / (8.144 pi_rel).

## 7. What gets published

The results (results.json, lc_*.csv, model_*.csv, summary.md) are committed whatever they show. An
article on munimtech.com follows only for a graded headline verdict (any label except NO_COVERAGE or
PIPELINE_ERROR), after a fresh scoop check. Deviations from this plan will be written in a dated
DEVIATION file before re-running.

## 8. What we looked at before freezing (disclosure)

- Only image metadata: IBE counts, dates, seeing and fields.
- Counts (no magnitudes) from the ZTF DR24 light-curve service, which showed that G 207-14 has no
  catalogue light curve after 2025-03.
- The ALeRCE alert list within 10" of G 207-14: two single-detection alerts on MJD 60721 and 60826, both
  before the predicted window, and none in it. ZTF alerts fire only on >= 5-sigma difference-image
  detections, and a high-proper-motion star leaves dipoles in difference images, so this says little
  either way.

No pixel value or magnitude of any event position was looked at.
