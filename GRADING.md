# Grading of foretold-lens v0.1.0 (job foretold-lens-v1)

Graded 2026-11-28 by AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com), against
`PREREGISTRATION.md` (frozen 2026-11-24, before any photometry of the event positions was measured).
No human expert has checked it yet. The job ran once on the Dell, 2026-11-27 12:01-12:26 UTC (unit
tests and analysis exit 0, peak memory 0.12 GB). Its output is committed unchanged in
`results/v0.1.0/`. (It finished on 27 November; grading was delayed by a day because our run checklist
missed it.)

## 1. Pre-registered verdicts

| event | lens | predicted peak | verdict | D_obs | slide threshold T | slides | p | recovery (median / pessimistic) |
|---|---|---|---|---|---|---|---|---|
| 2042703905028908928-32410 | **G 207-14** | 2025-07-19.7 +- 2.4 d | **MARGINAL** | 15.0 | 65.9 | 34 | 0.343 | 0.41 / 0.18 |
| 2030898139472914688-4054 | DR3 2030898139472914688 | 2025-09-07 | NO_COVERAGE | - | - | - | - | - |
| 4288183244850920192-10753 | DR3 4288183244850920192 | 2024-06-19 | NO_COVERAGE | - | - | - | - | - |
| 514864534672816640-10470 | DR3 514864534672816640 | 2025-04-18 | NO_COVERAGE | - | - | - | - | - |
| 516896775755948416-58094 | DR3 516896775755948416 | 2021-02-05 | NOT_SEEN_UNCALIBRATED | 7.1 | 25.0 | 7 | 0.125 | 0.00 / 0.00 |
| 4078382201351836288-8 | DR3 4078382201351836288 | 2025-01-04 | NO_COVERAGE | - | - | - | - | - |

## 2. G 207-14 (the headline event): not seen, and the test was not sensitive enough to decide

- 377 forced-photometry points from 5 ZTF groups (g, r and i; fields 684 and 685), 2024-05 to 2025-10.
- The best point-lens fit inside the search window has D = 15.0, which is `MARGINAL` (9 <= D < 25).
  It is not a convincing lensing signal: its peak sits at t0 = tca + 10.0 d, the very edge of the
  +-10 d window (4.1 times Kluter et al.'s timing error), its u0 = 3.0 is the largest value on the
  grid, and it implies a g-band bump of about 1.7%, not the predicted 4.5%. It is also not unusual:
  11 of the 34 time slides away from the prediction give D >= 15 (p = 0.34).
- **Sensitivity.** The slides show that this field's ZTF photometry has month-scale systematics
  (slide D values up to 65.9, at slides 50-100 days either side of the prediction), so the frozen threshold is
  T = 65.9. Injected at the slide centres, the median prediction (u0 = 0.86, tE = 9.1 d, source
  share 9.7% in g and 3.2% in r, i.e. +4.5% in g and +1.5% in r at peak) was recovered above T in 41%
  of slides, the pessimistic one in 18%. Under the frozen rules this test cannot say whether the
  predicted brightening happened.
- Reported, not a gate: the fit's g amplitude exceeds its r amplitude, as expected for a source bluer
  than the lens; t0 is outside Mustill et al.'s 2025-07-20/27 window.
- Post hoc (not pre-registered; read with care): at 28 of the 34 slide centres the injected median
  prediction gave a larger D than the 15.0 found at the predicted time (23 of 34 for the pessimistic
  prediction). That hints the median prediction would usually have left a bigger mark than what is
  seen, but the injected D values include each slide's own systematics, so this is not a measurement.

## 3. The other events

- Four events have too few ZTF images near their predicted times (`NO_COVERAGE`): ZTF's public
  survey does not revisit their fields often enough, or they fall in seasonal gaps.
- DR3 516896775755948416 (February 2021) is `NOT_SEEN_UNCALIBRATED`: only 7 slides, and the
  predicted signal (source share 7.5%) was recovered in none of them, so the test had no sensitivity.

## 4. Verdict

No Gaia-predicted photometric microlensing event is detected in public ZTF data. For G 207-14, the
best case, the predicted 4.5% g-band brightening in July 2025 is neither seen nor excluded: ZTF's
systematics in this field are as large as the expected signal. Ground photometry with a better
systematics floor (ATLAS forced photometry, or a dedicated difference-imaging reduction of the same
ZTF images) is the way to settle it. Per the pre-registration an article may follow for a graded
label; we judge this result too weak for one on its own and will fold it into a later write-up only
if a better reduction settles the event.

## 5. Credit

Predictions: A. J. Mustill, M. B. Davies and L. Lindegren (2018, A&A 617, A135; arXiv:1805.11638);
J. Kluter, U. Bastian, M. Demleitner and J. Wambsganss (2022, AJ 163, 176; arXiv:2112.12152). Data:
ZTF public survey (NASA/IPAC IRSA), ESA Gaia.
