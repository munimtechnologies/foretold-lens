# foretold-lens

**Did the microlensing events that Gaia foretold happen on schedule?**

In 2018, Mustill, Davies & Lindegren used Gaia DR2 to predict that the nearby red dwarf G 207-14 would
pass in front of a faint background star in late July 2025 and briefly brighten it by gravitational
microlensing. It was the most probable event in their 20-year list. Kluter et al. (2022) re-predicted it
from Gaia eDR3: closest approach 2025-07-19.7 +- 2.4 days.

The public Zwicky Transient Facility (ZTF) imaged the field about once per night in two filters right
through that window. This repository holds a pre-registered test, frozen before any photometry was
measured, of that event and of every other Gaia-predicted event in ZTF's sky and public era.

## Result (graded 2026-11-28, see [GRADING.md](GRADING.md); output in `results/v0.1.0/`)

- **G 207-14: `MARGINAL`, undecided.** The best fit in the window (D = 15.0) sits at the window's edge
  with implausible parameters and is matched by 11 of 34 time slides (p = 0.34). ZTF's month-scale
  systematics in this field are as large as the predicted +4.5% g-band brightening: the median
  prediction would have been recovered in only 41% of injections. The event is neither seen nor
  excluded.
- Four other events have no usable ZTF coverage; one (February 2021) was tested without sensitivity.
- Settling G 207-14 needs photometry with a lower systematics floor (ATLAS forced photometry or a
  dedicated difference-imaging reduction).

## Files

- `PREREGISTRATION.md`: question, frozen sample, photometry, statistics and verdict labels
- `events.py`: the six frozen events (Kluter et al. 2022 catalogue values, copied verbatim)
- `ztf_io.py`: IRSA IBE image search and cut-outs (public ZTF data)
- `phot.py`: minimal FITS reader and forced differential aperture photometry (stdlib only)
- `lensfit.py`: point-lens search, time-slide false-alarm calibration, injection-recovery, verdicts
- `job.py`: runs everything, `python job.py --data-dir data --out-dir out`
- `test_foretold.py`: unit tests (`python -m unittest -v test_foretold`)
- `notes/scoop-check-2026-11-24.md`: prior-work search
- `GRADING.md`, `results/v0.1.0/`: the graded run

Credits: the predictions are by A. J. Mustill, M. B. Davies and L. Lindegren (2018, A&A 617, A135) and
J. Kluter, U. Bastian, M. Demleitner and J. Wambsganss (2022, AJ 163, 176). The data are from ZTF (public
survey, via NASA/IPAC IRSA) and ESA Gaia.

Made with AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com). No human expert has
checked it yet. Licence: Apache-2.0 (code); text and results CC BY 4.0.
