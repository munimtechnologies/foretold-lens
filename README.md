# foretold-lens

**Did the microlensing events that Gaia foretold happen on schedule?**

In 2018, Mustill, Davies & Lindegren used Gaia DR2 to predict that the nearby red dwarf G 207-14 would
pass in front of a faint background star in late July 2025 and briefly brighten it by gravitational
microlensing. It was the most probable event in their 20-year list. Kluter et al. (2022) re-predicted it
from Gaia eDR3: closest approach 2025-07-19.7 +- 2.4 days.

The public Zwicky Transient Facility (ZTF) imaged the field about once per night in two filters right
through that window. This repository holds a pre-registered test, frozen before any photometry was
measured, of that event and of every other Gaia-predicted event in ZTF's sky and public era.

- `PREREGISTRATION.md`: question, frozen sample, photometry, statistics and verdict labels
- `events.py`: the six frozen events (Kluter et al. 2022 catalogue values, copied verbatim)
- `ztf_io.py`: IRSA IBE image search and cut-outs (public ZTF data)
- `phot.py`: minimal FITS reader and forced differential aperture photometry (stdlib only)
- `lensfit.py`: point-lens search, time-slide false-alarm calibration, injection-recovery, verdicts
- `job.py`: runs everything, `python job.py --data-dir data --out-dir out`
- `test_foretold.py`: unit tests (`python -m unittest -v test_foretold`)
- `notes/scoop-check-2026-11-24.md`: prior-work search

Status: pre-registered 2026-11-24; the job is queued. Results will appear in `results/`.

Credits: the predictions are by A. J. Mustill, M. B. Davies and L. Lindegren (2018, A&A 617, A135) and
J. Kluter, U. Bastian, M. Demleitner and J. Wambsganss (2022, AJ 163, 176). The data are from ZTF (public
survey, via NASA/IPAC IRSA) and ESA Gaia.

Made with AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com). No human expert has
checked it yet. Licence: Apache-2.0.
