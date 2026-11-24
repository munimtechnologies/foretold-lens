# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Frozen event list for the foretold-lens test (pre-registered 2026-11-24).

Selection rule (fixed before any photometry was measured): every event in the Kluter et al. (2022,
AJ 163, 176; arXiv:2112.12152) Gaia eDR3 prediction catalogue (GAVO table plc3.data) with
  magnification > 0.003 mag, ob_dec > -28 deg, 2018.21 <= tca <= 2025.80
(ZTF start 2018-03-17 to ZTF DR24's public cut-off 2025-10-22). Queried 2026-11-24 via
https://dc.g-vo.org/tap/sync ; values copied verbatim. Positions/proper motions are Gaia eDR3
(epoch J2016.0). Mustill, Davies & Lindegren (2018, A&A 617, A135; arXiv:1805.11638) predicted
two of these events independently from Gaia DR2 (windows given below).

Made with AI agents working in MV2's name (MV2, Munim, Inc.).
"""
import math

KAPPA = 8.144  # mas / Msun
GAIA_EPOCH = 2016.0
DR24_END_MJD = 60970.0  # ZTF DR24 public survey cut-off, 2025-10-22


def jyear_to_mjd(y):
    """Julian year -> MJD (J2000.0 = JD 2451545.0 = MJD 51544.5)."""
    return 51544.5 + (y - 2000.0) * 365.25


EVENTS = [
    dict(event_id="2042703905028908928-32410", name="G 207-14",
         lens_ra=286.431100010524, lens_dec=31.1767774645684, lens_pmra=-64.01847, lens_pmdec=-252.18665,
         ob_ra=286.430901991362, ob_dec=31.1761294208638, ob_pmra=-0.21782744, ob_pmdec=-5.919696,
         lens_parallax=27.205883, ob_parallax=0.51627,
         lens_g=14.6741, lens_bp=16.4641, lens_rp=13.404, ob_g=18.5897, ob_bp=18.8846, ob_rp=17.1151,
         mass=0.18661565, theta_e=6.369548986414099, err_theta_e=0.31889787,
         tca=2025.5475410761549, err_tca=0.006697983, dist=5.4784546, err_dist=1.6347889,
         u=0.8602497767307322, u_err_m=-0.25892273, u_err_p=0.2695651,
         magnification=0.01320609399127736, mag_err_m=-0.0055731363, mag_err_p=0.012331393,
         mustill2018=dict(window="2025-07-20/2025-07-27", mjd_lo=60876.0, mjd_hi=60883.0,
                          peak_net_magnification=0.0213, p_within_einstein=0.7677,
                          dr2_lens="2042703905028908928", dr2_source="2042703900727995008")),
    dict(event_id="2030898139472914688-4054", name="DR3 2030898139472914688",
         lens_ra=300.453087988931, lens_dec=31.3678957964694, lens_pmra=-36.323368, lens_pmdec=-305.56668,
         ob_ra=300.452973004965, ob_dec=31.3670913654108, ob_pmra=-1.4060253, ob_pmdec=-5.2472854,
         lens_parallax=23.409584, ob_parallax=0.33124173,
         lens_g=15.5855, lens_bp=17.4421, lens_rp=14.2807, ob_g=19.8842, ob_bp=20.8378, ob_rp=18.3211,
         mass=0.15646687, theta_e=5.42346723682085, err_theta_e=0.27485916,
         tca=2025.6821564607703, err_tca=0.013350332, dist=1.4879596, err_dist=3.5142505,
         u=0.2731994496630562, u_err_m=-0.65441084, u_err_p=0.6739963,
         magnification=0.025180824389252336, mag_err_m=-0.018281274, mag_err_p=0.097899824,
         mustill2018=dict(window="2025-09-01/2025-09-16", mjd_lo=60919.0, mjd_hi=60934.0,
                          peak_net_magnification=0.00619, p_within_einstein=0.4929,
                          dr2_lens="2030898139472914688", dr2_source="2030898105088410112")),
    dict(event_id="4288183244850920192-10753", name="DR3 4288183244850920192",
         lens_ra=295.128493067148, lens_dec=2.51198560283781, lens_pmra=-127.09552, lens_pmdec=-103.8133,
         ob_ra=295.128229338603, ob_dec=2.51176125556963, ob_pmra=-15.05553, ob_pmdec=-7.275136,
         lens_parallax=9.61684, ob_parallax=1.3468164,
         lens_g=16.822, lens_bp=None, lens_rp=None, ob_g=14.7685, ob_bp=15.1648, ob_rp=14.0593,
         mass=0.19291635, theta_e=3.6049680849853294, err_theta_e=0.18395709,
         tca=2024.463681849529, err_tca=0.019087374, dist=6.925381, err_dist=2.658706,
         u=1.918488878121712, u_err_m=-0.7462721, u_err_p=0.77420557,
         magnification=0.06269165832462897, mag_err_m=-0.039440528, mag_err_p=0.15247862,
         mustill2018=None),
    dict(event_id="514864534672816640-10470", name="DR3 514864534672816640",
         lens_ra=32.4160886604094, lens_dec=63.4212387316894, lens_pmra=113.56117, lens_pmdec=-55.752907,
         ob_ra=32.4166314345204, ob_dec=63.4211097060611, ob_pmra=19.022787, ob_pmdec=-7.8787284,
         lens_parallax=11.568782, ob_parallax=2.3908858,
         lens_g=16.6855, lens_bp=None, lens_rp=None, ob_g=12.6817, ob_bp=13.1067, ob_rp=12.0388,
         mass=0.17812964, theta_e=3.6492464289358018, err_theta_e=0.18522084,
         tca=2025.2975410761549, err_tca=0.018743085, dist=12.728736, err_dist=3.173371,
         u=3.5033915372278766, u_err_m=-0.88484865, u_err_p=0.9190582,
         magnification=0.01050248028679354, mag_err_m=-0.005818336, mag_err_p=0.01813863,
         mustill2018=None),
    dict(event_id="516896775755948416-58094", name="DR3 516896775755948416",
         lens_ra=37.0155937006301, lens_dec=65.8992118523984, lens_pmra=147.21892, lens_pmdec=-138.63506,
         ob_ra=37.0160966844952, ob_dec=65.8990132472717, ob_pmra=-1.2620691, ob_pmdec=0.27209398,
         lens_parallax=15.100111, ob_parallax=0.28659385,
         lens_g=17.053, lens_bp=18.9783, lens_rp=15.6736, ob_g=19.7847, ob_bp=None, ob_rp=None,
         mass=0.13396877, theta_e=4.020632687409326, err_theta_e=0.21551938,
         tca=2021.0982972341442, err_tca=0.02325712, dist=4.5665326, err_dist=4.984569,
         u=1.1216649290408582, u_err_m=-1.2265277, u_err_p=1.2948362,
         magnification=0.017964560507177152, mag_err_m=-0.015087696, mag_err_p=0.12726928,
         mustill2018=None),
    dict(event_id="4078382201351836288-8", name="DR3 4078382201351836288",
         lens_ra=280.059823894269, lens_dec=-23.3491144234504, lens_pmra=-216.46521, lens_pmdec=-179.44998,
         ob_ra=280.059242574598, ob_dec=-23.349549410231, ob_pmra=-3.9906933, ob_pmdec=-5.0438766,
         lens_parallax=6.547443, ob_parallax=0.007311601,
         lens_g=19.6559, lens_bp=None, lens_rp=None, ob_g=19.5816, ob_bp=20.2773, ob_rp=18.9664,
         mass=0.10565399, theta_e=2.3724642965326908, err_theta_e=0.15629895,
         tca=2025.0101512315396, err_tca=0.018550381, dist=9.217305, err_dist=4.828212,
         u=3.890480680472494, u_err_m=-2.025027, u_err_p=2.1177623,
         magnification=0.0038544041416052055, mag_err_m=-0.0030481992, mag_err_p=0.039502468,
         mustill2018=None),
]

HEADLINE = "2042703905028908928-32410"


def source_fraction(m_lens, m_src):
    """Fraction of the blended flux that comes from the (lensed) source."""
    if m_lens is None or m_src is None:
        return None
    return 1.0 / (1.0 + 10 ** (-0.4 * (m_lens - m_src)))


def derived(ev):
    """Quantities used by the frozen grading rules."""
    dpmra = ev["lens_pmra"] - ev["ob_pmra"]
    dpmdec = ev["lens_pmdec"] - ev["ob_pmdec"]
    mu_rel = math.hypot(dpmra, dpmdec)  # mas/yr
    te_days = ev["theta_e"] / mu_rel * 365.25
    fs_G = source_fraction(ev["lens_g"], ev["ob_g"])
    fs_bp = source_fraction(ev["lens_bp"], ev["ob_bp"])
    fs_rp = source_fraction(ev["lens_rp"], ev["ob_rp"])
    fs = {"g": fs_bp if fs_bp is not None else fs_G, "r": fs_rp if fs_rp is not None else fs_G}
    fs["i"] = fs["r"]
    tca_mjd = jyear_to_mjd(ev["tca"])
    err_days = ev["err_tca"] * 365.25
    w = max(10.0, 3.0 * err_days)
    u0 = ev["u"]
    u0_pess = ev["u"] + ev["u_err_p"]
    # cut-out centre: the brighter object (in G) at tca; lens and source are < 0.1" apart then
    if ev["lens_g"] <= ev["ob_g"]:
        ra0, dec0, pmra, pmdec = ev["lens_ra"], ev["lens_dec"], ev["lens_pmra"], ev["lens_pmdec"]
    else:
        ra0, dec0, pmra, pmdec = ev["ob_ra"], ev["ob_dec"], ev["ob_pmra"], ev["ob_pmdec"]
    dy = ev["tca"] - GAIA_EPOCH
    dec = dec0 + pmdec * dy / 3.6e6
    ra = ra0 + pmra * dy / 3.6e6 / math.cos(math.radians(dec0))
    pi_rel = ev["lens_parallax"] - ev["ob_parallax"]
    return dict(mu_rel=mu_rel, te_days=te_days, fs=fs, tca_mjd=tca_mjd, err_tca_days=err_days,
                window_half=w, u0=u0, u0_pess=u0_pess, center_ra=ra, center_dec=dec, pi_rel=pi_rel,
                exclusion=w + max(30.0, 4.0 * te_days))
