"""The UNB3m kernel and toolkit module (DESIGN 2.7 and 3.6).

Checked against the publication (Leandro, Santos and Langley 2006, staged as
refs/UNB3m_Leandro_2006): its look-up values (Tables 1 and 2) at the table
nodes and seasons of eqs. (1) to (3); eqs. (4) to (7) and (9) to (11)
evaluated here from their printed form; the claims of its text and figures
(Figure 4: UNB3 and UNB3m wet delays differ by almost nothing at 15 degrees
and below and by up to about 4 cm near 60 degrees, around day 100 or 300;
the hydrostatic delay is the same in both); its accuracy statement (Table
7). Independently: the hydrostatic height factor the troposphere kernel
already has for DESIGN Q7 (the same eq. (4)), TM-23 5.4.2.1 eq. (5.63) for
the UNB3 hydrostatic delay, and the saturation vapour pressure of water at
0 and 100 degrees Celsius.
"""

import math

import pytest

from ssrcore._kernels import unb3m as k
from ssrcore.toolkit import KernelError, troposphere, unb3m

RAD = math.pi / 180.0
NODES = [15.0, 30.0, 45.0, 60.0, 75.0]
WINTER = 28.0                     # UNB3m eq. (3), cos = 1: X = Avg - Amp
SUMMER = 28.0 + 365.25 / 2.0      # cos = -1: X = Avg + Amp
QUARTER = 28.0 + 365.25 / 4.0     # cos = 0: X = Avg

# UNB3m Table 1: pressure (mbar), WVP (mbar, UNB3 only), lambda; average then amplitude.
P_AVG = [1013.25, 1017.25, 1015.75, 1011.75, 1013.00]
P_AMP = [0.00, -3.75, -2.25, -1.75, -0.50]
WVP_AVG = [26.31, 21.79, 11.66, 6.78, 4.11]
WVP_AMP = [0.00, 8.85, 7.24, 5.36, 3.39]
LAMBDA_AVG = [2.77, 3.15, 2.57, 1.81, 1.55]
LAMBDA_AMP = [0.00, 0.33, 0.46, 0.74, 0.30]
# UNB3m Table 2: relative humidity (%).
RH_AVG = [75.0, 80.0, 76.0, 77.5, 82.5]
RH_AMP = [0.0, 0.0, -1.0, -2.5, 2.5]
# UNB3m eq. (5) text and eq. (4): R, k1, k'2, k3; TM-23 5.4.2.1: g.
R, K1, K2P, K3, G = 287.054, 77.60, 16.6, 377600.0, 9.80665


def g_m(lat, h):
    """UNB3m eq. (6)."""
    return 9.784 * (1.0 - 2.66e-3 * math.cos(2.0 * lat) - 2.8e-7 * h)


def interp(row, abs_lat):
    """UNB3m eqs. (1), (2)."""
    if abs_lat <= 15.0:
        return row[0]
    if abs_lat >= 75.0:
        return row[4]
    i = min(int((abs_lat - 15.0) / 15.0), 3)
    return row[i] + (row[i + 1] - row[i]) * (abs_lat - NODES[i]) / 15.0


# ---- eqs. (1) to (3), Tables 1 and 2 ------------------------------------------------------


@pytest.mark.parametrize("i", range(5), ids=[f"{n:.0f}deg" for n in NODES])
def test_table_values_at_the_nodes_and_their_seasons(i):
    lat = NODES[i] * RAD
    for doy, sign in ((WINTER, -1.0), (SUMMER, 1.0)):
        p, rh, lam = unb3m.table(lat, doy)
        assert p == pytest.approx(P_AVG[i] + sign * P_AMP[i], abs=1e-9)
        assert rh == pytest.approx(RH_AVG[i] + sign * RH_AMP[i], abs=1e-9)
        assert lam == pytest.approx(LAMBDA_AVG[i] + sign * LAMBDA_AMP[i], abs=1e-12)


def test_latitudes_between_nodes_are_linear_and_the_ends_hold():
    for lat_deg in (0.0, 7.5, 15.0, 22.5, 37.5, 52.5, 67.5, 75.0, 89.0):
        p, rh, lam = unb3m.table(lat_deg * RAD, QUARTER)
        assert p == pytest.approx(interp(P_AVG, lat_deg), abs=1e-9)
        assert rh == pytest.approx(interp(RH_AVG, lat_deg), abs=1e-9)
        assert lam == pytest.approx(interp(LAMBDA_AVG, lat_deg), abs=1e-12)
    assert unb3m.table(0.0, 100.0) == unb3m.table(15.0 * RAD, 100.0) == unb3m.table(10.0 * RAD, 300.0)


def test_southern_hemisphere_is_shifted_half_a_year():
    # TM-23 5.4.2.1 eq. (5.62): D_min 211 in the south, 28 in the north.
    for lat_deg in (20.0, 45.0, 62.0):
        for doy in (1.0, 100.5, 250.0):
            assert unb3m.table(-lat_deg * RAD, doy) == unb3m.table(lat_deg * RAD, doy - 183.0)


def test_meteorology_takes_t0_and_beta_from_the_troposphere_kernel():
    lat, doy = 51.0 * RAD, 140.0
    met = unb3m.meteorology(lat, doy)
    assert (met.temperature, met.lapse_rate) == troposphere.unb3m_temperature(lat, doy)
    assert (met.pressure, met.humidity, met.vapour_lapse) == unb3m.table(lat, doy)
    assert met.vapour_pressure == unb3m.vapour_pressure(met.humidity, met.temperature, met.pressure)


# ---- eqs. (9) to (11) ---------------------------------------------------------------------


def test_vapour_pressure_is_eqs_9_to_11():
    for rh, t, p in ((75.0, 299.65, 1013.25), (77.0, 272.15, 1018.0), (82.5, 249.15, 1013.5)):
        es = 0.01 * math.exp(1.2378847e-5 * t * t - 1.9121316e-2 * t + 33.93711047 - 6.3431645e3 / t)
        fw = 1.00062 + 3.14e-6 * p + 5.6e-7 * (t - 273.15) ** 2
        assert unb3m.vapour_pressure(rh, t, p) == pytest.approx(rh / 100.0 * es * fw, rel=1e-14)


def test_saturation_vapour_pressure_of_water():
    # Independent anchors: 6.112 hPa at 0 degC; the boiling point at 1013.25 hPa.
    def es(t):
        return unb3m.vapour_pressure(100.0, t, 0.0) / (1.00062 + 5.6e-7 * (t - 273.15) ** 2)
    assert es(273.15) == pytest.approx(6.112, abs=0.002)
    assert es(373.15) == pytest.approx(1013.25, abs=1.0)


# ---- eqs. (4) to (7) ----------------------------------------------------------------------


@pytest.mark.parametrize("lat_deg, h, doy", [(45.0, 0.0, 28.0), (12.0, 350.0, 200.0), (-33.0, 1200.0, 15.0),
                                             (60.0, 2500.0, 118.0), (80.0, -20.0, 300.0)])
def test_zenith_delays_are_eqs_4_to_7(lat_deg, h, doy):
    lat = lat_deg * RAD
    m = unb3m.meteorology(lat, doy)
    base = 1.0 - m.lapse_rate * h / m.temperature
    lam1 = m.vapour_lapse + 1.0
    gm = g_m(lat, h)
    tm = m.temperature * (1.0 - m.lapse_rate * R / (gm * lam1))
    dh = 1e-6 * K1 * R / gm * m.pressure * base ** (G / (R * m.lapse_rate))
    dw = (1e-6 * (tm * K2P + K3) * R / (gm * lam1 - m.lapse_rate * R) * m.vapour_pressure / m.temperature
          * base ** (lam1 * G / (R * m.lapse_rate) - 1.0))
    z = unb3m.zenith_delays(lat, h, doy)
    assert z.hydrostatic == pytest.approx(dh, rel=1e-13)
    assert z.wet == pytest.approx(dw, rel=1e-13)
    assert z.total == z.hydrostatic + z.wet


def test_hydrostatic_height_dependence_is_the_q7_factor():
    # The troposphere kernel's DESIGN Q7 factor is the same eq. (4) divided by its value at 0 m.
    for lat_deg, doy in ((45.0, 100.0), (-20.0, 10.0), (70.0, 200.0)):
        lat = lat_deg * RAD
        z0 = unb3m.zenith_delays(lat, 0.0, doy).hydrostatic
        for h in (30.0, 848.0, 2000.0, 4500.0):
            ratio = unb3m.zenith_delays(lat, h, doy).hydrostatic / z0
            assert ratio == pytest.approx(troposphere.hydrostatic_height_factor(lat, h, doy), rel=1e-14)


def test_hydrostatic_agrees_with_the_unb3_form_of_tm23():
    # TM-23 eq. (5.63): 1e-6 k1 R_d P / g_m with k1 = 77.604 and a constant g_m = 9.784.
    # At 45 degrees cos(2 phi) = 0, so UNB3m's g_m of eq. (6) is 9.784 at sea level and
    # only k1 differs (77.60 against 77.604).
    lat = 45.0 * RAD
    p = unb3m.meteorology(lat, 77.0).pressure
    tm23 = 1e-6 * 77.604 * R * p / 9.784
    assert unb3m.zenith_delays(lat, 0.0, 77.0).hydrostatic == pytest.approx(tm23 * 77.60 / 77.604, rel=1e-12)
    # Elsewhere they differ by the latitude term of eq. (6), up to 2.66e-3 of the delay,
    # and by k1, 0.004 / 77.6.
    for lat_deg in (0.0, 15.0, 75.0, 89.0):
        lat = lat_deg * RAD
        p = unb3m.meteorology(lat, 77.0).pressure
        tm23 = 1e-6 * 77.604 * R * p / 9.784
        bound = (2.66e-3 * abs(math.cos(2.0 * lat)) + 0.004 / 77.6) * 1.0001
        assert abs(unb3m.zenith_delays(lat, 0.0, 77.0).hydrostatic - tm23) < bound * tm23


def test_magnitudes_at_sea_level():
    # About 2.3 m hydrostatic; wet from under 2 cm (polar winter) to about 0.3 m (tropics).
    for lat_deg in range(-85, 90, 10):
        for doy in (1.0, 91.0, 182.0, 274.0):
            z = unb3m.zenith_delays(lat_deg * RAD, 0.0, doy)
            assert 2.27 < z.hydrostatic < 2.34
            assert 0.01 < z.wet < 0.30


def test_figure_4_unb3_and_unb3m_wet_delays():
    """UNB3 is the same eq. (5) with Table 1's water vapour pressure in place of
    Table 2's humidity; Figure 4: the difference is almost zero at 15 degrees and
    below and up to about 4 cm near 60 degrees, around day 100 or 300."""

    def difference(lat_deg, doy):
        lat = lat_deg * RAD
        m = unb3m.meteorology(lat, doy)
        e_unb3 = interp(WVP_AVG, lat_deg) - interp(WVP_AMP, lat_deg) * math.cos(2 * math.pi * (doy - 28) / 365.25)
        unb3 = unb3m.zenith_from(lat, 0.0, m._replace(vapour_pressure=e_unb3))
        new = unb3m.zenith_from(lat, 0.0, m)
        assert unb3.hydrostatic == new.hydrostatic   # "the hydrostatic delay will not be affected"
        return unb3.wet - new.wet

    for lat_deg in (0.0, 10.0, 15.0):
        assert max(abs(difference(lat_deg, d)) for d in range(1, 366)) < 0.005
    peak = {lat: max(range(1, 366), key=lambda d, lat=lat: abs(difference(lat, d))) for lat in (30.0, 45.0, 60.0, 75.0)}
    size = {lat: abs(difference(lat, d)) for lat, d in peak.items()}
    assert max(size, key=size.get) == 60.0
    assert 0.03 < size[60.0] < 0.05
    assert min(abs(peak[60.0] - 100.0), abs(peak[60.0] - 300.0)) <= 30


def test_stated_accuracy():
    # Table 7, UNB3m total delay minus ray trace, all years: mean -0.5 cm, std 4.9 cm.
    assert unb3m.TOTAL_SIGMA == 0.049 and unb3m.TOTAL_MEAN == -0.005


# ---- refusals -----------------------------------------------------------------------------


def test_height_beyond_the_profile_and_unphysical_parameters():
    lat = 45.0 * RAD
    with pytest.raises(KernelError) as e:
        unb3m.zenith_delays(lat, 60_000.0, 100.0)
    assert (e.value.kernel, e.value.status) == ("unb3m", k.HEIGHT)
    m = unb3m.meteorology(lat, 100.0)
    for bad in (m._replace(temperature=0.0), m._replace(lapse_rate=0.0), m._replace(vapour_lapse=-1.0)):
        with pytest.raises(KernelError) as e:
            unb3m.zenith_from(lat, 0.0, bad)
        assert (e.value.kernel, e.value.status) == ("unb3m", k.PARAMETERS)
    assert k.OK == 0 and k.HEIGHT > 0 and k.PARAMETERS > 0 and k.HEIGHT != k.PARAMETERS
