"""The troposphere kernel and toolkit module (DESIGN 2.7).

Checked against: the Niell coefficients of TM-23 5.4.2.2.1 Tables 5.2 and
5.3 at their table nodes; the properties TM-23 eqs. (5.67) to (5.69) imply
(unity at zenith, linear height term); an independent mapping function
(Black and Eisner, TM-23 eq. (5.61)) and the 1/sin E limit; SPARTN 8.12 on
synthetic coefficients (the polynomial and grid are tests/kernels/test_area.py's); the UNB3m table (Table 1, TM-23 Table 5.1) and, for
the height reduction of DESIGN Q7, the pressures of the International
Standard Atmosphere (ICAO Doc 7488, ISO 2533) and the exponential model of
TM-23 eq. (5.66).
"""

import math

import pytest

from ssrcore._kernels import troposphere as k
from ssrcore.toolkit import area
from ssrcore.toolkit import troposphere as tropo

RAD = math.pi / 180.0
HALF_YEAR = 365.25 / 2.0

# TM-23 Table 5.2 (hydrostatic average and amplitude) and Table 5.3 (wet), at 15..75 degrees.
HYD_AVG = [
    [1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3],
    [2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3],
    [62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3],
]
HYD_AMP = [
    [0.0, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5],
    [0.0, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5],
    [0.0, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5],
]
WET = [
    [5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4],
    [1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3],
    [4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2],
]
NODES = [15.0, 30.0, 45.0, 60.0, 75.0]


# ---- Niell coefficients, TM-23 5.4.2.2.1 ---------------------------------------------------

@pytest.mark.parametrize("j", range(5))
def test_niell_coefficients_at_table_nodes(j):
    lat = NODES[j] * RAD
    # doy 28: cos = 1, xi = avg - amp; half a year later: avg + amp (eq. 5.70).
    h, w = tropo.niell_coefficients(lat, 28.0)
    for i in range(3):
        assert h[i] == pytest.approx(HYD_AVG[i][j] - HYD_AMP[i][j], rel=1e-14, abs=1e-20)
        assert w[i] == pytest.approx(WET[i][j], rel=1e-15)
    h, _ = tropo.niell_coefficients(lat, 28.0 + HALF_YEAR)
    for i in range(3):
        assert h[i] == pytest.approx(HYD_AVG[i][j] + HYD_AMP[i][j], rel=1e-14, abs=1e-20)
    # A quarter year from T0, the seasonal term vanishes.
    h, _ = tropo.niell_coefficients(lat, 28.0 + 365.25 / 4.0)
    for i in range(3):
        assert h[i] == pytest.approx(HYD_AVG[i][j], rel=1e-12)


def test_niell_coefficients_interpolate_linearly_and_hold_at_the_ends():
    h, w = tropo.niell_coefficients(37.5 * RAD, 28.0 + 365.25 / 4.0)
    for i in range(3):
        assert h[i] == pytest.approx((HYD_AVG[i][1] + HYD_AVG[i][2]) / 2.0, rel=1e-12)
        assert w[i] == pytest.approx((WET[i][1] + WET[i][2]) / 2.0, rel=1e-14)
    assert tropo.niell_coefficients(5.0 * RAD, 100.0) == tropo.niell_coefficients(15.0 * RAD, 100.0)
    assert tropo.niell_coefficients(-89.0 * RAD, 100.0) == tropo.niell_coefficients(-80.0 * RAD, 100.0)


def test_niell_southern_hemisphere_is_half_a_year_out_of_phase():
    for lat in (20.0, 40.0, 70.0):
        for doy in (1.0, 28.0, 150.5, 300.0):
            south = tropo.niell_coefficients(-lat * RAD, doy)
            north = tropo.niell_coefficients(lat * RAD, doy + HALF_YEAR)
            for a, b in zip(south.hydrostatic + south.wet, north.hydrostatic + north.wet):
                assert a == pytest.approx(b, rel=1e-13)


# ---- Niell mapping ----------------------------------------------------------------------

def test_marini_is_one_at_zenith():
    for a, b, c in ((1.2e-3, 2.9e-3, 62.6e-3), (5.8e-4, 1.4e-3, 4.3e-2), (2.53e-5, 5.49e-3, 1.14e-3)):
        assert tropo.marini(1.0, a, b, c) == pytest.approx(1.0, abs=1e-15)


@pytest.mark.parametrize("lat", [-60.0, -10.0, 0.0, 30.0, 52.0, 80.0])
@pytest.mark.parametrize("height", [0.0, 500.0, 2500.0])
def test_mapping_is_one_at_zenith_at_any_height(lat, height):
    m = tropo.niell_mapping(lat * RAD, height, 123.4, math.pi / 2.0)
    assert m.hydrostatic == pytest.approx(1.0, abs=1e-14)
    assert m.wet == pytest.approx(1.0, abs=1e-14)


def test_height_term_is_linear_in_height_and_absent_from_wet():
    lat, doy, el = 45.0 * RAD, 200.0, 7.0 * RAD
    m0 = tropo.niell_mapping(lat, 0.0, doy, el)
    m1 = tropo.niell_mapping(lat, 1000.0, doy, el)
    m3 = tropo.niell_mapping(lat, 3000.0, doy, el)
    assert m1.wet == m0.wet == m3.wet
    assert m3.hydrostatic - m0.hydrostatic == pytest.approx(3.0 * (m1.hydrostatic - m0.hydrostatic), rel=1e-12)
    # eq. (5.67): the term per km is 1/sin E - m(E, a_ht, b_ht, c_ht), positive and small.
    per_km = 1.0 / math.sin(el) - tropo.marini(math.sin(el), 2.53e-5, 5.49e-3, 1.14e-3)
    assert m1.hydrostatic - m0.hydrostatic == pytest.approx(per_km, rel=1e-12)
    assert 0.0 < per_km < 0.05


@pytest.mark.parametrize("el_deg", [15.0, 20.0, 30.0, 45.0, 60.0])
def test_mapping_agrees_with_black_and_eisner(el_deg):
    # An independent mapping function, TM-23 eq. (5.61): M = 1.001 / sqrt(0.002001 + sin^2 E).
    s = math.sin(el_deg * RAD)
    black_eisner = 1.001 / math.sqrt(0.002001 + s * s)
    m = tropo.niell_mapping(45.0 * RAD, 0.0, 100.0, el_deg * RAD)
    assert m.hydrostatic == pytest.approx(black_eisner, rel=1e-2)
    assert m.wet == pytest.approx(black_eisner, rel=1e-2)


def test_mapping_tends_to_cosecant_and_wet_exceeds_hydrostatic_low_down():
    for el_deg in (3.0, 5.0, 10.0, 20.0, 40.0, 70.0):
        el = el_deg * RAD
        m = tropo.niell_mapping(30.0 * RAD, 0.0, 28.0, el)
        csc = 1.0 / math.sin(el)
        # A flat atmosphere gives 1/sin E; curvature makes the true value smaller.
        assert m.hydrostatic < csc and m.wet < csc
        assert m.wet > m.hydrostatic  # the wet layer is thinner, closer to the flat case
        if el_deg >= 20.0:
            assert m.hydrostatic == pytest.approx(csc, rel=1e-2)


@pytest.mark.parametrize("el", [0.0, -0.1, math.pi / 2.0 + 1e-9, float("nan")])
def test_mapping_refuses_elevations_outside_the_sky(el):
    status, _, _ = k.niell_mapping(0.5, 0.0, 10.0, el)
    assert status == k.ELEVATION
    with pytest.raises(ValueError):
        tropo.niell_mapping(0.5, 0.0, 10.0, el)


# ---- SPARTN 8.12 ------------------------------------------------------------------------

T = (0.123, -0.004, 0.0025, 0.0003)  # synthetic T00, T01, T10, T11


def test_slant_on_synthetic_coefficients():
    # SF043 encoded 100: 100 * 0.004 - 0.508 = -0.108 m, offset from 2.3 m (Table 6.2).
    t_h = 2.3 + (100 * 0.004 - 0.508)
    lat, lon_offset, h, doy = 47.0 * RAD, -0.75, 410.0, 250.0
    t_p = area.polynomial(2, T, 0.5, lon_offset)
    t_g = area.interpolate((0.4, 0.3, 0.2, 0.1), (0.01, 0.02, -0.01, 0.0))
    # At zenith both mappings are 1: the slant is the zenith total.
    mz = tropo.niell_mapping(lat, h, doy, math.pi / 2.0)
    assert tropo.slant(t_h, t_p, t_g, *mz) == pytest.approx(t_h + t_p + t_g, abs=1e-14)
    # At 12 degrees, each part scales with its own mapping (8.12).
    m = tropo.niell_mapping(lat, h, doy, 12.0 * RAD)
    slant = tropo.slant(t_h, t_p, t_g, m.hydrostatic, m.wet)
    assert slant == pytest.approx(t_h * m.hydrostatic + (t_p + t_g) * m.wet, abs=1e-13)
    assert tropo.slant(t_h, 0.0, 0.0, m.hydrostatic, m.wet) == pytest.approx(t_h * m.hydrostatic)
    assert tropo.slant(0.0, t_p, t_g, m.hydrostatic, m.wet) == pytest.approx((t_p + t_g) * m.wet)
    # Plausibility: about 2.4 m at zenith gives about 11 m at 12 degrees.
    assert 10.0 < slant < 12.5


# ---- UNB3m, DESIGN Q7 -------------------------------------------------------------------

@pytest.mark.parametrize("j,t_avg,t_amp,b_avg,b_amp", [
    (0, 299.65, 0.00, 6.30e-3, 0.00),
    (2, 283.15, 11.00, 5.58e-3, 0.32e-3),
    (4, 263.65, 14.50, 4.53e-3, 0.62e-3),
])
def test_unb3m_temperature_at_table_nodes(j, t_avg, t_amp, b_avg, b_amp):
    lat = NODES[j] * RAD
    t, b = tropo.unb3m_temperature(lat, 28.0)
    assert t == pytest.approx(t_avg - t_amp, rel=1e-14)
    assert b == pytest.approx(b_avg - b_amp, rel=1e-13)
    # TM-23 eq. (5.62): D_min = 211 in the south.
    assert tropo.unb3m_temperature(-lat, 211.0) == pytest.approx((t, b), rel=1e-14)


def test_unb3m_temperature_interpolates_between_nodes():
    t, b = tropo.unb3m_temperature(52.5 * RAD, 28.0 + 365.25 / 4.0)
    assert t == pytest.approx((283.15 + 272.15) / 2.0, rel=1e-12)
    assert b == pytest.approx((5.58e-3 + 5.39e-3) / 2.0, rel=1e-12)


# International Standard Atmosphere (ICAO Doc 7488; ISO 2533): T0 = 288.15 K, lapse 6.5 K/km;
# published pressures, hPa, at geopotential heights 1, 2, 3 km; 1013.25 hPa at 0 m.
ISA = [(1000.0, 898.746), (2000.0, 794.952), (3000.0, 701.085)]


@pytest.mark.parametrize("height,pressure", ISA)
def test_height_factor_follows_the_standard_atmosphere(height, pressure):
    # UNB3m eq. (4) reduced to a ratio: pressure ratio of the lapse-rate profile, times g_m(0)/g_m(h).
    lat = 45.0 * RAD
    gm0 = 9.784 * (1.0 - 2.66e-3 * math.cos(2.0 * lat))
    gmh = 9.784 * (1.0 - 2.66e-3 * math.cos(2.0 * lat) - 2.8e-7 * height)
    f = tropo.unb3m_height_factor(lat, 288.15, 0.0065, height)
    assert f == pytest.approx(pressure / 1013.25 * gm0 / gmh, rel=2e-5)


def test_height_factor_against_the_exponential_model():
    # TM-23 eq. (5.66): T_dry = 2.3 exp(-0.116e-3 H). Independent model; agreement to 2 %
    # of the factor up to 1 km (5 cm of a 2.3 m delay): the seasonal temperature of UNB3m
    # matters (a cold column decays faster, 1.8 % off at 60 degrees in winter), which the
    # global model lacks. The isothermal test below is the tighter check.
    for lat_deg in (-45.0, 15.0, 45.0, 60.0):
        for doy in (28.0, 200.0):
            for height in (300.0, 1000.0):
                f = tropo.hydrostatic_height_factor(lat_deg * RAD, height, doy)
                assert f == pytest.approx(math.exp(-0.116e-3 * height), rel=2e-2)
    assert tropo.hydrostatic_height_factor(0.7, 0.0, 99.0) == 1.0


def test_height_factor_against_the_isothermal_scale_height():
    # Hydrostatic equilibrium of an isothermal layer at T0: ln f = -g h / (R_d T0). The lapse rate
    # changes this by a relative beta h / (2 T0) to first order; g_m(h) adds 2.8e-7 h.
    for lat_deg in (-60.0, 15.0, 45.0, 75.0):
        for doy in (28.0, 120.0, 211.0):
            t0, beta = tropo.unb3m_temperature(lat_deg * RAD, doy)
            for height in (100.0, 500.0, 1000.0):
                f = tropo.hydrostatic_height_factor(lat_deg * RAD, height, doy)
                isothermal = -9.80665 * height / (287.054 * t0)
                assert math.log(f) == pytest.approx(isothermal, rel=beta * height / t0)


def test_height_factor_refuses_heights_beyond_the_profile():
    status, _ = k.unb3m_height_factor(0.5, 288.15, 0.0065, 288.15 / 0.0065)
    assert status == k.HEIGHT
    with pytest.raises(ValueError):
        tropo.unb3m_height_factor(0.5, 288.15, 0.0065, 50000.0)
