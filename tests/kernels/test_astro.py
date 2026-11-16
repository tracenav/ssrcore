"""The astro kernel against published positions (DESIGN Q11).

- VSOP87: the check values of the distribution file ``vsop87.chk`` (VSOP87D,
  Earth, JD 2451545.0), which the truncated series must meet within its
  stated accuracy.
- ELP 2000-82B: the apparent geocentric positions of the Moon for 2005
  published by IMCCE (G. Francou, ftp.imcce.fr/pub/ephem/moon/geopos,
  file MOON.Y05; dynamical time; computed from the complete ELP 2000-82B).
  Distance is frame-free; right ascension and declination are on the true
  equator and equinox of date, so they also test the nutation of TN36 ch. 5.
  The apparent place differs from the geometric by the Moon's light time,
  under 1 arcsec.
- Earth rotation: ERA at J2000.0 is theta_0 of TN36 ch. 1 Table 1.1; GMST at
  J2000.0 UT1 is checked against a second formulation, GMST1982 (Aoki et al.
  1982, 24110.54841 s at 0h UT1 of 2000-01-01), which TN36 eq. (5.32) offsets
  by 0.014506 arcsec.
- Sun in the terrestrial frame: at the equinoxes of 2005 published by the
  US Naval Observatory (Earth's Seasons, aa.usno.navy.mil), the apparent Sun
  crosses the equator; the geometric Sun is then 20.5 arcsec of annual
  aberration ahead along the ecliptic, i.e. at declination +-20.5" sin(eps)
  = +-8.2 arcsec.
"""

import math

import pytest

from ssrcore._kernels import astro
from ssrcore.toolkit import astro as toolkit_astro

ARCSEC = math.pi / 648000.0
GPS_EPOCH_JD = 2444244.5


def gps_of_tt_jd(jd_tt: float) -> tuple[int, float]:
    seconds = (jd_tt - GPS_EPOCH_JD) * 86400.0 - 51.184
    week = int(seconds // 604800)
    return week, seconds - week * 604800


def gps_of_utc(y: int, mo: int, d: int, h: int, mi: int, leap: int) -> tuple[int, float]:
    a = (14 - mo) // 12
    yy, mm = y + 4800 - a, mo + 12 * a - 3
    jdn = d + (153 * mm + 2) // 5 + 365 * yy + yy // 4 - yy // 100 + yy // 400 - 32045
    seconds = (jdn - 0.5 - GPS_EPOCH_JD) * 86400.0 + h * 3600 + mi * 60 + leap
    week = int(seconds // 604800)
    return week, seconds - week * 604800


def test_tt_centuries_at_j2000():
    # J2000.0 TT is 2000-01-01 11:58:55.816 GPS: week 1042, 561548.816 s
    assert abs(astro.tt_centuries(1042, 561548.816)) < 1e-15
    week, sow = gps_of_tt_jd(2451545.0 + 36525.0)
    assert abs(astro.tt_centuries(week, sow) - 1.0) < 1e-12


def test_vsop87_check_values_at_j2000():
    # vsop87.chk: VSOP87D EARTH JD2451545.0: l 1.7519238681 rad, b -.0000039656 rad, r .9833276819 au
    lon, lat, r = astro.earth_heliocentric(0.0)
    assert abs(lon - 1.7519238681) < 2.0 * ARCSEC
    assert abs(lat - -0.0000039656) < 2.0 * ARCSEC
    assert abs(r - 0.9833276819) * 149597870.7 < 1000.0      # km


def test_sun_is_the_earth_reversed():
    for t in (-0.1, 0.0, 0.123, 0.5):
        lon, lat, r = astro.earth_heliocentric(t)
        slon, slat, sr = astro.sun_ecliptic(t)
        assert abs(math.remainder(slon - lon - math.pi, 2 * math.pi)) < 1e-12
        assert slat == -lat
        assert sr == r * 149597870700.0


# IMCCE geopos, MOON.Y05: (month, day, hour TT, RA h m s, Dec sign d m s, distance km)
MOON_2005 = [
    (1, 7, 12, (16, 9, 49.757), (-1, 24, 3, 25.54), 364560.3858),
    (1, 15, 0, (23, 52, 36.084), (-1, 3, 36, 24.10), 374141.5371),
    (2, 7, 12, (20, 6, 43.946), (-1, 25, 23, 44.83), 358731.9666),
    (2, 15, 0, (2, 47, 14.453), (1, 17, 58, 21.00), 390812.3767),
    (3, 15, 0, (3, 18, 39.145), (1, 20, 58, 38.92), 390138.4454),
    (4, 7, 12, (0, 4, 19.537), (-1, 1, 43, 58.66), 372536.6366),
    (4, 15, 0, (6, 33, 57.519), (1, 28, 18, 53.47), 402267.7941),
    (5, 15, 0, (8, 52, 27.056), (1, 22, 44, 8.54), 404485.8954),
    (6, 7, 12, (5, 30, 37.616), (1, 27, 45, 51.39), 399628.0556),
    (7, 15, 0, (13, 38, 30.760), (-1, 11, 6, 19.56), 388087.0752),
    (8, 15, 0, (16, 57, 54.557), (-1, 27, 10, 55.76), 371727.4862),
    (9, 7, 12, (13, 37, 35.381), (-1, 11, 25, 4.00), 393421.7710),
    (9, 15, 0, (20, 51, 56.927), (-1, 22, 31, 10.36), 362377.9888),
    (10, 15, 0, (23, 19, 51.195), (-1, 6, 48, 7.05), 365553.2591),
    (11, 15, 0, (2, 26, 3.990), (1, 16, 54, 21.66), 377848.5469),
    (12, 15, 0, (4, 53, 10.315), (1, 27, 6, 9.26), 390458.0384),
]


def julian_day(y: int, mo: int, d: int, hour: float) -> float:
    a = (14 - mo) // 12
    yy, mm = y + 4800 - a, mo + 12 * a - 3
    jdn = d + (153 * mm + 2) // 5 + 365 * yy + yy // 4 - yy // 100 + yy // 400 - 32045
    return jdn - 0.5 + hour / 24.0


@pytest.mark.parametrize("month,day,hour,ra,dec,distance", MOON_2005)
def test_moon_against_imcce_2005(month, day, hour, ra, dec, distance):
    t = (julian_day(2005, month, day, hour) - 2451545.0) / 36525.0
    lbr = astro.moon_ecliptic(t)
    assert abs(lbr[2] / 1000.0 - distance) < 3.5                      # km; truncation 3.3 km
    x, y, z = toolkit_astro.true_equatorial(lbr, t)
    r = math.sqrt(x * x + y * y + z * z)
    ra_expected = (ra[0] + ra[1] / 60 + ra[2] / 3600) * 15 * math.pi / 180
    dec_expected = dec[0] * (dec[1] + dec[2] / 60 + dec[3] / 3600) * math.pi / 180
    d_ra = math.remainder(math.atan2(y, x) - ra_expected, 2 * math.pi) * math.cos(dec_expected)
    d_dec = math.asin(z / r) - dec_expected
    assert abs(d_ra) < 4.0 * ARCSEC
    assert abs(d_dec) < 4.0 * ARCSEC


def test_earth_rotation_angle_at_j2000_is_theta0():
    # JD 2451545.0 UT1 is GPS second 630763200 on the UT1 scale (week 1042, 561600 s).
    assert abs(astro.earth_rotation_angle(1042, 561600.0, 0.0) - 2 * math.pi * 0.7790572732640) < 1e-12


def test_gmst_at_j2000_against_gmst1982():
    gmst1982 = (24110.54841 + 43200.0) / 86400.0 * 2 * math.pi
    t = astro.tt_centuries(1042, 561600.0)
    expected = gmst1982 + (0.014506 + 4612.156534 * t) * ARCSEC
    assert abs(astro.gmst(1042, 561600.0, 0.0) - expected) < 0.001 * ARCSEC


@pytest.mark.parametrize("when,sign", [((2005, 3, 20, 12, 33), 1), ((2005, 9, 22, 22, 23), -1)])
def test_sun_on_the_equator_at_the_usno_equinoxes(when, sign):
    week, sow = gps_of_utc(*when, leap=13)          # GPS - UTC = 13 s in 2005
    sun, _ = astro.sun_moon_ecef(week, sow, -13.0)   # UT1 - UTC < 0.9 s, left out (the kernel)
    dec = math.asin(sun[2] / math.sqrt(sun[0] ** 2 + sun[1] ** 2 + sun[2] ** 2))
    aberration = 20.5 * ARCSEC * math.sin(23.44 * math.pi / 180)
    # the minute of the published time leaves 0.5 arcsec; truncation 2 arcsec
    assert abs(dec - sign * aberration) < 3.0 * ARCSEC


def test_toolkit_sun_and_moon_distances():
    week, sow = gps_of_utc(2024, 6, 1, 0, 0, leap=18)
    s = toolkit_astro.sky(week, sow, -18.0)
    assert 1.471e11 < math.dist(s.sun, (0, 0, 0)) < 1.521e11
    assert 3.56e8 < math.dist(s.moon, (0, 0, 0)) < 4.07e8
    assert (s.sun, s.moon) == astro.sun_moon_ecef(week, sow, -18.0)
    assert s.tt_centuries == astro.tt_centuries(week, sow)
    assert s.gmst == astro.gmst(week, sow, -18.0)
