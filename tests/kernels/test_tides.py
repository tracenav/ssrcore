"""The tides kernel: IERS Conventions (2010), TN36 section 7.1.1.

- Worked examples: the test cases published with the IERS Conventions
  software for 7.1.1 (TN36 footnote 1, DEHANTTIDEINEL.F header, cases of
  2009, 2012 and 2015), within 0.1 mm per component. Their Sun and Moon are
  given; the times are 0h UTC, here with UT1 = UTC (|UT1 - UTC| < 0.9 s
  moves the step 2 terms by under 0.001 mm). The software sums the full
  Table 7.3a series of its own; TN36 prints the terms >= 0.05 mm, which is
  what the kernel uses.
- Eq. (5) at three geometries where it reduces to a number by hand.
- The magnitudes TN36 7.1.1 states for each part, with the Moon and the Sun
  at their mean distances: latitude dependence not exceeding 0.4 mm radial
  and 0.2 mm transverse; degree 3 reaching 1.7 mm radial, transverse not
  exceeding 0.2 mm; the l(1) terms up to 0.8 mm (diurnal) and 1.0 mm
  (semidiurnal); "can reach" or "up to" is checked within 0.1 mm.
- The permanent deformation of 7.1.1.2: about -12 cm at the poles and +6 cm
  at the equator, radially.
"""

import math

import pytest

from ssrcore._kernels import astro, tides
from ssrcore.toolkit import astro as toolkit_astro
from ssrcore.toolkit import tides as toolkit_tides

MS, MM = tides.SUN_MASS_RATIO, tides.MOON_MASS_RATIO
RE = 6378136.6
MOON_MEAN, SUN_MEAN = 384400e3, 149597870700.0
GPS_EPOCH_JD = 2444244.5

# (site, sun, moon, (year, month, day) at 0h UTC, TAI - UTC, expected displacement)
IERS_CASES = [
    ((4075578.385, 931852.890, 4801570.154),
     (137859926952.015, 54228127881.4350, 23509422341.6960),
     (-179996231.920342, -312468450.131567, -169288918.592160),
     (2009, 4, 13), 34,
     (0.07700420357108125891, 0.06304056321824967613, 0.05516568152597246810)),
    ((1112189.660, -4842955.026, 3985352.284),
     (-54537460436.2357, 130244288385.279, 56463429031.5996),
     (300396716.912, 243238281.451, 120548075.939),
     (2012, 7, 13), 35,
     (-0.2036831479592075833e-01, 0.5658254776225972449e-01, -0.7597679676871742227e-01)),
    ((1112200.5696, -4842957.8511, 3985345.9122),
     (100210282451.6279, 103055630398.3160, 56855096480.4475),
     (369817604.4348, 1897917.5258, 120804980.8284),
     (2015, 7, 15), 36,
     (0.00509570869172363845, 0.0828663025983528700, -0.0636634925404189617)),
]


def gps_of_utc_midnight(y: int, mo: int, d: int, tai_utc: int) -> tuple[int, float]:
    a = (14 - mo) // 12
    yy, mm = y + 4800 - a, mo + 12 * a - 3
    jdn = d + (153 * mm + 2) // 5 + 365 * yy + yy // 4 - yy // 100 + yy // 400 - 32045
    seconds = (jdn - 0.5 - GPS_EPOCH_JD) * 86400.0 + (tai_utc - 19)
    week = int(seconds // 604800)
    return week, seconds - week * 604800


@pytest.mark.parametrize("site,sun,moon,date,tai_utc,expected", IERS_CASES)
def test_iers_worked_examples(site, sun, moon, date, tai_utc, expected):
    week, sow = gps_of_utc_midnight(*date, tai_utc)
    t = astro.tt_centuries(week, sow)
    gmst = astro.gmst(week, sow, -(tai_utc - 19.0))
    got = toolkit_tides.displacement(site, sun, moon, t, gmst)
    for g, e in zip(got, expected):
        assert abs(g - e) < 0.0001


def radial_transverse(site, v):
    n = math.dist(site, (0, 0, 0))
    radial = sum(a * b for a, b in zip(site, v)) / n
    return radial, math.sqrt(max(0.0, sum(a * a for a in v) - radial * radial))


def test_degree2_by_hand():
    site = (RE, 0.0, 0.0)                      # equator: geocentric latitude 0
    k = MM * RE**4 / MOON_MEAN**3
    # body at the zenith: h2 K radially (P2 = 1), nothing across
    r, tr = radial_transverse(site, tides.degree2(site, (MOON_MEAN, 0, 0), MM, 0.6078, 0.0, 0.0847, 0.0))
    assert abs(r - 0.6078 * k) < 1e-12 and tr < 1e-12
    # on the horizon: -h2 K / 2
    r, tr = radial_transverse(site, tides.degree2(site, (0, MOON_MEAN, 0), MM, 0.6078, 0.0, 0.0847, 0.0))
    assert abs(r + 0.6078 * k / 2) < 1e-12 and tr < 1e-12
    # at 45 degrees: h2 K / 4 radially, 3 l2 K (1/sqrt 2)(1/sqrt 2) = 1.5 l2 K across
    c = MOON_MEAN / math.sqrt(2)
    r, tr = radial_transverse(site, tides.degree2(site, (c, c, 0), MM, 0.6078, 0.0, 0.0847, 0.0))
    assert abs(r - 0.6078 * k / 4) < 1e-12 and abs(tr - 1.5 * 0.0847 * k) < 1e-12
    # with the latitude terms at the equator, P2(sin 0) = -1/2: h = 0.6078 + 0.0003
    r, _ = radial_transverse(site, tides.degree2(site, (MOON_MEAN, 0, 0), MM, 0.6078, -0.0006, 0.0847, 0.0002))
    assert abs(r - 0.6081 * k) < 1e-12


def scan(fn, distance, max_dec):
    """Largest radial and transverse parts over sites and body directions."""
    best_r = best_t = 0.0
    for lat in range(-90, 91, 5):
        la = math.radians(lat)
        site = (RE * math.cos(la), 0.0, RE * math.sin(la))
        for dec in range(-max_dec, max_dec + 1, 5):
            de = math.radians(dec)
            for ha in range(0, 360, 5):
                h = math.radians(ha)
                body = (distance * math.cos(de) * math.cos(h), distance * math.cos(de) * math.sin(h),
                        distance * math.sin(de))
                r, t = radial_transverse(site, fn(site, body))
                best_r, best_t = max(best_r, abs(r)), max(best_t, t)
    return best_r * 1000.0, best_t * 1000.0     # mm


MOON_DEC, SUN_DEC = 30, 25    # bounds of the declinations, degrees


def test_documented_latitude_dependence():
    total_r = total_t = 0.0
    for ratio, distance, dec in ((MM, MOON_MEAN, MOON_DEC), (MS, SUN_MEAN, SUN_DEC)):
        def part(site, body):
            a = tides.degree2(site, body, ratio, 0.6078, -0.0006, 0.0847, 0.0002)
            b = tides.degree2(site, body, ratio, 0.6078, 0.0, 0.0847, 0.0)
            return [x - y for x, y in zip(a, b)]
        r, t = scan(part, distance, dec)
        total_r, total_t = total_r + r, total_t + t
    assert total_r <= 0.4 and total_t <= 0.2


def test_documented_degree3():
    r, t = scan(lambda s, b: tides.degree3(s, b, MM), MOON_MEAN, 90)
    assert abs(r - 1.7) <= 0.1 and t <= 0.2
    # the Sun's degree 3 is negligible, which is why the kernel leaves it out
    r, t = scan(lambda s, b: tides.degree3(s, b, MS), SUN_MEAN, 90)
    assert r < 0.01 and t < 0.01


@pytest.mark.parametrize("diurnal,semidiurnal,documented", [(0.0012, 0.0, 0.8), (0.0, 0.0024, 1.0)])
def test_documented_l1(diurnal, semidiurnal, documented):
    total = 0.0
    for ratio, distance, dec in ((MM, MOON_MEAN, MOON_DEC), (MS, SUN_MEAN, SUN_DEC)):
        r, t = scan(lambda s, b: tides.l1(s, b, ratio, diurnal, semidiurnal), distance, dec)
        assert r < 1e-9
        total += t
    assert abs(total - documented) <= 0.1


def test_out_of_phase_is_small():
    # hI, lI of a few 1e-3 against h2 = 0.6: under 1 % of the degree 2 tide
    r, t = scan(lambda s, b: tides.out_of_phase(s, b, MM, -0.0025, -0.0007, -0.0022, -0.0007),
                MOON_MEAN, MOON_DEC)
    assert 0.0 < r < 1.0 and 0.0 < t < 1.0


@pytest.mark.parametrize("lat,radial", [(90, -0.1205), (-90, -0.1205), (0, 0.060325)])
def test_permanent_tide(lat, radial):
    la = math.radians(lat)
    site = (RE * math.cos(la), 0.0, RE * math.sin(la))
    r, t = radial_transverse(site, toolkit_tides.permanent_tide(site))
    assert abs(r - radial) < 1e-6 and t < 1e-9


def test_permanent_tide_transverse_at_45():
    site = (RE * math.sqrt(0.5), 0.0, RE * math.sqrt(0.5))
    north = (-math.sqrt(0.5), 0.0, math.sqrt(0.5))
    d = tides.permanent(site)
    p2 = 0.25
    assert abs(sum(a * b for a, b in zip(d, north)) - (-0.0252 - 0.0001 * p2)) < 1e-12


def test_step2_bounded_by_its_tables():
    # |dr| <= sum of the |Delta R| of Tables 7.3a and 7.3b, in mm
    bound = (12.00 + 0.78 + 1.73 + 0.12 + 1.23 + 0.07 + 0.51 + 0.50 + 0.01 + 0.22 + 0.01
             + 0.11 + 0.01 + 0.10 + 0.08 + 0.06 + 0.06
             + 0.47 + 0.16 + 0.20 + 0.11 + 0.11 + 0.09 + 0.13 + 0.15 + 0.05 + 0.06)
    biggest = 0.0
    for hour in range(0, 48):
        t = hour / 24.0 / 36525.0
        for lat in (-45, 0, 45):
            la = math.radians(lat)
            site = (RE * math.cos(la), 0.0, RE * math.sin(la))
            r, _ = radial_transverse(site, tides.step2(site, t, astro.gmst(1042, 561600.0 + hour * 3600.0, 0.0)))
            biggest = max(biggest, abs(r) * 1000.0)
    assert 10.0 < biggest <= bound


def test_toolkit_end_to_end():
    site = (4075578.385, 931852.890, 4801570.154)
    d = toolkit_tides.solid_earth_tide(site, toolkit_astro.sky(2310, 345600.0, -18.0))
    assert 0.0 < math.dist(d, (0, 0, 0)) < 0.5
