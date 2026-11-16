"""A second evaluation of the solid Earth tide, apart from the `tides` kernel.

    uv run python tools/tide_check.py refs/Q11_sun_moon_tides

Displacement of a site by the solid Earth tide, IERS TN36 7.1.1 steps 1 and
2, conventional tide-free (TN36 7.1.1.2), written from the TN36 text in
Python. The Sun and Moon are those of the complete ELP 2000-82B and VSOP87D
series: `astro.c` rebuilt in a temporary directory with every term of the
sources (tools/astro_tables.py with no truncation), the only part shared with
the engine being the summation of a series. They are rotated to the
Earth-fixed frame here by an equinox-based path (mean obliquity TN36 eq.
5.40, nutation in longitude and obliquity, GAST = GMST of TN36 eqs. 5.15 and
5.32 plus the equation of the equinoxes), where the engine takes a
CIO-based path; polar motion is left out (under 0.5 arcsec).

It prints, for the DEHANTTIDEINEL header cases (their Sun and Moon given),
the record's and this evaluation's difference from the published values;
then the largest and the rms difference between the record's displacement
(`toolkit.tides.solid_earth_tide` with `toolkit.astro.sky`) and this one
over a grid of sites and times (`grid`, GPS − UTC 18 s). So it checks the
truncation of the series in `astro.c`, the engine's Earth-fixed rotation
and the `tides` kernel together.

[design: an evaluation that shares neither the truncated series nor the
kernel checks both; a tool, not a test, because it needs the uncommitted
sources of refs/ and a C compiler]
"""

from __future__ import annotations

import ctypes
import datetime
import math
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
import astro_tables  # noqa: E402

from ssrcore.toolkit import astro, tides  # noqa: E402

AS = math.pi / 648000.0
# [spec: TN36 ch. 1 Table 1.1: GM_Moon / GM_Earth, GM_Sun, GM_Earth, a_E]
MU_MOON = 0.0123000371
MU_SUN = 1.32712442099e20 / 3.986004418e14
RE = 6378136.6
J2000 = datetime.datetime(2000, 1, 1, 12)
GPS0 = datetime.datetime(1980, 1, 6)
LEAP = 18

# [spec: DEHANTTIDEINEL.F header (TN36 footnote 1): site, Sun, Moon (m, ECEF),
#  date at 0h UTC, TAI − UTC (s), displacement (m)]
CASES = [
    ((4075578.385, 931852.890, 4801570.154), (137859926952.015, 54228127881.4350, 23509422341.6960),
     (-179996231.920342, -312468450.131567, -169288918.592160), datetime.datetime(2009, 4, 13), 34,
     (0.7700420357108125891e-01, 0.6304056321824967613e-01, 0.5516568152597246810e-01)),
    ((1112189.660, -4842955.026, 3985352.284), (-54537460436.2357, 130244288385.279, 56463429031.5996),
     (300396716.912, 243238281.451, 120548075.939), datetime.datetime(2012, 7, 13), 35,
     (-0.2036831479592075833e-01, 0.5658254776225972449e-01, -0.7597679676871742227e-01)),
    ((1112200.5696, -4842957.8511, 3985345.9122), (100210282451.6279, 103055630398.3160, 56855096480.4475),
     (369817604.4348, 1897917.5258, 120804980.8284), datetime.datetime(2015, 7, 15), 36,
     (0.00509570869172363845, 0.0828663025983528700, -0.0636634925404189617)),
]

# [spec: TN36 Table 7.3a, the terms it prints: Doodson multipliers of τ, s, h,
#  p, N', ps; ΔR in and out of phase; ΔT in and out of phase (mm)]
DIURNAL = [
    ((1, -2, 0, 1, 0, 0), -0.08, 0.00, -0.01, 0.01),
    ((1, -1, 0, 0, -1, 0), -0.10, 0.00, 0.00, 0.00),
    ((1, -1, 0, 0, 0, 0), -0.51, 0.00, -0.02, 0.03),
    ((1, 0, 0, 1, 0, 0), 0.06, 0.00, 0.00, 0.00),
    ((1, 1, -3, 0, 0, 1), -0.06, 0.00, 0.00, 0.00),
    ((1, 1, -2, 0, 0, 0), -1.23, -0.07, 0.06, 0.01),
    ((1, 1, 0, 0, -1, 0), -0.22, 0.01, 0.01, 0.00),
    ((1, 1, 0, 0, 0, 0), 12.00, -0.78, -0.67, -0.03),
    ((1, 1, 0, 0, 1, 0), 1.73, -0.12, -0.10, 0.00),
    ((1, 1, 1, 0, 0, -1), -0.50, -0.01, 0.03, 0.00),
    ((1, 1, 2, 0, 0, 0), -0.11, 0.01, 0.01, 0.00),
]
# [spec: TN36 Table 7.3b, as Table 7.3a]
LONG = [
    ((0, 0, 0, 0, 1, 0), 0.47, 0.16, 0.23, 0.07),
    ((0, 0, 2, 0, 0, 0), -0.20, -0.11, -0.12, -0.05),
    ((0, 1, 0, -1, 0, 0), -0.11, -0.09, -0.08, -0.04),
    ((0, 2, 0, 0, 0, 0), -0.13, -0.15, -0.11, -0.07),
    ((0, 2, 0, 0, 1, 0), -0.05, -0.06, -0.05, -0.03),
]


def build(sources: Path, directory: Path) -> ctypes.CDLL:
    """[design: astro.c with every term of its sources, compiled as the kernels
    are (C99, no contraction), loaded with ctypes]"""
    astro_tables.ELP_MIN = 0.0
    astro_tables.VSOP_MIN = 0.0
    tables = astro_tables.generate(sources)
    text = astro_tables.TARGET.read_text()
    head, rest = text.split(astro_tables.BEGIN, 1)
    _, tail = rest.split(astro_tables.END, 1)
    src = directory / "astro_full.c"
    src.write_text(head + astro_tables.BEGIN + tables + astro_tables.END + tail)
    lib = directory / "astro_full.so"
    subprocess.run(["cc", "-O2", "-std=c99", "-ffp-contract=off", "-shared", "-fPIC",
                    "-I", str(ROOT / "kernels" / "include"), str(src), "-o", str(lib), "-lm"], check=True)
    out = ctypes.CDLL(str(lib))
    vec = ctypes.POINTER(ctypes.c_double * 3)
    out.ssr_astro_moon_ecliptic.argtypes = [ctypes.c_double, vec]
    out.ssr_astro_sun_ecliptic.argtypes = [ctypes.c_double, vec]
    out.ssr_astro_nutation.argtypes = [ctypes.c_double, ctypes.POINTER(ctypes.c_double),
                                       ctypes.POINTER(ctypes.c_double)]
    return out


def mean_obliquity(t: float) -> float:
    """[spec: TN36 eq. (5.40)] Mean obliquity of the ecliptic (rad), t in TT centuries."""
    return (84381.406 - 46.836769 * t - 0.0001831 * t * t + 0.00200340 * t ** 3
            - 0.000000576 * t ** 4 - 0.0000000434 * t ** 5) * AS


def gmst(tu_days: float, t: float) -> float:
    """[spec: TN36 eqs. (5.15), (5.32)] GMST (rad); tu_days = Julian UT1 date − 2451545.0."""
    frac = tu_days % 1.0
    era = 2 * math.pi * ((frac + 0.7790572732640 + 0.00273781191135448 * tu_days) % 1.0)
    return era + (0.014506 + 4612.156534 * t + 1.3915817 * t * t - 0.00000044 * t ** 3
                  - 0.000029956 * t ** 4 - 0.0000000368 * t ** 5) * AS


def _ecliptic(fn, t: float) -> tuple[float, float, float]:
    """[design: longitude, latitude (rad) and distance (m) of date from the series]"""
    out = (ctypes.c_double * 3)()
    fn(t, ctypes.byref(out))
    return tuple(out)


def to_earth_fixed(lbr, gast: float, dpsi: float, eps: float) -> tuple[float, float, float]:
    """[spec: TN36 5.4.4, the equinox-based transformation: true ecliptic of date
    (longitude plus Δψ) to the true equator by the true obliquity, then to the
    Earth-fixed frame by GAST; polar motion left out]"""
    lon, lat, r = lbr
    lon += dpsi
    x, y, z = r * math.cos(lat) * math.cos(lon), r * math.cos(lat) * math.sin(lon), r * math.sin(lat)
    ye, ze = y * math.cos(eps) - z * math.sin(eps), y * math.sin(eps) + z * math.cos(eps)
    c, s = math.cos(gast), math.sin(gast)
    return (c * x + s * ye, -s * x + c * ye, ze)


def sun_moon(lib: ctypes.CDLL, tt_days: float, ut1_days: float):
    """[spec: TN36 eq. (5.39) for GAST = GMST + Δψ cos ε] Sun and Moon (m, Earth-fixed),
    t (TT centuries) and GMST; times in days from J2000.0."""
    t = tt_days / 36525.0
    dpsi, deps = ctypes.c_double(), ctypes.c_double()
    lib.ssr_astro_nutation(t, ctypes.byref(dpsi), ctypes.byref(deps))
    eps = mean_obliquity(t) + deps.value
    g = gmst(ut1_days, t)
    gast = g + dpsi.value * math.cos(eps)
    moon = to_earth_fixed(_ecliptic(lib.ssr_astro_moon_ecliptic, t), gast, dpsi.value, eps)
    sun = to_earth_fixed(_ecliptic(lib.ssr_astro_sun_ecliptic, t), gast, dpsi.value, eps)
    return sun, moon, t, g


def _dot(a, b) -> float:
    """[design: vector helper]"""
    return sum(x * y for x, y in zip(a, b))


def _unit(a):
    """[design: vector helper: unit vector and norm]"""
    n = math.sqrt(_dot(a, a))
    return tuple(x / n for x in a), n


def _add(a, b, k: float = 1.0):
    """[design: vector helper: a + k b]"""
    return tuple(x + k * y for x, y in zip(a, b))


def doodson(t: float, gmst_rad: float) -> tuple[float, ...]:
    """[spec: TN36 eq. (5.43) for the Delaunay arguments; TN36 ch. 6, the relations
    s = F + Ω, h = s − D, p = s − l, N' = −Ω, ps = s − D − l', τ = GMST + π − s]"""
    def p(c):
        """[design: arcsec polynomial in t, to radians]"""
        return sum(ci * t ** i for i, ci in enumerate(c)) * AS
    l = p((485868.249036, 1717915923.2178, 31.8792, 0.051635, -0.00024470))
    lp = p((1287104.79305, 129596581.0481, -0.5532, 0.000136, -0.00001149))
    f = p((335779.526232, 1739527262.8478, -12.7512, -0.001037, 0.00000417))
    d = p((1072260.70369, 1602961601.2090, -6.3706, 0.006593, -0.00003169))
    om = p((450160.398036, -6962890.5431, 7.4722, 0.007702, -0.00005939))
    s = f + om
    return (gmst_rad + math.pi - s, s, s - d, s - l, -om, s - d - lp)


def displacement(site, sun, moon, t: float, gmst_rad: float) -> tuple[float, float, float]:
    """[spec: TN36 7.1.1 step 1, eqs. (5), (6), (8) to (11); step 2, eqs. (12),
    (13) with Tables 7.3a, 7.3b] Displacement (m, Earth-fixed), conventional
    tide-free."""
    rhat, r = _unit(site)
    sinphi = site[2] / r
    cosphi = math.hypot(site[0], site[1]) / r
    lam = math.atan2(site[1], site[0])
    east = (-math.sin(lam), math.cos(lam), 0.0)
    north = (-sinphi * math.cos(lam), -sinphi * math.sin(lam), cosphi)
    phi = math.atan2(sinphi, cosphi)
    p2 = (3 * sinphi * sinphi - 1) / 2
    h2 = 0.6078 - 0.0006 * p2
    l2 = 0.0847 + 0.0002 * p2
    d = (0.0, 0.0, 0.0)
    for body, mu in ((moon, MU_MOON), (sun, MU_SUN)):
        unit, rn = _unit(body)
        f = mu * RE ** 4 / rn ** 3
        c = _dot(unit, rhat)
        d = _add(d, rhat, f * h2 * (3 * c * c - 1) / 2)                     # eq. (5)
        d = _add(d, _add(unit, rhat, -c), f * 3 * l2 * c)
        x, y, z = body
        big_phi = math.asin(z / rn)
        lam_j = math.atan2(y, x)
        p21 = 3 * math.sqrt(x * x + y * y) * z / rn ** 2                     # 3 sinΦ cosΦ
        p22 = 3 * (x * x + y * y) / rn ** 2                                  # 3 cos²Φ
        dl = lam - lam_j
        dn = -0.0012 * sinphi * f * p21 * sinphi * math.cos(dl)               # eq. (8)
        de = 0.0012 * sinphi * f * p21 * math.cos(2 * phi) * math.sin(dl)
        dn += -0.5 * 0.0024 * sinphi * cosphi * f * p22 * math.cos(2 * dl)   # eq. (9)
        de += -0.5 * 0.0024 * sinphi * cosphi * f * p22 * sinphi * math.sin(2 * dl)
        s2 = math.sin(2 * big_phi)
        dr = -0.75 * (-0.0025) * f * s2 * math.sin(2 * phi) * math.sin(dl)   # eq. (10)
        dn += -1.5 * (-0.0007) * f * s2 * math.cos(2 * phi) * math.sin(dl)
        de += -1.5 * (-0.0007) * f * s2 * sinphi * math.cos(dl)
        c2 = math.cos(big_phi) ** 2
        dr += -0.75 * (-0.0022) * f * c2 * cosphi ** 2 * math.sin(2 * dl)    # eq. (11)
        dn += 0.75 * (-0.0007) * f * c2 * math.sin(2 * phi) * math.sin(2 * dl)
        de += 0.75 * (-0.0007) * f * c2 * (-2 * cosphi) * math.cos(2 * dl)
        d = _add(_add(_add(d, rhat, dr), north, dn), east, de)
    unit, rn = _unit(moon)                                                    # eq. (6), the Moon
    f3 = MU_MOON * RE ** 5 / rn ** 4
    c = _dot(unit, rhat)
    d = _add(d, rhat, f3 * 0.292 * (2.5 * c ** 3 - 1.5 * c))
    d = _add(d, _add(unit, rhat, -c), f3 * 0.015 * (7.5 * c * c - 1.5))
    beta = doodson(t, gmst_rad)
    dr = dn = de = 0.0
    for n, rip, rop, tip, top in DIURNAL:                                     # eq. (12)
        th = _dot(n, beta) + lam
        dr += (rip * math.sin(th) + rop * math.cos(th)) * math.sin(2 * phi)
        de += (tip * math.cos(th) - top * math.sin(th)) * sinphi
        dn += (tip * math.sin(th) + top * math.cos(th)) * math.cos(2 * phi)
    for n, rip, rop, tip, top in LONG:                                        # eq. (13)
        th = _dot(n, beta)
        dr += (1.5 * sinphi ** 2 - 0.5) * (rip * math.cos(th) + rop * math.sin(th))
        dn += (tip * math.cos(th) + top * math.sin(th)) * math.sin(2 * phi)
    return _add(_add(_add(d, rhat, dr * 1e-3), north, dn * 1e-3), east, de * 1e-3)


def header_cases() -> None:
    """[spec: DEHANTTIDEINEL.F header] The record's and this evaluation's
    difference from the published displacements (mm)."""
    for site, sun, moon, when, tai_utc, expected in CASES:
        ut = (when - J2000).total_seconds() / 86400.0
        t = (ut + (tai_utc + 32.184) / 86400.0) / 36525.0
        mine = displacement(site, sun, moon, t, gmst(ut, t))
        gps = when + datetime.timedelta(seconds=tai_utc - 19)
        week = (gps - GPS0).days // 7
        sec = (gps - GPS0 - datetime.timedelta(weeks=week)).total_seconds()
        ours = tides._tides.displacement(site, sun, moon, astro.tt_centuries(week, sec),
                                         astro.gmst(week, sec, -(tai_utc - 19.0)))
        print(f"{when:%Y-%m-%d}: record − published "
              f"{[round((a - b) * 1e3, 4) for a, b in zip(ours, expected)]} mm, "
              f"this − published {[round((a - b) * 1e3, 4) for a, b in zip(mine, expected)]} mm")


def grid():
    """[design: sites every 30° of latitude from −75° to 75° and every 45° of
    longitude, on the ellipsoid, every 3 h of 2024-01-01 and 2024-07-01: both
    seasons, both hemispheres and every phase of the diurnal and semidiurnal
    tides] (GPS week, seconds of week, site ECEF, label) of each point."""
    from ssrcore.toolkit import geodesy
    for day in (datetime.datetime(2024, 1, 1), datetime.datetime(2024, 7, 1)):
        for hour in range(0, 24, 3):
            gps = day + datetime.timedelta(hours=hour, seconds=LEAP) - GPS0
            week, sec = gps.days // 7, (gps - datetime.timedelta(weeks=gps.days // 7)).total_seconds()
            for lat in range(-75, 76, 30):
                for lon in range(0, 360, 45):
                    site = geodesy.geodetic_to_ecef((math.radians(lat), math.radians(lon), 0.0))
                    yield week, sec, tuple(site), f"{day:%Y-%m-%d} {hour:02d}h {lat:+d} {lon:d}"


def compare_points(lib: ctypes.CDLL, points) -> None:
    """[design: the largest and the rms difference over the points, each
    (week, seconds, site, label), GPS − UTC taken as LEAP]
    The record's displacement against this one at every distinct point."""
    seen: set = set()
    worst, at, sq = 0.0, None, 0.0
    for week, sec, site, label in points:
        if (week, sec, site) in seen:
            continue
        seen.add((week, sec, site))
        ours = tides.solid_earth_tide(site, astro.sky(week, sec, -float(LEAP)))
        gps = GPS0 + datetime.timedelta(weeks=week, seconds=sec)
        ut = (gps - datetime.timedelta(seconds=LEAP) - J2000).total_seconds() / 86400.0
        tt = (gps + datetime.timedelta(seconds=51.184) - J2000).total_seconds() / 86400.0
        sun, moon, t, g = sun_moon(lib, tt, ut)
        mine = displacement(site, sun, moon, t, g)
        sq += sum((a - b) ** 2 for a, b in zip(ours, mine))
        d = max(abs(a - b) for a, b in zip(ours, mine))
        if d > worst:
            worst, at = d, label
    n = len(seen)
    print(f"{n} times and sites: largest |record − this| in a component {worst * 1e3:.4f} mm "
          f"({at}), rms {1e3 * math.sqrt(sq / (3 * n)) if n else 0.0:.4f} mm")


def main(argv: list[str]) -> int:
    """[design: command line of the module docstring]"""
    if len(argv) != 2:
        print(__doc__)
        return 2
    header_cases()
    with tempfile.TemporaryDirectory() as tmp:
        compare_points(build(Path(argv[1]), Path(tmp)), grid())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
