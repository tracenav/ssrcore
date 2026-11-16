"""Conventions the toolkit modules share (DESIGN 3.1; ``ssrcore.toolkit``).

The modules were written apart; these checks hold them to one pattern: a
kernel status other than success becomes a :class:`KernelError` naming the
kernel and carrying the kernel's own code, and the functions of a line of
sight take the satellite before the receiver.
"""

import math

import pytest

from ssrcore._kernels import area as k_area
from ssrcore._kernels import correction as k_correction
from ssrcore._kernels import geodesy as k_geodesy
from ssrcore._kernels import glonass as k_glonass
from ssrcore._kernels import ionosphere as k_ionosphere
from ssrcore._kernels import kepler as k_kepler
from ssrcore._kernels import relativity as k_relativity
from ssrcore._kernels import timescale as k_timescale
from ssrcore._kernels import troposphere as k_troposphere
from ssrcore._kernels import windup as k_windup
from ssrcore.toolkit import (KernelError, area, geodesy, ionosphere, orbits, relativity, ssr, timescale,
                             troposphere, windup)

R_SAT = (15_600e3, 7_540e3, 20_140e3)
R_RCV = (4_027_894.0, 307_045.6, 4_919_474.9)
ORIGIN = (0.0, 0.0, 0.0)
AREA = area.Area(50.0, 10.0, 5, 6, 0.5, 0.75)

FAILURES = [
    ("geodesy", k_geodesy.NO_SOLUTION, lambda: geodesy.ecef_to_geodetic(ORIGIN)),
    ("geodesy", k_geodesy.ZERO_RANGE, lambda: geodesy.line_of_sight(R_RCV, R_RCV)),
    ("area", k_area.INVALID, lambda: area.extent(area.Area(50.0, 10.0, 9, 6, 0.5, 0.75))),
    ("area", k_area.OUTSIDE, lambda: area.grid(AREA, 40.0, 0.0)),
    ("area", k_area.EQUATION_TYPE, lambda: area.polynomial(5, (1.0, 0.0, 0.0, 0.0), 0.0, 0.0)),
    ("ionosphere", k_ionosphere.BAD_FREQUENCY, lambda: ionosphere.delay(1.0, 0.0)),
    ("troposphere", k_troposphere.ELEVATION, lambda: troposphere.niell_mapping(0.8, 0.0, 100.0, -0.1)),
    ("troposphere", k_troposphere.HEIGHT, lambda: troposphere.unb3m_height_factor(0.8, 288.0, 6.5e-3, 50_000.0)),
    ("relativity", k_relativity.THROUGH_CENTRE, lambda: relativity.shapiro((2.6e7, 0.0, 0.0), (-6.4e6, 0.0, 0.0))),
    ("windup", k_windup.DEGENERATE, lambda: windup.satellite_frame(ORIGIN, (1.0, 0.0, 0.0))),
    ("correction", k_correction.DEGENERATE, lambda: ssr.orbit_frame(R_SAT, ORIGIN)),
    ("timescale", k_timescale.RANGE, lambda: timescale.days_from_date(2015, 2, 29)),
    ("timescale", k_timescale.RANGE, lambda: timescale.spartn_resolve_half_day(43200, 0)),
    ("kepler", k_kepler.NOT_CONVERGED, lambda: orbits.eccentric_anomaly(1.0, float("nan"))),
    ("glonass", k_glonass.TOO_FAR, lambda: orbits.glonass_state(
        orbits.GlonassOrbit(2300, 0.0, R_SAT, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), 2301, 0.0)),
]


@pytest.mark.parametrize("kernel, status, call", FAILURES, ids=[f"{f[0]}-{f[1]}" for f in FAILURES])
def test_a_kernel_failure_is_a_kernel_error(kernel, status, call):
    assert status > 0
    with pytest.raises(KernelError) as raised:
        call()
    assert raised.value.kernel == kernel
    assert raised.value.status == status
    assert isinstance(raised.value, ValueError)


def test_valid_outcomes_are_values_not_errors():
    # A special value, a position outside an area and an invalid grid node
    # are data the caller acts on (DESIGN 2.3, rule A7).
    assert ssr.field(52, 0x3F).status is ssr.Status.SPECIAL
    assert area.contains(AREA, 40.0, 0.0) is False
    assert area.interpolate((0.25,) * 4, (1.0, None, 1.0, 1.0)) is None


def test_satellite_comes_before_receiver():
    # Each function of a line of sight, given (satellite, receiver), sees the
    # satellite above the receiver; the geometric range is the same.
    elevation, _ = geodesy.elevation_azimuth(R_SAT, R_RCV)
    assert elevation > 0.0
    unit, _, _ = geodesy.line_of_sight(R_SAT, R_RCV)
    d = [s - r for s, r in zip(R_SAT, R_RCV)]
    assert all(math.isclose(u, c / math.dist(R_SAT, R_RCV), abs_tol=1e-15) for u, c in zip(unit, d))
    assert relativity.shapiro(R_SAT, R_RCV) > 0.0
    assert relativity.sagnac(R_SAT, R_RCV) == pytest.approx(
        geodesy.WGS84_OMEGA_E * (R_SAT[0] * R_RCV[1] - R_SAT[1] * R_RCV[0]) / ssr.SPEED_OF_LIGHT, abs=1e-3)
