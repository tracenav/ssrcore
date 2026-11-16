"""Solid Earth tide displacement of a site.

DESIGN.md section 2.7. IERS Conventions (2010), Technical Note 36, section
7.1.1, steps 1 and 2, and 7.1.1.2 for the permanent deformation. The Sun
and the Moon are arguments, an :class:`ssrcore.toolkit.astro.Sky` of the
epoch (one evaluation of the astronomy serves the tide and the wind-up).
The displacement is the conventional tide-free one of TN36 7.1.1,
permanent part included, as the coordinates of an ITRF datum need.
"""

from ssrcore._kernels import tides as _tides
from ssrcore.toolkit import Vector
from ssrcore.toolkit.astro import Sky


def displacement(site: Vector, sun: Vector, moon: Vector, t: float, gmst: float) -> Vector:
    """[spec: TN36 7.1.1 steps 1 and 2] Displacement (ECEF, m) of ``site``
    (ECEF, m) with the geocentric Sun and Moon (ECEF, m), TT ``t`` in Julian
    centuries from J2000.0 and Greenwich mean sidereal time ``gmst`` (rad)."""
    return _tides.displacement(site, sun, moon, t, gmst)


def solid_earth_tide(site: Vector, sky: Sky) -> Vector:
    """[spec: TN36 7.1.1 steps 1 and 2; Sun and Moon of DESIGN Q11] Displacement
    (ECEF, m) of ``site`` at the epoch of ``sky`` (``toolkit.astro.sky``)."""
    return _tides.displacement(site, sky.sun, sky.moon, sky.tt_centuries, sky.gmst)


def permanent_tide(site: Vector) -> Vector:
    """[spec: TN36 7.1.1.2 eqs. (14a), (14b)] The permanent deformation (ECEF,
    m), to add to a conventional tide-free position for a mean-tide one."""
    return _tides.permanent(site)
