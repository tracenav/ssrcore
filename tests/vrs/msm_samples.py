"""Minimal synthetic MSMs for the station and schedule tests (no observations:
those are test_observations.py's). One satellite, one signal, one cell of an MSM4 of each
system, with a GLONASS header where the system needs it (RTCM 10403.3 Table
3.5-78)."""

from __future__ import annotations

from gnsswire import rtcm3

from ssrcore.model.types import System

PREFIX = {System.GPS: 107, System.GLONASS: 108, System.GALILEO: 109, System.BEIDOU: 112, System.QZSS: 111}


def msm4(system: System, satellite: int = 1, epoch_ms: int = 0, station: int = 4095, multiple: int = 1) -> rtcm3.Msm:
    """An MSM4 of one satellite with one signal (signal ID 2), arbitrary DF003 and DF393."""
    glonass = system == System.GLONASS
    return rtcm3.Msm(number=PREFIX[system] * 10 + 4, station=station,
                     epoch=None if glonass else epoch_ms,
                     glonass_day=0 if glonass else None, glonass_epoch=epoch_ms if glonass else None,
                     multiple=multiple, iods=0, reserved=0, clock_steering=1, external_clock=0,
                     smoothing=0, smoothing_interval=0,
                     sat_mask=1 << (64 - satellite), sig_mask=1 << (32 - 2), cell_mask=1,
                     rough_ms=(70,), extended_info=None, rough_mod=(512,), rough_rate=None,
                     fine_range=(100,), fine_phase=(-200,), lock=(5,), half_cycle=(0,), cnr=(40,),
                     fine_rate=None)
