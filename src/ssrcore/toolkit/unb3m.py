"""UNB3m: a-priori hydrostatic and wet zenith delays from latitude, height and day of year.

Typed wrappers of ``ssrcore._kernels.unb3m``; no state, no clock. After R.
Leandro, M. Santos, R. B. Langley, "UNB Neutral Atmosphere Models:
Development and Performance" (2006), staged as ``refs/UNB3m_Leandro_2006``:

1. the five meteorological parameters at sea level for the latitude and day
   of year (eqs. (1) to (3) over Table 1, with relative humidity from Table
   2 in place of Table 1's water vapour pressure); temperature T0 and lapse
   rate beta are ``toolkit.troposphere.unb3m_temperature`` (the same Table
   1 rows, held once), pressure, humidity and lambda this module's kernel;
2. the water vapour pressure from the humidity (eqs. (9) to (11));
3. the zenith delays at the height (eqs. (4) to (7)).

The slant delay is m_h d_zh + m_nh d_znh with the Niell mapping functions
(eq. (8)), which are ``toolkit.troposphere.niell_mapping``.

Conventions as ``toolkit.troposphere``: latitude in radians, heights in
metres, ``doy`` counted from January 0.0. UNB3m's H is orthometric height
(text after eq. (5)); the functions take whatever height the caller passes.

Accuracy the publication claims (UNB3m minus ray tracing of 703 711
radiosonde profiles, 223 stations in North America, 1990 to 1996): total
zenith delay mean -0.5 cm, standard deviation 4.9 cm, RMS 4.9 cm (Table 7);
the non-hydrostatic part carries almost all of it (4.8 cm, Table 3), with a
bias growing with height to -6.1 cm above 2000 m (Table 5); the hydrostatic
part 1.2 to 1.9 cm standard deviation with no height trend (Table 6). Not
validated in the southern hemisphere (Conclusions).
"""

from typing import NamedTuple

from ssrcore._kernels import unb3m as _k
from ssrcore.toolkit import KernelError
from ssrcore.toolkit import troposphere as _troposphere

#: [spec: UNB3m Table 7, UNB3m total zenith delay minus ray trace, all years] standard
#: deviation, m; the provider's stated 1-sigma of the model's total zenith delay.
TOTAL_SIGMA = 0.049

#: [spec: UNB3m Table 7, all years] mean of UNB3m minus ray trace, m.
TOTAL_MEAN = -0.005


class Meteorology(NamedTuple):
    """Sea-level parameters of UNB3m at one latitude and day.
    [spec: UNB3m eqs. (1) to (3), (9) to (11), Tables 1 and 2]"""

    pressure: float          # P0, mbar
    temperature: float       # T0, K
    humidity: float          # RH, %
    vapour_pressure: float   # e0, mbar
    lapse_rate: float        # beta, K/m
    vapour_lapse: float      # lambda, unitless


class ZenithDelays(NamedTuple):
    """Hydrostatic d_zh and non-hydrostatic d_znh zenith delays, m. [spec: UNB3m eqs. (4), (5)]"""

    hydrostatic: float
    wet: float

    @property
    def total(self) -> float:
        """[spec: UNB3m eq. (8) at zenith, both mappings 1] d_zh + d_znh, m."""
        return self.hydrostatic + self.wet


def table(latitude: float, doy: float) -> tuple[float, float, float]:
    """[spec: UNB3m eqs. (1) to (3), Tables 1 and 2; TM-23 5.4.2.1 eq. (5.62)] (P0 mbar, RH %, lambda).

    [design: south of the equator the season is shifted by half a year (D_min
    211), as TM-23 shifts the UNB3 table and ``troposphere.unb3m_temperature``
    does; UNB3m itself was validated in the north only]
    """
    return _k.parameters(latitude, doy)


def vapour_pressure(humidity: float, temperature: float, pressure: float) -> float:
    """[spec: UNB3m eqs. (9) to (11)] Water vapour pressure e0 (mbar) from RH (%), T0 (K), P0 (mbar)."""
    return _k.vapour_pressure(humidity, temperature, pressure)


def meteorology(latitude: float, doy: float) -> Meteorology:
    """[spec: UNB3m eqs. (1) to (3), (9) to (11), Tables 1 and 2] The sea-level parameters."""
    temperature, lapse_rate = _troposphere.unb3m_temperature(latitude, doy)
    pressure, humidity, lam = _k.parameters(latitude, doy)
    e0 = _k.vapour_pressure(humidity, temperature, pressure)
    return Meteorology(pressure, temperature, humidity, e0, lapse_rate, lam)


def zenith_from(latitude: float, height: float, met: Meteorology) -> ZenithDelays:
    """[spec: UNB3m eqs. (4) to (7)] Zenith delays at ``height`` (m) for given sea-level parameters.

    Raises KernelError where 1 - beta H / T0 <= 0 or the parameters are not physical.
    """
    status, dh, dw = _k.zenith(latitude, height, met.pressure, met.temperature, met.vapour_pressure,
                               met.lapse_rate, met.vapour_lapse)
    if status == _k.HEIGHT:
        raise KernelError("unb3m", status, f"height {height!r} m is beyond the lapse-rate profile")
    if status != _k.OK:
        raise KernelError("unb3m", status, "temperature, lapse rate or lambda not physical")
    return ZenithDelays(dh, dw)


def zenith_delays(latitude: float, height: float, doy: float) -> ZenithDelays:
    """[spec: UNB3m eqs. (1) to (11)] Hydrostatic and wet zenith delays (m) at a latitude (rad),
    height (m) and day of year. Raises as ``zenith_from``."""
    return zenith_from(latitude, height, meteorology(latitude, doy))
