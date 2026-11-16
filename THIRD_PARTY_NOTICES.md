# Third-party notices

The code of ssrcore is under the MIT licence (`LICENSE`); the items below are
numerical values taken from published works.

## Compiled into the package

Truncated tables in `kernels/src/astro.c`, compiled into every wheel:

- **ELP 2000-82B lunar theory**: 315 of the 37 872 terms of ELP1 to ELP36,
  every term of at least 0.2″ or 0.2 km, and the coefficients of its mean
  arguments. M. Chapront-Touzé and J. Chapront, "The lunar ephemeris
  ELP 2000", Astron. Astrophys. 124, 50 (1983); "ELP 2000-85: a
  semi-analytical lunar ephemeris adequate for historical times", Astron.
  Astrophys. 190, 342 (1988). Distributed by IMCCE – Observatoire de Paris
  (https://ftp.imcce.fr/pub/ephem/moon/elp82b/) and by CDS as catalogue VI/79.
- **VSOP87, version D**: 69 of the 2 425 terms of the Earth file VSOP87D.ear,
  every term of at least 1e-6 rad or au. P. Bretagnon and G. Francou,
  "Planetary theories in rectangular and spherical variables. VSOP 87
  solutions", Astron. Astrophys. 202, 309 (1988). Distributed by IMCCE –
  Observatoire de Paris (https://ftp.imcce.fr/pub/ephem/planets/vsop87/) and
  by CDS as catalogue VI/81.

The distributions state no licence or conditions of use; the tables are
included with this attribution.

## In tests and examples only (source distribution, never a wheel)

- Test cases of DEHANTTIDEINEL.F, software of the IERS Conventions (2010),
  IERS Technical Note 36 (https://iers-conventions.obspm.fr/):
  `tests/kernels/test_tides.py`, `tools/tide_check.py`.
- IMCCE geocentric Moon positions, MOON.Y05
  (https://ftp.imcce.fr/pub/ephem/moon/geopos/): `tests/kernels/test_astro.py`.
- US Naval Observatory season times of 2005 (https://aa.usno.navy.mil/):
  `tests/kernels/test_astro.py`.
- IGS broadcast ephemerides (via BKG,
  https://igs.bkg.bund.de/root_ftp/IGS/BRDC/) and GFZ MGEX rapid orbits
  (ftp://ftp.gfz-potsdam.de/pub/GNSS/products/mgex/), both of 2024-01-01:
  `tests/kernels/test_kepler.py`, `test_glonass.py`,
  `tests/test_joins.py`, `test_engine.py`, `examples/`, `tools/query/scene.c`.
  Openly available under the IGS data policy (https://igs.org/data-access/),
  with acknowledgement of the IGS and its contributing agencies: data provided
  by the International GNSS Service, BKG Frankfurt and GFZ Potsdam.
