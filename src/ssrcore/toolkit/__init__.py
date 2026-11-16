"""Stateless public functions over the C kernels: one module per kernel.

DESIGN.md section 2.7. Each module is a thin, typed wrapper of
``ssrcore._kernels.<name>``; none holds state or reads a clock.

Conventions shared by every module (DESIGN 3.1, 4):

- vectors are ``Vector`` tuples (x, y, z), ECEF unless a name says otherwise;
- a function of a line of sight takes the satellite first, then the
  receiver;
- angles are radians, except in ``area``, whose SPARTN fields are degrees;
- a kernel status other than success is raised as :class:`KernelError`;
  an argument the Python layer refuses before calling the kernel (a wrong
  count of values) is a plain ``ValueError``; an outcome that is valid data
  (a field's special value, a position outside an area, a grid with an
  invalid node) is returned as a value, not raised.
"""

from __future__ import annotations

Vector = tuple[float, float, float]


class KernelError(ValueError):
    """[design: DESIGN 3.1; a kernel returns a status, the toolkit is the public
    edge where a failure becomes an exception] A kernel refused its input.

    ``kernel`` names the kernel, ``status`` is its status code (positive,
    ``SSR_<KERNEL>_<CONDITION>``). A subclass of ``ValueError``: the input
    was outside the domain of the function.
    """

    def __init__(self, kernel: str, status: int, message: str) -> None:
        super().__init__(f"{kernel}: {message} (status {status})")
        self.kernel = kernel
        self.status = status
