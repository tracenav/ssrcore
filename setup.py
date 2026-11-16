"""Build of the kernel extensions; all other metadata is in pyproject.toml.

Every ``src/ssrcore/_kernels/<name>.pyx`` is one extension module,
``ssrcore._kernels.<name>``, built from that binding and the C kernel
``kernels/src/<name>.c``. Extensions are discovered, not listed, so that
adding a kernel does not edit this file.
"""

from pathlib import Path

from Cython.Build import cythonize
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

ROOT = Path(__file__).parent
KERNEL_INCLUDE = "kernels/include"

# [design: reproducible floating point] The same inputs must give the same
# bits on Linux (gcc or clang, x86-64 or arm64) and macOS (clang, arm64 or
# x86-64), so that results agree across hosts well within 0.1 mm and a
# microcontroller build can be checked against the host.
#   -std=c99                    the kernels' language, also the embedded one
#   -ffp-contract=off           no fused multiply-add: clang fuses a*b+c by
#                               default and gcc does with -O2 on arm64, which
#                               changes the last bit on one host but not on
#                               the other
#   -fno-fast-math              no reassociation, no reciprocal, no
#                               flush-to-zero; IEEE 754 semantics throughout
# x86-64 evaluates double in SSE2 registers, so no x87 excess precision
# arises on any host target. The platform's libm is linked; its
# transcendental functions may differ in the last ulp between glibc and
# Apple's libm, which is far below the 0.1 mm tolerance (docs/DESIGN.md,
# section 4).
CFLAGS = ["-std=c99", "-O2", "-ffp-contract=off", "-fno-fast-math"]

# [design: the same semantics under MSVC, which takes none of the flags
# above] /fp:precise is IEEE 754 semantics without contraction (the default
# from Visual Studio 2022, stated); the C runtime holds the math functions,
# so there is no separate libm to link.
MSVC_CFLAGS = ["/fp:precise"]


class BuildKernels(build_ext):
    """[design: the flags of the compiler setuptools chose] GCC and clang take
    CFLAGS and link libm, as each extension states; MSVC takes MSVC_CFLAGS and
    links no libm."""

    def build_extensions(self):
        """[design: the compiler is known only here, once setuptools has chosen it]"""
        if self.compiler.compiler_type == "msvc":
            for extension in self.extensions:
                extension.extra_compile_args = list(MSVC_CFLAGS)
                extension.libraries = [lib for lib in extension.libraries if lib != "m"]
        super().build_extensions()


def kernel_extensions() -> list[Extension]:
    """[design: discovery] One extension per binding; its kernel has the same name."""
    extensions = []
    for binding in sorted((ROOT / "src/ssrcore/_kernels").glob("*.pyx")):
        name = binding.stem
        extensions.append(Extension(
            f"ssrcore._kernels.{name}",
            sources=[str(binding.relative_to(ROOT)), f"kernels/src/{name}.c"],
            include_dirs=[KERNEL_INCLUDE],
            libraries=["m"],
            extra_compile_args=CFLAGS,
        ))
    return extensions


setup(ext_modules=cythonize(kernel_extensions(), language_level=3),
      cmdclass={"build_ext": BuildKernels})
