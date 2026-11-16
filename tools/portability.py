"""The kernels on microcontrollers: builds, sizes, symbols, and one query in C.

    uv run python tools/portability.py report        sizes and needed routines per target
    uv run python tools/portability.py query         the host build of tools/query, its output
    uv run python tools/portability.py time [N]      host seconds per query and per epoch
    uv run python tools/portability.py count         double operations per part, kind, kernel
    uv run python tools/portability.py archive DIR   the ESP32-S3 archive and object of a query

Targets (``TARGETS``): the host C compiler, and each cross compiler found:
Arm Cortex-M4 with its single-precision FPU (``arm-none-eabi-gcc``, the
flags of DESIGN 4), Espressif's ESP32-S3 (Xtensa LX7,
single-precision FPU) and ESP32 (Xtensa LX6, single-precision FPU), and
ESP32-C3 and ESP32-C6 (RISC-V, no FPU). The Espressif compilers are looked
for on the PATH and under ``~/.espressif/tools``; they are only run.

Every kernel is compiled unchanged with the strict flags of the host build
(``STRICT``: C99, pedantic, every warning an error, no contraction, no fast
math). For an ELF target ``size -A`` gives each object's code (``.text``,
``.literal``) and constant data (``.rodata``, ``.srodata``) and ``nm -u``
the routines it needs: the math library's, and the compiler's run-time
routines, which on a target without a double-precision FPU are the
software double arithmetic (``__adddf3``, ``__muldf3``, ``__divdf3``,
comparisons, conversions). tests/test_portability.py runs the same builds
as the gate. [design: the kernels' portability is measured, not assumed
(DESIGN 4, Q17); a tool, so that the tables of DESIGN 4 can be rebuilt]
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import sysconfig
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
KERNEL_SOURCES = ROOT / "kernels" / "src"
KERNEL_INCLUDE = ROOT / "kernels" / "include"
QUERY = ROOT / "tools" / "query"

# [design: the flags the kernels are held to everywhere: those of setup.py's
# host build (C99, -O2, no contraction, no fast math) with every warning of
# DESIGN 4, as errors]
STRICT = ["-std=c99", "-pedantic", "-Wall", "-Wextra", "-Wconversion", "-Wdouble-promotion", "-Werror",
          "-ffp-contract=off", "-fno-fast-math", "-O2"]

# [spec: C99 7.12, the functions of <math.h> the kernels may call]
MATH_LIBRARY = {
    "acos", "asin", "atan", "atan2", "cos", "sin", "tan", "cosh", "sinh", "tanh", "exp", "log", "log10",
    "pow", "sqrt", "cbrt", "hypot", "fabs", "floor", "ceil", "fmod", "round", "trunc", "lround", "modf",
    "frexp", "ldexp", "copysign", "fmin", "fmax", "remainder", "isfinite", "__fpclassifyd", "__isfinited",
    "finite",
}
# [spec: the GCC internals manual, "The GCC low-level runtime library", soft
# float routines: names ending in df2 or df3 (binary64) or converting to or
# from df]
SOFT_DOUBLE = re.compile(r"^__\w*df\w*$")
# [spec: GCC manual, "Language Standards Supported by GCC": "GCC requires the
# freestanding environment provide memcpy, memmove, memset and memcmp"; it
# emits calls to them for copies of arrays and structures]
COMPILER_REQUIRED = {"memcpy", "memmove", "memset", "memcmp"}


@dataclass(frozen=True)
class Target:
    """[design: one compiler and its flags; ``elf`` when its binutils can size it]"""

    name: str
    what: str
    compiler: str | None
    flags: tuple[str, ...]
    elf: bool = True
    absent: str = ""

    def tool(self, name: str) -> str:
        """[design: the binutils beside the compiler, ``<prefix>-size``, ``<prefix>-nm``]"""
        assert self.compiler is not None
        return re.sub(r"gcc$", name, self.compiler)


def _espressif(triple: str, program: str) -> str | None:
    """[design: an Espressif compiler on the PATH, else the newest under ~/.espressif/tools]"""
    found = shutil.which(program)
    if found:
        return found
    home = Path(os.environ.get("IDF_TOOLS_PATH", Path.home() / ".espressif")) / "tools" / triple
    candidates = sorted(home.glob(f"*/{triple}/bin/{program}"))
    return str(candidates[-1]) if candidates else None


def _host_compiler() -> str:
    """[design: the compiler setuptools builds the extensions with on Linux and macOS ($CC,
    else the interpreter's sysconfig CC, else cc), so that the C program and the engine
    are compiled alike and can agree to the last bit]"""
    return (os.environ.get("CC") or sysconfig.get_config_var("CC") or "cc").split()[0]


def targets() -> list[Target]:
    """[design: the host and every microcontroller target of DESIGN 4]"""
    xtensa, riscv = "xtensa-esp-elf", "riscv32-esp-elf"
    s3 = _espressif(xtensa, "xtensa-esp32s3-elf-gcc")
    esp32 = _espressif(xtensa, "xtensa-esp32-elf-gcc")
    rv = _espressif(riscv, "riscv32-esp-elf-gcc")
    arm = shutil.which("arm-none-eabi-gcc")
    host = shutil.which(_host_compiler())
    return [
        Target("host", f"the host C compiler ({_host_compiler()})", host, (), elf=False,
               absent=f"no host C compiler {_host_compiler()!r} on the PATH"),
        Target("cortex-m4", "Arm Cortex-M4, single-precision FPU (arm-none-eabi-gcc)", arm,
               ("-mcpu=cortex-m4", "-mthumb", "-mfpu=fpv4-sp-d16", "-mfloat-abi=hard"),
               absent="arm-none-eabi-gcc is not installed"),
        Target("esp32s3", "ESP32-S3, Xtensa LX7, single-precision FPU (xtensa-esp32s3-elf-gcc)", s3,
               ("-mlongcalls",), absent="xtensa-esp32s3-elf-gcc is not installed (Espressif toolchain)"),
        Target("esp32", "ESP32, Xtensa LX6, single-precision FPU (xtensa-esp32-elf-gcc)", esp32,
               ("-mlongcalls",), absent="xtensa-esp32-elf-gcc is not installed (Espressif toolchain)"),
        Target("esp32c3", "ESP32-C3, RISC-V RV32IMC, no FPU (riscv32-esp-elf-gcc)", rv,
               ("-march=rv32imc_zicsr_zifencei", "-mabi=ilp32"),
               absent="riscv32-esp-elf-gcc is not installed (Espressif toolchain)"),
        Target("esp32c6", "ESP32-C6, RISC-V RV32IMAC, no FPU (riscv32-esp-elf-gcc)", rv,
               ("-march=rv32imac_zicsr_zifencei", "-mabi=ilp32"),
               absent="riscv32-esp-elf-gcc is not installed (Espressif toolchain)"),
    ]


def kernels() -> list[str]:
    """[design: every kernel of kernels/src, by name]"""
    return sorted(p.stem for p in KERNEL_SOURCES.glob("*.c"))


def run(command: list[str], **kw) -> subprocess.CompletedProcess:
    """[design: a build step; its output kept for the failure message]"""
    return subprocess.run(command, capture_output=True, text=True, **kw)


def compile_kernel(target: Target, kernel: str, out: Path, optimise: str = "-O2") -> subprocess.CompletedProcess:
    """[design: one kernel, unchanged, to one object with STRICT and the target's flags]"""
    flags = [f for f in STRICT if not f.startswith("-O")] + [optimise]
    return run([target.compiler, *flags, *target.flags, f"-I{KERNEL_INCLUDE}", "-c",
                str(KERNEL_SOURCES / f"{kernel}.c"), "-o", str(out)])


def included_headers(target: Target, kernel: str) -> list[tuple[str, str]]:
    """(includer, included) of every header a kernel's source and kernel headers include,
    from the compiler's own trace (-H). [design: what the preprocessor really opens]"""
    r = run([target.compiler, *STRICT, *target.flags, f"-I{KERNEL_INCLUDE}", "-H", "-fsyntax-only",
             str(KERNEL_SOURCES / f"{kernel}.c")])
    stack = [str(KERNEL_SOURCES / f"{kernel}.c")]
    pairs = []
    for line in r.stderr.splitlines():
        m = re.match(r"^(\.+) (.+)$", line)
        if not m:
            continue
        depth, header = len(m[1]), m[2].strip()
        if not header.endswith(".h"):
            continue  # Apple clang also lists the SDK's settings file, which is not included
        del stack[depth:]
        pairs.append((stack[depth - 1], header))
        stack.append(header)
    return pairs


@dataclass
class Footprint:
    """[design: what one object costs: code and constant bytes, mutable bytes, needed symbols]"""

    code: int = 0
    constant: int = 0
    mutable: int = 0
    undefined: list[str] = field(default_factory=list)


def footprint(target: Target, obj: Path) -> Footprint:
    """[design: ``size -A`` sections and ``nm -u`` of an ELF object]"""
    f = Footprint()
    for line in run([target.tool("size"), "-A", str(obj)]).stdout.splitlines():
        parts = line.split()
        if len(parts) < 2 or not parts[1].isdigit():
            continue
        section, size = parts[0], int(parts[1])
        if section.startswith((".text", ".literal")):
            f.code += size
        elif section.startswith((".rodata", ".srodata")):
            f.constant += size
        elif section.startswith((".data", ".sdata", ".bss", ".sbss")):
            f.mutable += size
    nm = run([target.tool("nm"), "-u", str(obj)]).stdout.split()
    f.undefined = sorted({s for s in nm if s not in ("U", "w")})
    return f


def classify(symbols: set[str]) -> tuple[list[str], list[str], list[str], list[str]]:
    """(math library, software double, other compiler run-time, anything else) of undefined
    symbols. [design: MATH_LIBRARY, SOFT_DOUBLE; other ``__`` names are libgcc's, and
    COMPILER_REQUIRED are the compiler's own calls]"""
    math_ = sorted(s for s in symbols if s in MATH_LIBRARY)
    soft = sorted(s for s in symbols if SOFT_DOUBLE.match(s) and s not in MATH_LIBRARY)
    runtime = sorted(s for s in symbols if (s.startswith("__") or s in COMPILER_REQUIRED)
                     and s not in math_ and s not in soft)
    other = sorted(symbols - set(math_) - set(soft) - set(runtime))
    return math_, soft, runtime, other


def report(directory: Path, optimise: str = "-O2") -> str:
    """[design: the table of DESIGN 4, Measured on microcontrollers: per target and kernel the
    code and constant bytes, the total, and the routines the objects need]"""
    lines = []
    names = kernels()
    for t in targets():
        if not t.elf:
            continue
        if t.compiler is None:
            lines.append(f"{t.name}: skipped, {t.absent}")
            continue
        out = directory / t.name
        out.mkdir(parents=True, exist_ok=True)
        lines.append(f"{t.name}: {t.what}, {optimise}")
        lines.append(f"  {'kernel':<12}{'code':>8}{'const':>8}")
        total_code = total_const = 0
        needed: set[str] = set()
        users: dict[str, list[str]] = {}
        for k in names:
            obj = out / f"{k}.o"
            r = compile_kernel(t, k, obj, optimise)
            if r.returncode != 0:
                lines.append(f"  {k:<12} FAILED: {r.stderr.strip()}")
                continue
            f = footprint(t, obj)
            total_code += f.code
            total_const += f.constant
            needed |= set(f.undefined)
            for s in f.undefined:
                users.setdefault(s, []).append(k)
            lines.append(f"  {k:<12}{f.code:>8}{f.constant:>8}")
        lines.append(f"  {'total':<12}{total_code:>8}{total_const:>8}")
        math_, soft, runtime, other = classify(needed)
        lines.append(f"  math library ({len(math_)}): {' '.join(math_)}")
        lines.append(f"  software double ({len(soft)}): {' '.join(soft) or 'none'}")
        lines.append(f"  other run-time ({len(runtime)}): "
                     + (" ".join(f"{s} ({', '.join(users[s])})" for s in runtime) or "none"))
        if other:
            lines.append(f"  other ({len(other)}): {' '.join(other)}")
    return "\n".join(lines)


# ---- the query program ------------------------------------------------------------------

QUERY_SOURCES = ("query.c", "scene.c")


def build_query(directory: Path) -> Path:
    """[design: the host build of tools/query with STRICT and the math library, linked last]"""
    exe = directory / "one_query"
    sources = [str(QUERY / "main.c")] + [str(QUERY / s) for s in QUERY_SOURCES] \
        + [str(KERNEL_SOURCES / f"{k}.c") for k in kernels()]
    r = run([_host_compiler(), *STRICT, f"-I{KERNEL_INCLUDE}", f"-I{QUERY}", "-o", str(exe), *sources, "-lm"])
    if r.returncode != 0:
        raise RuntimeError(f"host build of tools/query failed:\n{r.stderr}")
    return exe


def query_values(exe: Path) -> dict[str, float]:
    """[design: main.c's "name value" lines, as Python floats (17 digits are exact)]"""
    r = run([str(exe)])
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    return {name: float(value) for name, value in (line.split() for line in r.stdout.splitlines())}


def time_query(exe: Path, n: int = 20000) -> dict[str, float]:
    """[design: main.c's ``time N``, the best of five runs of N, seconds]"""
    r = run([str(exe), "time", str(n)])
    return {name: float(value) for name, value in (line.split() for line in r.stdout.splitlines())}


def cxx_compiler() -> str | None:
    """[design: $CXX, else c++, for the counting build]"""
    return shutil.which(os.environ.get("CXX", "c++").split()[0])


def build_count(directory: Path) -> Path:
    """[design: count.hpp: every source as C++ with ``double`` counted, one unit per file]"""
    cxx = cxx_compiler()
    if cxx is None:
        raise RuntimeError("no C++ compiler for the counting build")
    objects = []
    units = [(KERNEL_SOURCES / f"{k}.c", k) for k in kernels()] \
        + [(QUERY / s, "query") for s in QUERY_SOURCES] + [(QUERY / "count.cpp", "count")]
    for source, unit in units:
        obj = directory / f"{source.stem}.count.o"
        r = run([cxx, "-x", "c++", "-std=c++17", "-O2", "-ffp-contract=off", "-fno-fast-math", "-w",
                 "-include", str(QUERY / "count.hpp"), f"-DSSRQ_UNIT=\"{unit}\"",
                 f"-I{KERNEL_INCLUDE}", f"-I{QUERY}", "-c", str(source), "-o", str(obj)])
        if r.returncode != 0:
            raise RuntimeError(f"counting build of {source.name} failed:\n{r.stderr}")
        objects.append(str(obj))
    exe = directory / "count_query"
    r = run([cxx, "-o", str(exe), *objects, "-lm"])
    if r.returncode != 0:
        raise RuntimeError(f"counting link failed:\n{r.stderr}")
    return exe


def counts(exe: Path) -> tuple[dict[tuple[str, str, str], int], dict[str, float]]:
    """[design: count.cpp's lines: (part, unit, kind) to count, and the values]"""
    r = run([str(exe)])
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    n, values = {}, {}
    for line in r.stdout.splitlines():
        parts = line.split()
        if parts[0] == "count":
            n[(parts[1], parts[2], parts[3])] = int(parts[4])
        else:
            values[parts[1]] = float(parts[2])
    return n, values


KIND_ORDER = ("add", "mul", "div", "sqrt", "cmp", "neg", "conv", "sin", "cos", "tan", "asin", "acos", "atan",
              "atan2", "exp", "log", "pow", "hypot", "fabs", "floor", "ceil", "fmod", "round", "isfinite")


def count_table(n: dict[tuple[str, str, str], int]) -> str:
    """[design: per part, the operations by kind, then by kernel]"""
    lines = []
    for part in ("receiver", "choice", "epoch", "query"):
        kinds: dict[str, int] = {}
        units: dict[str, dict[str, int]] = {}
        for (p, unit, kind), v in n.items():
            if p == part:
                kinds[kind] = kinds.get(kind, 0) + v
                units.setdefault(unit, {})[kind] = units.get(unit, {}).get(kind, 0) + v
        lines.append(f"{part}: " + ", ".join(f"{k} {kinds[k]}" for k in KIND_ORDER if k in kinds))
        for unit in sorted(units):
            u = units[unit]
            lines.append(f"  {unit:<12}" + ", ".join(f"{k} {u[k]}" for k in KIND_ORDER if k in u))
    return "\n".join(lines)


def build_archive(target: Target, directory: Path) -> dict[str, object]:
    """The query and the kernels for a target, as far as its toolchain goes alone: the
    objects, a static archive, and a relocatable object with the math library and the
    compiler's run-time routines it needs linked in (``-r``), so that its size is what a
    firmware image gains. [design: a full link needs the vendor's start-up code and linker
    script (a board support package), which ssrcore does not take]"""
    directory.mkdir(parents=True, exist_ok=True)
    objects = []
    for source in [QUERY / s for s in QUERY_SOURCES] + [KERNEL_SOURCES / f"{k}.c" for k in kernels()]:
        obj = directory / f"{source.stem}.o"
        r = run([target.compiler, *STRICT, *target.flags, f"-I{KERNEL_INCLUDE}", f"-I{QUERY}", "-c",
                 str(source), "-o", str(obj)])
        if r.returncode != 0:
            raise RuntimeError(f"{target.name} build of {source.name} failed:\n{r.stderr}")
        objects.append(obj)
    archive = directory / "libssrcore_query.a"
    archive.unlink(missing_ok=True)
    r = run([target.tool("ar"), "rcs", str(archive), *map(str, objects)])
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    linked = directory / "ssrcore_query_linked.o"
    r = run([target.compiler, *target.flags, "-r", "-nostdlib", "-o", str(linked),
             "-Wl,--whole-archive", str(archive), "-Wl,--no-whole-archive", "-lm", "-lgcc"])
    if r.returncode != 0:
        raise RuntimeError(f"relocatable link failed:\n{r.stderr}")
    own = Footprint()
    needed: set[str] = set()
    for obj in objects:
        f = footprint(target, obj)
        own.code, own.constant, own.mutable = own.code + f.code, own.constant + f.constant, own.mutable + f.mutable
        needed |= set(f.undefined)
    whole = footprint(target, linked)
    # the math library and run-time routines alone, as the same link pulls them in
    routines = directory / "routines_linked.o"
    forced = [f"-Wl,-u,{s}" for s in sorted(needed) if not s.startswith("ssr")]
    r = run([target.compiler, *target.flags, "-r", "-nostdlib", "-o", str(routines), *forced, "-lm", "-lgcc"])
    runtime = footprint(target, routines) if r.returncode == 0 else None
    return {"archive": archive, "linked": linked, "own": own, "whole": whole, "routines": runtime,
            "needed": sorted(needed)}


def main(argv: list[str]) -> int:
    """[design: the module docstring]"""
    command = argv[1] if len(argv) > 1 else "report"
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        if command == "report":
            print(report(directory))
            print()
            print(report(directory / "Os", "-Os"))
        elif command == "query":
            for name, value in query_values(build_query(directory)).items():
                print(f"{name} {value!r}")
        elif command == "time":
            n = int(argv[2]) if len(argv) > 2 else 20000
            for name, value in time_query(build_query(directory), n).items():
                print(f"{name} {value * 1e6:.3f} us")
        elif command == "count":
            n, _ = counts(build_count(directory))
            print(count_table(n))
        elif command == "archive":
            out = Path(argv[2]) if len(argv) > 2 else directory
            s3 = next(t for t in targets() if t.name == "esp32s3")
            if s3.compiler is None:
                print(f"esp32s3: skipped, {s3.absent}")
                return 0
            a = build_archive(s3, out)
            own, whole = a["own"], a["whole"]
            print(f"archive {a['archive']}")
            print(f"query and kernels: code {own.code}, constant {own.constant}, mutable {own.mutable} bytes"
                  " (objects, before the link relaxes their long calls)")
            print(f"with libm and libgcc linked in (-r): code {whole.code}, constant {whole.constant}, "
                  f"mutable {whole.mutable} bytes")
            if a["routines"] is not None:
                rt = a["routines"]
                print(f"of which the math library and run-time routines alone: code {rt.code}, "
                      f"constant {rt.constant} bytes")
            print(f"still undefined: {' '.join(whole.undefined) or 'none'}")
        else:
            print(__doc__)
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
