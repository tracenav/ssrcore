"""The rules of DESIGN.md that can be checked mechanically (DESIGN 5.1).

Checked over the repository as it stands, not over a list of expected files,
so that every kernel and module added later is held to the same rules:

- every C function, in a kernel header and in a kernel source, has a comment
  carrying ``[spec: …]`` or ``[design: …]`` directly before it, and every
  public Python function (``.py`` and Cython ``.pyx``) has one in its
  docstring (DESIGN 5.1, Specification-driven; DESIGN 3.1);
- kernels follow DESIGN 3.1: one header and one source per kernel, includes
  only ``<stdint.h>``, ``<math.h>`` and other kernel headers, exported
  functions named ``ssr_<kernel>_<what>``, status codes 0 for success and
  distinct positive otherwise, explicit integer widths, no ``float``, no
  allocation, no I/O, no ``errno``, no static mutable state;
- each Cython binding declares its own kernel's header and imports nothing
  from ``ssrcore``;
- no module imports upward across the layers of DESIGN 3.1, and only the
  toolkit imports the private ``ssrcore._kernels``;
- the core opens no file or socket and reads no clock (DESIGN 2.1, I/O at
  the edge);
- no ``REVIEW`` marker, a basis still to be found, remains (DESIGN 5.1,
  Specification-driven: there is no third kind), and one fails the gate.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "src" / "ssrcore"
KERNEL_SOURCES = ROOT / "kernels" / "src"
KERNEL_HEADERS = ROOT / "kernels" / "include" / "ssrcore"
BINDINGS = PACKAGE / "_kernels"

# A REVIEW marker fails the gate (DESIGN 5.1, Specification-driven).
REVIEW_MARKERS_FAIL = True

# Directories that are never source: refs/ holds staged specifications and
# build/ local build output and data; the rest are tool state.
NOT_SOURCE = {"refs", "build", "dist", ".git", ".venv", ".pytest_cache", "__pycache__"}

MARKER = re.compile(r"\[(?:spec|design):\s*[^\]\s][^\]]*\]")


def source_files(*suffixes: str, root: Path = ROOT) -> list[Path]:
    """Files of the repository with these suffixes, outside ``NOT_SOURCE``."""
    found = []
    for path in root.rglob("*"):
        rel = path.relative_to(root)
        if path.suffix in suffixes and path.is_file() and not NOT_SOURCE & set(rel.parts) \
                and not any(p.endswith(".egg-info") for p in rel.parts):
            found.append(path)
    return sorted(found)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


# ---- C ----------------------------------------------------------------------------------


def c_items(text: str):
    """Top-level items of a C file as (kind, text, line, preceding comment).

    ``kind`` is ``"code"`` for a declaration or definition head (the text
    before its ``;`` or ``{`` at brace depth 0) or ``"pp"`` for a
    preprocessor line. The preceding comment is the text of the comments
    between the previous item and this one, so a function's own comment is
    the one directly above it. Bodies (brace depth > 0) are returned as
    ``"body"`` items for the checks that look inside functions.
    """
    items = []
    i, n, line = 0, len(text), 1
    comment: list[str] = []
    code, code_line = [], None
    depth, body, body_line = 0, [], None
    while i < n:
        c = text[i]
        if text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            chunk = text[i:j]
            if depth == 0 and not "".join(code).strip():
                comment.append(chunk)
            line += chunk.count("\n")
            i = j
            continue
        if text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            if depth == 0 and not "".join(code).strip():
                comment.append(text[i:j])
            i = j
            continue
        if c in "\"'":
            j = i + 1
            while j < n and text[j] != c:
                j += 2 if text[j] == "\\" else 1
            chunk = text[i:j + 1]
            (body if depth else code).append(chunk)
            i = j + 1
            continue
        if c == "#" and depth == 0 and not "".join(code).strip():
            j = i
            while True:  # a preprocessor line, with its continuations
                k = text.find("\n", j)
                k = n if k < 0 else k
                if k > 0 and text[k - 1] == "\\" and k < n:
                    j = k + 1
                    continue
                break
            items.append(("pp", text[i:k], line, "\n".join(comment)))
            comment = []
            line += text[i:k].count("\n")
            i = k
            continue
        if depth == 0:
            if c == "{":
                items.append(("code", "".join(code), code_line or line, "\n".join(comment)))
                code, code_line, comment = [], None, []
                depth, body, body_line = 1, [], line
            elif c == ";":
                items.append(("code", "".join(code), code_line or line, "\n".join(comment)))
                code, code_line, comment = [], None, []
            else:
                if not c.isspace() and code_line is None:
                    code_line = line
                code.append(c)
        else:
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    items.append(("body", "".join(body), body_line, ""))
                    i += 1
                    continue
            body.append(c)
        if c == "\n":
            line += 1
        i += 1
    return items


FUNCTION_HEAD = re.compile(
    r"^(?P<specifiers>(?:\w+[\s*]+)+?)(?P<name>[A-Za-z_]\w*)\s*\((?P<params>.*)\)\s*$", re.S)
NOT_FUNCTION = re.compile(r"^\s*(typedef|struct|enum|union)\b|=")


def c_functions(path: Path):
    """(name, line, is_static, preceding comment) of every function declared or defined."""
    for kind, code, line, comment in c_items(path.read_text()):
        head = " ".join(code.split())
        if kind != "code" or not head or NOT_FUNCTION.search(head):
            continue
        m = FUNCTION_HEAD.match(head)
        if m:
            yield m["name"], line, "static" in m["specifiers"].split(), comment


def kernel_c_files() -> list[Path]:
    return sorted(KERNEL_SOURCES.glob("*.c")) + sorted(KERNEL_HEADERS.glob("*.h"))


def test_every_c_function_has_a_basis():
    """DESIGN 5.1, Specification-driven; DESIGN 3.1: a comment directly above
    each function, in the header and in the source, carries the basis."""
    missing = [f"{rel(p)}:{line} {name}" for p in kernel_c_files()
               for name, line, _, comment in c_functions(p) if not MARKER.search(comment)]
    assert not missing, "C functions without [spec: …] or [design: …] in the comment above them:\n" \
        + "\n".join(missing)


def test_exported_c_functions_are_prefixed():
    """DESIGN 3.1: functions prefixed ``ssr_`` (a static helper is private to its source)."""
    wrong = [f"{rel(p)}:{line} {name}" for p in kernel_c_files()
             for name, line, static, _ in c_functions(p) if not static and not name.startswith("ssr_")]
    assert not wrong, "exported kernel functions without the ssr_ prefix:\n" + "\n".join(wrong)


def test_exported_c_functions_name_their_kernel():
    """DESIGN 3.1: an exported function is ``ssr_<kernel>_<what>``, <kernel> the
    file's name, since every kernel shares one C namespace."""
    wrong = [f"{rel(p)}:{line} {name}" for p in kernel_c_files()
             for name, line, static, _ in c_functions(p)
             if not static and not name.startswith(f"ssr_{p.stem}_")]
    assert not wrong, "exported kernel functions not named ssr_<kernel>_<what>:\n" + "\n".join(wrong)


STATUS_DEFINE = re.compile(r"^\s*#\s*define\s+SSR_(?P<kernel>[A-Z0-9]+)_(?P<name>\w+)\s+(?P<value>-?\d+)\b", re.M)


def test_status_codes_are_zero_for_success_and_distinct_positive_otherwise():
    """DESIGN 3.1: ``SSR_<KERNEL>_OK`` is 0 and every other condition of the
    kernel has its own positive code; the ``SSR_<KERNEL>_`` prefix is kept for
    status codes."""
    wrong = []
    for header in sorted(KERNEL_HEADERS.glob("*.h")):
        codes = {m["name"]: int(m["value"]) for m in STATUS_DEFINE.finditer(header.read_text())
                 if m["kernel"] == header.stem.upper()}
        if not codes:
            continue
        if codes.get("OK") != 0:
            wrong.append(f"{rel(header)}: SSR_{header.stem.upper()}_OK is not 0")
        failures = [v for k, v in codes.items() if k != "OK"]
        if any(v <= 0 for v in failures) or len(set(failures)) != len(failures):
            wrong.append(f"{rel(header)}: failure codes not distinct and positive: {codes}")
    assert not wrong, "\n".join(wrong)


def test_one_header_and_one_source_per_kernel():
    """DESIGN 3.1: ``kernels/include/ssrcore/<name>.h`` and ``kernels/src/<name>.c``,
    and a binding ``_kernels/<name>.pyx`` only for a kernel that exists."""
    sources = {p.stem for p in KERNEL_SOURCES.glob("*.c")}
    headers = {p.stem for p in KERNEL_HEADERS.glob("*.h")}
    bindings = {p.stem for p in BINDINGS.glob("*.pyx")}
    assert sources == headers, f"kernels without a header {sorted(sources - headers)}, " \
        f"headers without a source {sorted(headers - sources)}"
    assert bindings <= sources, f"bindings without a kernel: {sorted(bindings - sources)}"


def test_header_functions_are_defined_in_their_source():
    """DESIGN 3.1: a kernel's header declares what its own source defines."""
    wrong = []
    for header in sorted(KERNEL_HEADERS.glob("*.h")):
        source = KERNEL_SOURCES / f"{header.stem}.c"
        declared = {name for name, _, static, _ in c_functions(header) if not static}
        defined = {name for name, _, _, _ in c_functions(source)} if source.is_file() else set()
        wrong += [f"{rel(header)}: {name} not in {rel(source)}" for name in sorted(declared - defined)]
    assert not wrong, "\n".join(wrong)
    assert "ssr_scale_decode" in {name for name, *_ in c_functions(KERNEL_HEADERS / "scale.h")}


INCLUDE = re.compile(r'^\s*#\s*include\s*(?P<what><[^>]+>|"[^"]+")')
SYSTEM_HEADERS = {"<stdint.h>", "<math.h>"}


def test_kernels_include_only_allowed_headers():
    """DESIGN 3.1: ``<stdint.h>``, ``<math.h>`` and other kernel headers (``"ssrcore/<name>.h"``)."""
    wrong = []
    for path in kernel_c_files():
        for kind, text, line, _ in c_items(path.read_text()):
            m = INCLUDE.match(text) if kind == "pp" else None
            if not m:
                continue
            what = m["what"]
            if what in SYSTEM_HEADERS:
                continue
            local = re.fullmatch(r'"ssrcore/(\w+)\.h"', what)
            if local and (KERNEL_HEADERS / f"{local[1]}.h").is_file():
                continue
            wrong.append(f"{rel(path)}:{line} {what}")
    assert not wrong, "kernel includes outside <stdint.h>, <math.h> and kernel headers:\n" + "\n".join(wrong)


def strip_c(text: str) -> str:
    """The text without comments and string or character literals."""
    return re.sub(r'/\*.*?\*/|//[^\n]*|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'', " ", text, flags=re.S)


KERNEL_FORBIDDEN = {
    # DESIGN 3.1, Layers: explicit integer widths, double precision.
    r"\b(?:int|long|short|char|signed|unsigned)\b": "an integer type without explicit width (use <stdint.h>)",
    r"\bfloat\b": "single precision (DESIGN 4 and Q17: double unless a specification says otherwise)",
    # DESIGN 3.1: no allocation, no I/O, no errno reliance.
    r"\b(?:malloc|calloc|realloc|free|alloca)\s*\(": "allocation",
    r"\b(?:printf|fprintf|sprintf|snprintf|puts|fopen|fread|fwrite|fclose|FILE)\b": "I/O",
    r"\berrno\b": "errno",
}
KERNEL_FORBIDDEN_INT = next(iter(KERNEL_FORBIDDEN))
STATIC_VARIABLE = re.compile(r"\bstatic\b(?![^;{(]*\bconst\b)[^;{(]*[;=\[]")


def test_kernels_are_plain_arithmetic():
    """DESIGN 3.1, Layers and the kernel conventions: what a kernel may not contain."""
    wrong = []
    for path in kernel_c_files():
        code = strip_c(path.read_text())
        lines = code.split("\n")
        for pattern, what in KERNEL_FORBIDDEN.items():
            for k, text in enumerate(lines, 1):
                if not text.lstrip().startswith("#") and re.search(pattern, text):
                    wrong.append(f"{rel(path)}:{k} {what}: {text.strip()}")
        for m in STATIC_VARIABLE.finditer(code):
            k = code.count("\n", 0, m.start()) + 1
            wrong.append(f"{rel(path)}:{k} static mutable state: {lines[k - 1].strip()}")
    assert not wrong, "\n".join(wrong)


# ---- Python and Cython ------------------------------------------------------------------


def python_files() -> list[Path]:
    """The package and the build script; tests are not public functions."""
    return [p for p in source_files(".py") if not rel(p).startswith("tests/")]


def public_functions(tree: ast.Module):
    """(qualified name, node) of module functions and of methods of public classes,
    public meaning not starting with ``_``; dunder methods are not checked."""
    def walk(body, prefix):
        for node in body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith("_"):
                yield prefix + node.name, node
            elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
                yield from walk(node.body, f"{prefix}{node.name}.")
            elif isinstance(node, (ast.If, ast.Try)):  # definitions under a module-level condition
                yield from walk(node.body + node.orelse, prefix)
    yield from walk(tree.body, "")


def test_every_public_python_function_has_a_basis():
    """DESIGN 5.1, Specification-driven: the docstring carries the basis."""
    missing = []
    for path in python_files():
        tree = ast.parse(path.read_text(), filename=str(path))
        for name, node in public_functions(tree):
            if not MARKER.search(ast.get_docstring(node) or ""):
                missing.append(f"{rel(path)}:{node.lineno} {name}")
    assert not missing, "public functions without [spec: …] or [design: …] in the docstring:\n" \
        + "\n".join(missing)


PYX_DEF = re.compile(r"^(?P<indent>[ \t]*)(?:async\s+)?(?:def|cpdef)\s+(?:[\w\[\]\s\*]+?\s+)?"
                     r"(?P<name>[A-Za-z_]\w*)\s*\(", re.M)
DOCSTRING = re.compile(r'\s*[rRuU]?("""|\'\'\')(?P<doc>.*?)\1', re.S)


def pyx_functions(text: str):
    """(name, line, docstring) of each ``def`` and ``cpdef`` of a Cython binding."""
    for m in PYX_DEF.finditer(text):
        depth, j = 0, m.end() - 1  # past the parameter list, to the colon
        while j < len(text):
            depth += {"(": 1, ")": -1}.get(text[j], 0)
            if depth == 0 and text[j] == ":":
                break
            j += 1
        doc = DOCSTRING.match(text, j + 1)
        yield m["name"], text.count("\n", 0, m.start()) + 1, doc["doc"] if doc else ""


def test_every_cython_function_has_a_basis():
    """DESIGN 3.1: each binding function has a docstring repeating the citation."""
    missing = [f"{rel(p)}:{line} {name}" for p in sorted(BINDINGS.glob("*.pyx"))
               for name, line, doc in pyx_functions(p.read_text())
               if not name.startswith("_") and not MARKER.search(doc)]
    assert not missing, "binding functions without [spec: …] or [design: …] in the docstring:\n" \
        + "\n".join(missing)


def test_each_binding_declares_its_own_kernel_and_nothing_else():
    """DESIGN 3.1: ``_kernels/<name>.pyx`` declares ``ssrcore/<name>.h``, and a binding
    imports nothing from ssrcore (the kernels are the lowest layer)."""
    wrong = []
    for path in sorted(BINDINGS.glob("*.pyx")):
        text = path.read_text()
        headers = set(re.findall(r'cdef\s+extern\s+from\s+"([^"]+)"', text))
        if f"ssrcore/{path.stem}.h" not in headers:
            wrong.append(f"{rel(path)}: no cdef extern from \"ssrcore/{path.stem}.h\"")
        for m in re.finditer(r"^\s*(?:from\s+(\S+)\s+c?import|c?import\s+(\S+))", text, re.M):
            module = m[1] or m[2]
            if module.split(".")[0] in ("ssrcore", ".", "..") or module.startswith("."):
                wrong.append(f"{rel(path)}: imports {module}")
    assert not wrong, "\n".join(wrong)


# ---- layers -----------------------------------------------------------------------------

KERNELS, TOOLKIT, MODEL, INTERFACE = range(4)
LAYER_NAMES = {KERNELS: "kernels", TOOLKIT: "toolkit", MODEL: "model", INTERFACE: "interface"}


def layer(module: str) -> int | None:
    """DESIGN 3.1: the layer of an ssrcore module, None outside the package."""
    parts = module.split(".")
    if parts[0] != "ssrcore":
        return None
    if len(parts) == 1:
        return INTERFACE  # the package's public names
    return {"_kernels": KERNELS, "toolkit": TOOLKIT, "model": MODEL, "adapters": MODEL,
            "engine": INTERFACE, "record": INTERFACE, "vrs": INTERFACE}.get(parts[1], INTERFACE)


def module_name(path: Path) -> tuple[str, bool]:
    """Dotted name of a package file and whether it is a package ``__init__``."""
    parts = list(path.relative_to(ROOT / "src").with_suffix("").parts)
    package = parts[-1] == "__init__"
    if package:
        parts.pop()
    return ".".join(parts), package


def imports(path: Path):
    """(line, imported module) for every import of a module, nested ones included,
    relative imports resolved, and ``from a import b`` read as both a and a.b
    (only as a.b when a is the package root)."""
    name, package = module_name(path)
    base = name.split(".") if package else name.split(".")[:-1]
    for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield node.lineno, alias.name
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                anchor = base[:len(base) - node.level + 1]
                module = ".".join(anchor + ([node.module] if node.module else []))
            else:
                module = node.module or ""
            if module != "ssrcore":  # ``from ssrcore import x`` names x, read below
                yield node.lineno, module
            for alias in node.names:
                yield node.lineno, f"{module}.{alias.name}"


def test_no_module_imports_upward():
    """DESIGN 3.1: dependencies point down only (interface, model, toolkit, kernels);
    ``ssrcore._kernels`` is private to the toolkit; gnsswire is for the model and
    the interface, never for the toolkit or the kernels."""
    wrong = []
    for path in sorted(PACKAGE.rglob("*.py")):
        name, _ = module_name(path)
        own = layer(name)
        for line, module in imports(path):
            target = layer(module)
            where = f"{rel(path)}:{line} ({LAYER_NAMES[own]}) imports {module}"
            if target is not None and target > own:
                wrong.append(f"{where} ({LAYER_NAMES[target]})")
            elif target == KERNELS and own not in (KERNELS, TOOLKIT):
                wrong.append(f"{where}: the kernels are reached through the toolkit")
            elif module.split(".")[0] == "gnsswire" and own in (KERNELS, TOOLKIT):
                wrong.append(f"{where}: the toolkit takes numbers, not messages")
    assert not wrong, "\n".join(sorted(set(wrong)))


# Modules and calls by which code reaches files, sockets, processes or clocks.
IO_MODULES = {"socket", "select", "selectors", "ssl", "asyncio", "subprocess", "os", "io", "pathlib",
              "shutil", "tempfile", "time", "random", "secrets", "urllib", "http", "ftplib", "logging"}
IO_CALLS = {"open", "input", "print", "breakpoint"}
CLOCK_ATTRIBUTES = {"now", "today", "utcnow"}


def test_core_does_no_io_and_reads_no_clock():
    """DESIGN 2.1, I/O at the edge: the package opens no file or socket and reads no
    clock, so that every run can be replayed from recorded inputs."""
    wrong = []
    for path in sorted(PACKAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for line, module in imports(path):
            if module.split(".")[0] in IO_MODULES:
                wrong.append(f"{rel(path)}:{line} imports {module}")
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                f = node.func
                if isinstance(f, ast.Name) and f.id in IO_CALLS:
                    wrong.append(f"{rel(path)}:{node.lineno} calls {f.id}()")
                elif isinstance(f, ast.Attribute) and f.attr in CLOCK_ATTRIBUTES:
                    wrong.append(f"{rel(path)}:{node.lineno} calls .{f.attr}()")
    assert not wrong, "\n".join(wrong)


# ---- REVIEW markers ---------------------------------------------------------------------

# Built so that this file does not match itself.
REVIEW = re.compile(r"(?:#|//|/\*|\*)\s*" + "REV" + r"IEW\s*:")


def review_markers(root: Path = ROOT) -> list[str]:
    found = []
    for path in source_files(".py", ".pyx", ".pxd", ".c", ".h", root=root):
        for k, text in enumerate(path.read_text().split("\n"), 1):
            if REVIEW.search(text):
                found.append(f"{path.relative_to(root)}:{k}")
    return found


def check_review_markers(root: Path = ROOT) -> None:
    """DESIGN 5.1, Specification-driven: a REVIEW marker is neither kind, so none remains."""
    found = review_markers(root)
    if not found:
        return
    message = f"{len(found)} REVIEW marker(s) without a basis: " + ", ".join(found)
    if REVIEW_MARKERS_FAIL:
        pytest.fail(message)
    pytest.skip(message + " (reported, not failed: REVIEW_MARKERS_FAIL is off)")


def test_review_markers():
    """DESIGN 5.1, Specification-driven: no REVIEW marker in the repository, and one fails the gate."""
    assert REVIEW_MARKERS_FAIL
    check_review_markers()


def test_a_review_marker_fails_the_gate(tmp_path):
    """The check fails, not skips, on a marker in any source file."""
    (tmp_path / "src").mkdir()
    (tmp_path / "src" / "m.py").write_text("x = 1  # " + "REV" + "IEW: not taken from a specification\n")
    with pytest.raises(pytest.fail.Exception, match="1 REVIEW marker"):
        check_review_markers(tmp_path)


# ---- the checks themselves --------------------------------------------------------------


SAMPLE_C = '''/* File comment. */
#include <stdint.h>
#include "ssrcore/nothing.h"

/* [spec: DOC 1.2] Documented. */
double ssr_good(double x, int32_t n);

/* No basis here. */
static double helper(double x)
{
    if (x > 0) { return x; }
    return -x;
}

static const double TABLE[2] = {1.0, 2.0};
static double cache;
typedef struct { double a; } thing;
/* [design: reason] */ int32_t ssr_two(const double v[3],
                                      double out[3]);
'''


def test_c_scanner_on_a_sample(tmp_path):
    """The C checks find what they should: functions with their comments, includes,
    static state; and not struct typedefs or constant tables."""
    path = tmp_path / "sample.c"
    path.write_text(SAMPLE_C)
    functions = {name: (static, bool(MARKER.search(comment)))
                 for name, _, static, comment in c_functions(path)}
    assert functions == {"ssr_good": (False, True), "helper": (True, False), "ssr_two": (False, True)}
    includes = [INCLUDE.match(t)["what"] for k, t, _, _ in c_items(SAMPLE_C) if k == "pp" and INCLUDE.match(t)]
    assert includes == ["<stdint.h>", '"ssrcore/nothing.h"']
    statics = [m.group(0) for m in STATIC_VARIABLE.finditer(strip_c(SAMPLE_C))]
    assert statics == ["static double cache;"]
    assert re.search(KERNEL_FORBIDDEN_INT, "    int i = 0;") and not re.search(KERNEL_FORBIDDEN_INT, "int32_t i;")



def test_python_checks_on_a_sample():
    """Public functions and methods are found, private ones and dunders are not;
    a marker needs a non-empty basis."""
    tree = ast.parse('def f():\n    "[spec: X 1]"\n'
                     'def _g(): pass\n'
                     'class C:\n    def m(self): pass\n    def __add__(self, o): pass\n'
                     'class _D:\n    def m(self): pass\n')
    assert [n for n, _ in public_functions(tree)] == ["f", "C.m"]
    assert MARKER.search("[design: a reason]") and not MARKER.search("[spec: ]")
    text = 'def plain(double x):\n    """[spec: SPARTN 6.3] x."""\n    return x\n' \
           'cpdef double typed(double x) nogil:\n    return x\n'
    assert [(n, bool(MARKER.search(d))) for n, _, d in pyx_functions(text)] == [("plain", True), ("typed", False)]


def test_layer_assignment():
    """DESIGN 3.1: the four layers and where each package belongs."""
    assert [layer(m) for m in ("ssrcore._kernels.scale", "ssrcore.toolkit.geodesy", "ssrcore.model.types",
                               "ssrcore.adapters.rtcm3", "ssrcore.engine", "ssrcore.record",
                               "ssrcore.vrs.station", "ssrcore", "gnsswire.spartn")] \
        == [KERNELS, TOOLKIT, MODEL, MODEL, INTERFACE, INTERFACE, INTERFACE, INTERFACE, None]


def test_review_pattern():
    """The marker is found in Python and C comments, not in prose that names it."""
    word = "REV" + "IEW"
    assert REVIEW.search(f"# {word}: not taken from a specification")
    assert REVIEW.search(f"/* {word}: x */") and REVIEW.search(f"    // {word}: x")
    assert not REVIEW.search(f"a `{word}` marker")
