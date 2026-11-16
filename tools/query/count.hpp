/* Counting build of the kernels and of query.c (C++ only; not part of the
 * kernels or of the firmware artefact).
 *
 * Every source is compiled as C++ with this header forced in first
 * (-include) and SSRQ_UNIT set to its name. After the system headers, the
 * token double is defined to ssrq_counted, a class that holds one double and
 * does each operation exactly as the C build does (same operation, same
 * order, IEEE 754 binary64, no contraction), counting it by kind in the
 * counters of its translation unit: so the number of double-precision
 * operations a query makes, by kind and by kernel, is what a target without
 * a double-precision FPU executes as software routines (tools/portability.py,
 * docs/DESIGN.md Q17).
 *
 * Counted: add and subtract, multiply, divide, compare, negate, conversion
 * of an integer to double and of a double to an integer, and each math
 * library call by name. A double literal becomes an ssrq_counted without a
 * count; an operation between two literals is done by the compiler. An
 * integer literal mixed with a double, which a C compiler converts at
 * compile time, is counted as a conversion, so that count is an upper bound.
 * [design: an exact count of what the code asks for, independent of any
 * target's compiler] */
#ifndef SSRCORE_QUERY_COUNT_HPP
#define SSRCORE_QUERY_COUNT_HPP

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <cmath>

#ifndef SSRQ_UNIT
#define SSRQ_UNIT "unknown"
#endif

enum ssrq_kind {
    SSRQ_ADD, SSRQ_MUL, SSRQ_DIV, SSRQ_CMP, SSRQ_NEG, SSRQ_CONV,
    SSRQ_SQRT, SSRQ_SIN, SSRQ_COS, SSRQ_TAN, SSRQ_ASIN, SSRQ_ACOS, SSRQ_ATAN, SSRQ_ATAN2,
    SSRQ_EXP, SSRQ_LOG, SSRQ_POW, SSRQ_HYPOT, SSRQ_FABS, SSRQ_FLOOR, SSRQ_CEIL, SSRQ_FMOD,
    SSRQ_ROUND, SSRQ_ISFINITE, SSRQ_KINDS
};

/* One translation unit's counters; the registry is defined in count.cpp. */
struct ssrq_counts {
    const char *unit;
    unsigned long n[SSRQ_KINDS];
};
ssrq_counts *ssrq_register(const char *unit);

/* The counters of this translation unit (internal linkage, one per unit). */
static ssrq_counts *const ssrq_here = ssrq_register(SSRQ_UNIT);

typedef double ssrq_f64;

struct ssrq_counted {
    ssrq_f64 v;
    constexpr ssrq_counted() : v(0.0) {}
    constexpr ssrq_counted(ssrq_f64 x) : v(x) {}
    ssrq_counted(int32_t x) : v((ssrq_f64)x) { ssrq_here->n[SSRQ_CONV]++; }
    ssrq_counted(int64_t x) : v((ssrq_f64)x) { ssrq_here->n[SSRQ_CONV]++; }
    ssrq_counted(uint32_t x) : v((ssrq_f64)x) { ssrq_here->n[SSRQ_CONV]++; }
    explicit operator int32_t() const { ssrq_here->n[SSRQ_CONV]++; return (int32_t)v; }
    explicit operator int64_t() const { ssrq_here->n[SSRQ_CONV]++; return (int64_t)v; }
    explicit operator uint32_t() const { ssrq_here->n[SSRQ_CONV]++; return (uint32_t)v; }
    ssrq_counted &operator+=(ssrq_counted b) { ssrq_here->n[SSRQ_ADD]++; v = v + b.v; return *this; }
    ssrq_counted &operator-=(ssrq_counted b) { ssrq_here->n[SSRQ_ADD]++; v = v - b.v; return *this; }
    ssrq_counted &operator*=(ssrq_counted b) { ssrq_here->n[SSRQ_MUL]++; v = v * b.v; return *this; }
    ssrq_counted &operator/=(ssrq_counted b) { ssrq_here->n[SSRQ_DIV]++; v = v / b.v; return *this; }
};

static inline ssrq_counted operator+(ssrq_counted a, ssrq_counted b) { ssrq_here->n[SSRQ_ADD]++; return a.v + b.v; }
static inline ssrq_counted operator-(ssrq_counted a, ssrq_counted b) { ssrq_here->n[SSRQ_ADD]++; return a.v - b.v; }
static inline ssrq_counted operator*(ssrq_counted a, ssrq_counted b) { ssrq_here->n[SSRQ_MUL]++; return a.v * b.v; }
static inline ssrq_counted operator/(ssrq_counted a, ssrq_counted b) { ssrq_here->n[SSRQ_DIV]++; return a.v / b.v; }
static inline ssrq_counted operator-(ssrq_counted a) { ssrq_here->n[SSRQ_NEG]++; return -a.v; }
static inline ssrq_counted operator+(ssrq_counted a) { return a; }
static inline bool operator<(ssrq_counted a, ssrq_counted b) { ssrq_here->n[SSRQ_CMP]++; return a.v < b.v; }
static inline bool operator>(ssrq_counted a, ssrq_counted b) { ssrq_here->n[SSRQ_CMP]++; return a.v > b.v; }
static inline bool operator<=(ssrq_counted a, ssrq_counted b) { ssrq_here->n[SSRQ_CMP]++; return a.v <= b.v; }
static inline bool operator>=(ssrq_counted a, ssrq_counted b) { ssrq_here->n[SSRQ_CMP]++; return a.v >= b.v; }
static inline bool operator==(ssrq_counted a, ssrq_counted b) { ssrq_here->n[SSRQ_CMP]++; return a.v == b.v; }
static inline bool operator!=(ssrq_counted a, ssrq_counted b) { ssrq_here->n[SSRQ_CMP]++; return a.v != b.v; }
static inline bool operator!(ssrq_counted a) { ssrq_here->n[SSRQ_CMP]++; return !a.v; }

#define SSRQ_MATH1(name, kind) \
    static inline ssrq_counted name(ssrq_counted a) { ssrq_here->n[kind]++; return ::name(a.v); }
#define SSRQ_MATH2(name, kind) \
    static inline ssrq_counted name(ssrq_counted a, ssrq_counted b) { ssrq_here->n[kind]++; return ::name(a.v, b.v); }
SSRQ_MATH1(sqrt, SSRQ_SQRT)
SSRQ_MATH1(sin, SSRQ_SIN)
SSRQ_MATH1(cos, SSRQ_COS)
SSRQ_MATH1(tan, SSRQ_TAN)
SSRQ_MATH1(asin, SSRQ_ASIN)
SSRQ_MATH1(acos, SSRQ_ACOS)
SSRQ_MATH1(atan, SSRQ_ATAN)
SSRQ_MATH2(atan2, SSRQ_ATAN2)
SSRQ_MATH1(exp, SSRQ_EXP)
SSRQ_MATH1(log, SSRQ_LOG)
SSRQ_MATH2(pow, SSRQ_POW)
SSRQ_MATH2(hypot, SSRQ_HYPOT)
SSRQ_MATH1(fabs, SSRQ_FABS)
SSRQ_MATH1(floor, SSRQ_FLOOR)
SSRQ_MATH1(ceil, SSRQ_CEIL)
SSRQ_MATH2(fmod, SSRQ_FMOD)
SSRQ_MATH1(round, SSRQ_ROUND)
static inline bool isfinite(ssrq_counted a) { ssrq_here->n[SSRQ_ISFINITE]++; return std::isfinite(a.v); }

#define double ssrq_counted

#endif
