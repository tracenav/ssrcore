/* Driver of the counting build (count.hpp): the operations of each part of
 * one query, by kind and by translation unit, and the values, which must be
 * those of the C build bit for bit.
 *
 * Output, one per line: "count <part> <unit> <kind> <n>" for every non-zero
 * counter, then "value <name> <17 digits>" for the results main.c prints
 * under the same names. [design: tools/portability.py reads it] */
#include "query.h"

static const char *const KIND_NAMES[SSRQ_KINDS] = {
    "add", "mul", "div", "cmp", "neg", "conv",
    "sqrt", "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
    "exp", "log", "pow", "hypot", "fabs", "floor", "ceil", "fmod",
    "round", "isfinite"};

static ssrq_counts registry[64];
static int registered = 0;

/* [design: one block of counters per translation unit, named for it] */
ssrq_counts *ssrq_register(const char *unit)
{
    for (int k = 0; k < registered; k++) {
        if (strcmp(registry[k].unit, unit) == 0) {
            return &registry[k];
        }
    }
    if (registered == 64) {
        abort();
    }
    registry[registered].unit = unit;
    return &registry[registered++];
}

static void reset()
{
    for (int k = 0; k < registered; k++) {
        memset(registry[k].n, 0, sizeof registry[k].n);
    }
}

static void report(const char *part)
{
    for (int k = 0; k < registered; k++) {
        for (int j = 0; j < SSRQ_KINDS; j++) {
            if (registry[k].n[j]) {
                printf("count %s %s %s %lu\n", part, registry[k].unit, KIND_NAMES[j], registry[k].n[j]);
            }
        }
    }
}

static void value(const char *name, ssrq_counted x)
{
    printf("value %s %.17g\n", name, x.v);
}

static void vector(const char *name, const ssrq_counted *x, int n)
{
    char buffer[64];
    for (int k = 0; k < n; k++) {
        snprintf(buffer, sizeof buffer, "%s[%d]", name, k);
        value(buffer, x[k]);
    }
}

int main()
{
    static ssrq_inputs in;
    static ssrq_receiver_state r;
    static ssrq_choice_state c;
    static ssrq_epoch_state e;
    static ssrq_query_state q;
    ssrq_scene(&in);
    reset();
    if (ssrq_receiver(&in, &r) != SSRQ_OK) return 1;
    report("receiver");
    reset();
    if (ssrq_choice(&in, &c) != SSRQ_OK) return 1;
    report("choice");
    reset();
    if (ssrq_epoch(&in, &r, &e) != SSRQ_OK) return 1;
    report("epoch");
    reset();
    if (ssrq_query(&in, &r, &c, &e, &q) != SSRQ_OK) return 1;
    report("query");
    vector("epoch.sun", e.sun, 3);
    vector("epoch.moon", e.moon, 3);
    vector("epoch.tide_displacement", e.tide_displacement, 3);
    value("epoch.hydrostatic_zenith", e.hydrostatic_zenith);
    vector("query.position", q.position, 3);
    value("query.range", q.range);
    value("query.clock", q.clock);
    vector("query.unit", q.unit, 3);
    value("query.shapiro", q.shapiro);
    value("query.relativity_clock", q.relativity_clock);
    value("query.solid_tide", q.solid_tide);
    value("query.windup", q.windup);
    value("query.tropo_slant", q.tropo_slant);
    value("query.code_delay[0]", q.code_delay[0]);
    value("query.code_delay[1]", q.code_delay[1]);
    return 0;
}
