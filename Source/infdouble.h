#pragma once

static inline double infplus(double a, double b)
{
    return (a < 0 || b < 0)? -1 : a + b;
}

static inline double infmin(double a, double b)
{
    if (a < 0) return b;
    if (b < 0) return a;
    return (a < b)? a : b;
}

static inline int inflt(double a, double b)
{
    if (a < 0) return 0; /* if a == inf, then surely a >= b */
    if (b < 0) return 1; /* if a < inf and b == inf, then a < b */
    return (!!(a < b));
}
