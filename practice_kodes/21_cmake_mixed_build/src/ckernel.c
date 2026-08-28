/* A C kernel called from Fortran via bind(C) -- the ICON/YAC pattern from
 * Exercise 04, here to make the build genuinely mixed-language. */
#include <math.h>

void c_scale_and_sum(int n, double factor, const double *x, double *out)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += factor * x[i];
    *out = s;
}
