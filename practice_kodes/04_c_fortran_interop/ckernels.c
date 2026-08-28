/* ==========================================================================
 * EXERCISE 04 — the C side of the Fortran/C boundary
 * ==========================================================================
 *
 * This file plays the role of YAC: a C library that an application written
 * in Fortran (ICON) calls into. Everything awkward about that boundary is
 * represented here — array layout, communicator handles, strings, structs.
 *
 * Fill in the TODOs marked C1..C5. See interop.f90 for the Fortran side.
 * ========================================================================== */

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* --------------------------------------------------------------------------
 * C1: a plain kernel called from Fortran.
 *
 * Fortran passes arrays by reference, so `x` and `y` arrive as pointers with
 * no copy — provided the Fortran side declares the interface with bind(C)
 * and passes a contiguous array. Compute y = a*x + y.
 * ------------------------------------------------------------------------ */
void c_axpy(int n, double a, const double *x, double *y)
{
    /* TODO C1: y[i] = a * x[i] + y[i] */
    (void)n; (void)a; (void)x; (void)y;
}

/* --------------------------------------------------------------------------
 * C2: THE classic interop bug — array index order.
 *
 * Fortran stores m(nrow, ncol) COLUMN-major: m(1,1), m(2,1), m(3,1), ...
 * C stores m[nrow][ncol] ROW-major:          m[0][0], m[0][1], m[0][2], ...
 *
 * The same flat buffer therefore means different things on each side. This
 * function receives the raw pointer from Fortran and must sum each Fortran
 * COLUMN, indexing the flat array the way Fortran laid it out:
 *
 *     Fortran m(i, j)  ==  flat[(j-1)*nrow + (i-1)]
 *
 * Write colsum[j] = sum over i of m(i,j). Getting this wrong does not
 * crash — it silently returns transposed garbage, which is why this bug
 * survives so long in real coupled models.
 * ------------------------------------------------------------------------ */
void c_column_sums(int nrow, int ncol, const double *m, double *colsum)
{
    /* TODO C2: loop j = 0..ncol-1, i = 0..nrow-1, colsum[j] += m[j*nrow + i] */
    (void)nrow; (void)ncol; (void)m; (void)colsum;
}

/* --------------------------------------------------------------------------
 * C3: receiving an MPI communicator from Fortran.
 *
 * Fortran and C use different handle representations. Fortran handles are
 * integers; C handles are opaque pointers/structs. The standard conversion
 * functions are MPI_Comm_f2c() and MPI_Comm_c2f().
 *
 * With `use mpi_f08` the Fortran integer handle lives in comm%MPI_VAL.
 * This is exactly how YAC's C API accepts a communicator from ICON.
 *
 * Do an MPI_Allreduce on the converted communicator and return the sum.
 * ------------------------------------------------------------------------ */
double c_sum_over_comm(int fortran_comm, double local_value)
{
    MPI_Comm comm;
    double total = 0.0;

    /* TODO C3a: comm = MPI_Comm_f2c(fortran_comm); */
    /* TODO C3b: MPI_Allreduce(&local_value, &total, 1, MPI_DOUBLE,
     *                         MPI_SUM, comm);                            */

    (void)comm; (void)fortran_comm;
    total = local_value;   /* placeholder: no reduction */
    return total;
}

/* --------------------------------------------------------------------------
 * C4: strings.
 *
 * C strings are NUL-terminated; Fortran strings are fixed-length and blank-
 * padded with the length passed separately (invisibly, in the old calling
 * convention). Under bind(C) the Fortran side must append c_null_char
 * itself and declare the dummy as character(kind=c_char), dimension(*).
 *
 * YAC's API is full of these: yac_cdef_field("sea_surface_temperature", ...).
 * Return the length of the field name and stash it so Fortran can verify
 * the string survived the crossing intact.
 * ------------------------------------------------------------------------ */
static char last_name[256] = {0};

int c_register_field(const char *name)
{
    /* TODO C4: copy `name` into last_name (bounded!), return (int)strlen(name).
     *          Use snprintf, not strcpy — the buffer is fixed size.         */
    (void)name;
    return -1;
}

/* Returns 1 if the last registered name matches `expect`, else 0.
 * Used by the Fortran side to prove the string crossed correctly. */
int c_check_last_field(const char *expect)
{
    return strcmp(last_name, expect) == 0 ? 1 : 0;
}

/* --------------------------------------------------------------------------
 * C5: interoperable structs.
 *
 * A bind(C) derived type in Fortran and a struct in C must agree on member
 * order, types, AND padding. This one is deliberately laid out so that
 * naive ordering would introduce padding on most ABIs — check that
 * sizeof() matches what Fortran's c_sizeof() reports.
 * ------------------------------------------------------------------------ */
typedef struct {
    int    grid_id;      /* 4 bytes                                  */
    int    n_cells;      /* 4 bytes  -- pairs with grid_id, no pad    */
    double dx;           /* 8 bytes, 8-aligned                        */
    double dy;           /* 8 bytes                                   */
} grid_desc_t;

/* Fortran passes the struct by reference; fill in a derived quantity. */
double c_grid_area(const grid_desc_t *g)
{
    /* TODO C5: return n_cells * dx * dy */
    (void)g;
    return 0.0;
}

/* Lets Fortran compare sizeof() across the boundary. */
int c_sizeof_grid_desc(void) { return (int)sizeof(grid_desc_t); }

/* --------------------------------------------------------------------------
 * C6 (read-only): C calling BACK into Fortran.
 *
 * This is the callback direction — YAC does this for user-supplied
 * interpolation. The Fortran routine is declared bind(C, name="f_scale")
 * and we call it here through a normal prototype.
 * ------------------------------------------------------------------------ */
extern void f_scale(int n, double factor, double *v);

void c_calls_fortran(int n, double factor, double *v)
{
    printf("    [C] calling back into Fortran f_scale(n=%d, factor=%.1f)\n",
           n, factor);
    f_scale(n, factor, v);
}
