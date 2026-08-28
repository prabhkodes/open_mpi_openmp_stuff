/* ==========================================================================
 * EXERCISE 11 — Spatial search and interpolation on the sphere
 * ==========================================================================
 *
 * GOAL
 *   Given a source grid of points on a sphere and a set of target points,
 *   find each target's nearest source point (and its k nearest), then
 *   interpolate. Do it three ways and measure the cost of each:
 *     1. brute force        O(n_tgt * n_src)   -- the correctness reference
 *     2. lat-lon bucketing  O(n_tgt * bucket)  -- simple, and usually enough
 *     3. k-d tree           O(n_tgt * log n)   -- what libraries actually use
 *
 * WHY (DKRZ)
 *   Exercise 10 gave you the weights once the overlap was known. This is the
 *   step BEFORE that: finding which source cells are anywhere near a target
 *   cell at all. For two unstructured grids with millions of cells each,
 *   the naive answer is 10^12 distance evaluations, which is why coupler
 *   SETUP — not the timestep — is often the thing that blows up when a
 *   scientist doubles the resolution.
 *
 *   "Experience in algorithms and data structures" is explicitly on the
 *   posting's wanted list. This is that bullet point, in the exact context
 *   the job cares about. It is also the exercise most likely to turn into a
 *   whiteboard question.
 *
 *   Written in C on purpose — YAC is C, and the posting asks for C/C++.
 *
 * TASKS
 *   TODO 1  great_circle_distance  — and why not Euclidean
 *   TODO 2  brute_force_nearest    — the reference
 *   TODO 3  bucket_build / bucket_nearest — lat-lon binning
 *   TODO 4  kdtree_build / kdtree_nearest — 3D k-d tree on the unit sphere
 *   TODO 5  idw_interpolate        — inverse-distance weighting from k
 *   TODO 6  answer the questions at the bottom
 *
 * ACCEPTANCE
 *   - bucket and k-d tree return EXACTLY the same nearest neighbour as
 *     brute force for every target point (this is not approximate --
 *     an off-by-one in a bucket boundary shows up here immediately)
 *   - a timing table showing the asymptotic difference
 *   - you can explain why the bucket method degrades near the poles and
 *     what you would do about it
 *
 * HINTS
 *   - Work in 3D Cartesian on the unit sphere. Then "nearest in great-circle
 *     distance" is exactly "nearest in Euclidean chord distance" (the map
 *     between them is monotone), so a plain 3D k-d tree is correct with no
 *     spherical geometry inside it at all. That trick is the whole reason
 *     this is tractable.
 *   - Do NOT use acos() in the inner loop. Compare squared chord distances
 *     and convert only at the end.
 *   - Bucket sizing: aim for ~1-4 points per bucket. Near the poles the
 *     lat-lon cells get narrow, so buckets there hold too few points and
 *     you must search more neighbouring buckets -- that is TODO 6b.
 *   - k-d tree build: median split on the widest dimension. nth_element is
 *     the textbook approach; a simple quickselect is fine here.
 * ========================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double x, y, z; } vec3;

typedef struct {
    int     n;
    double *lon, *lat;     /* radians */
    vec3   *xyz;           /* unit vectors -- the representation that matters */
} pointset;

static double wtime(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static vec3 to_xyz(double lon, double lat)
{
    vec3 v;
    v.x = cos(lat) * cos(lon);
    v.y = cos(lat) * sin(lon);
    v.z = sin(lat);
    return v;
}

/* --------------------------------------------------------------------------
 * TODO 1: great-circle distance between two unit vectors.
 *
 * The obvious formula is  d = acos(dot(a,b)).  It is also numerically
 * terrible for nearby points: as dot -> 1, acos loses almost all its
 * significant digits, and "nearby points" is precisely the case you care
 * about in a nearest-neighbour search.
 *
 * Use the atan2 form instead:
 *     d = atan2(|a x b|, a.b)
 * which is well conditioned over the whole range.
 *
 * Then note that for SEARCH you do not need this at all: chord distance
 * |a-b|^2 is a monotone function of great-circle distance, so ranking by
 * chord^2 gives an identical ordering with no transcendental functions.
 * Use chord^2 in the hot loops, and this function only for reporting.
 * ------------------------------------------------------------------------ */
double great_circle_distance(vec3 a, vec3 b)
{
    /* TODO 1: implement the atan2 form. */
    (void)a; (void)b;
    return 0.0;
}

static inline double chord2(vec3 a, vec3 b)
{
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx*dx + dy*dy + dz*dz;
}

/* --------------------------------------------------------------------------
 * TODO 2: brute force. The reference implementation.
 *
 * For each target point, scan every source point and keep the closest.
 * Write it simply and correctly -- everything else is validated against it,
 * so a bug here poisons the whole exercise.
 *
 * Fill nearest[j] with the index of the closest source point to target j.
 * ------------------------------------------------------------------------ */
void brute_force_nearest(const pointset *src, const pointset *tgt, int *nearest)
{
    /* TODO 2 */
    for (int j = 0; j < tgt->n; j++) nearest[j] = -1;
    (void)src;
}

/* ==========================================================================
 * TODO 3: lat-lon bucket grid.
 *
 * Bin every source point into a coarse (nlon x nlat) lat-lon grid. To query,
 * look in the target's own bucket first, then expand outward ring by ring
 * until the closest point found is provably closer than the nearest possible
 * point in the next ring.
 *
 * That termination test is the part people get wrong. Stopping as soon as
 * you find ANY point gives a wrong answer near bucket boundaries, and it is
 * wrong rarely enough that a weak test will not catch it -- which is exactly
 * why the acceptance criterion here is EXACT agreement with brute force.
 * ========================================================================== */
typedef struct {
    int nlon, nlat;
    int *bucket_start;   /* size nlon*nlat + 1, CSR-style */
    int *bucket_items;   /* size src->n, point indices sorted by bucket */
    const pointset *src;
} bucketgrid;

void bucket_build(bucketgrid *bg, const pointset *src, int nlon, int nlat)
{
    bg->src = src; bg->nlon = nlon; bg->nlat = nlat;
    bg->bucket_start = calloc((size_t)nlon*nlat + 1, sizeof(int));
    bg->bucket_items = malloc((size_t)src->n * sizeof(int));

    /* TODO 3a: counting sort into buckets.
     *   pass 1: count points per bucket
     *   prefix sum -> bucket_start
     *   pass 2: place each point index into bucket_items
     * This CSR layout is worth the extra pass: it is contiguous, has no
     * per-bucket allocation, and the query loop stays cache-friendly. */
    for (int i = 0; i < src->n; i++) bg->bucket_items[i] = i;
}

void bucket_nearest(const bucketgrid *bg, const pointset *tgt, int *nearest)
{
    /* TODO 3b: ring search with a correct termination test. */
    for (int j = 0; j < tgt->n; j++) nearest[j] = -1;
    (void)bg;
}

void bucket_free(bucketgrid *bg)
{
    free(bg->bucket_start); free(bg->bucket_items);
}

/* ==========================================================================
 * TODO 4: 3D k-d tree.
 *
 * Because the points live on the unit sphere and chord distance ranks the
 * same as great-circle distance, a plain 3D k-d tree over (x,y,z) is exactly
 * correct. No spherical geometry inside the tree at all.
 *
 * Build: recursively split the point set at the median of its widest
 * dimension. Store the tree in a flat array -- no per-node malloc.
 *
 * Query: descend to the leaf containing the query point, then unwind,
 * checking at each level whether the splitting plane is closer than the best
 * distance so far. If it is, the other subtree could contain something
 * better and must be searched too.
 *
 * That backtracking step is what makes it exact. Skip it and you have an
 * approximate nearest-neighbour search -- which is sometimes what you want,
 * but you must know which one you built.
 * ========================================================================== */
typedef struct {
    int   *idx;        /* permutation of point indices                 */
    int   *left, *right, *split_dim;
    double *split_val;
    int    nnodes, root;
    const pointset *src;
} kdtree;

void kdtree_build(kdtree *t, const pointset *src)
{
    t->src = src;
    t->idx = malloc((size_t)src->n * sizeof(int));
    for (int i = 0; i < src->n; i++) t->idx[i] = i;
    t->nnodes = 0; t->root = -1;
    t->left = t->right = t->split_dim = NULL;
    t->split_val = NULL;
    /* TODO 4a: build the tree. */
}

void kdtree_nearest(const kdtree *t, const pointset *tgt, int *nearest)
{
    /* TODO 4b: query with backtracking. */
    for (int j = 0; j < tgt->n; j++) nearest[j] = -1;
    (void)t;
}

void kdtree_free(kdtree *t)
{
    free(t->idx); free(t->left); free(t->right);
    free(t->split_dim); free(t->split_val);
}

/* --------------------------------------------------------------------------
 * TODO 5: inverse-distance weighting from the k nearest sources.
 *
 *     F_j = sum_k ( f_k / d_k^p ) / sum_k ( 1 / d_k^p )      p = 2 typically
 *
 * Two things to get right:
 *   - d = 0 exactly (target coincides with a source point) must not produce
 *     a NaN. Return that source value directly.
 *   - The weights sum to 1 by construction, so IDW is CONSISTENT (a constant
 *     field survives). It is NOT conservative -- nothing here knows about
 *     cell areas. Compare with Exercise 10 and be able to say which coupled
 *     fields you would and would not use this for.
 * ------------------------------------------------------------------------ */
void idw_interpolate(const pointset *src, const pointset *tgt,
                     const double *f_src, double *f_tgt, int k)
{
    /* TODO 5: for each target, find k nearest, weight by 1/d^2. */
    for (int j = 0; j < tgt->n; j++) f_tgt[j] = 0.0;
    (void)src; (void)f_src; (void)k;
}

/* ========================================================================== */

static void make_random_sphere(pointset *p, int n, unsigned seed)
{
    p->n = n;
    p->lon = malloc((size_t)n * sizeof(double));
    p->lat = malloc((size_t)n * sizeof(double));
    p->xyz = malloc((size_t)n * sizeof(vec3));
    srand(seed);
    for (int i = 0; i < n; i++) {
        double u = (double)rand() / RAND_MAX;
        double v = (double)rand() / RAND_MAX;
        /* Uniform on the sphere: lat = asin(2v-1), NOT lat = pi*(v-0.5).
           The naive version clusters points at the poles -- a classic bug
           in test-grid generators, and worth recognising on sight. */
        p->lon[i] = 2.0 * M_PI * u - M_PI;
        p->lat[i] = asin(2.0 * v - 1.0);
        p->xyz[i] = to_xyz(p->lon[i], p->lat[i]);
    }
}

static void free_pointset(pointset *p)
{
    free(p->lon); free(p->lat); free(p->xyz);
}

static int compare_results(const int *a, const int *b, int n, const char *label)
{
    int bad = 0, first = -1;
    for (int j = 0; j < n; j++) {
        if (a[j] != b[j]) { bad++; if (first < 0) first = j; }
    }
    if (bad == 0) {
        printf("     %-16s MATCHES brute force  (all %d targets)\n", label, n);
        return 1;
    }
    printf("     %-16s DIFFERS on %d / %d targets (first at j=%d)\n",
           label, bad, n, first);
    return 0;
}

int main(int argc, char **argv)
{
    int n_src = 20000, n_tgt = 20000, k = 4;
    if (argc > 1) n_src = atoi(argv[1]);
    if (argc > 2) n_tgt = atoi(argv[2]);

    pointset src, tgt;
    make_random_sphere(&src, n_src, 12345);
    make_random_sphere(&tgt, n_tgt, 67890);

    printf("=== Exercise 11: spatial search on the sphere ===\n");
    printf("source points %d   target points %d\n", n_src, n_tgt);
    printf("brute force would be %.2e distance evaluations\n",
           (double)n_src * (double)n_tgt);
    printf("\n");

    int *nn_brute  = malloc((size_t)n_tgt * sizeof(int));
    int *nn_bucket = malloc((size_t)n_tgt * sizeof(int));
    int *nn_kd     = malloc((size_t)n_tgt * sizeof(int));

    double t0, t_brute, t_bucket_build, t_bucket, t_kd_build, t_kd;

    t0 = wtime();
    brute_force_nearest(&src, &tgt, nn_brute);
    t_brute = wtime() - t0;
    printf("  brute force      search %8.4f s\n", t_brute);

    bucketgrid bg;
    int nb = (int)(sqrt((double)n_src) / 2.0); if (nb < 4) nb = 4;
    t0 = wtime();
    bucket_build(&bg, &src, 2*nb, nb);
    t_bucket_build = wtime() - t0;
    t0 = wtime();
    bucket_nearest(&bg, &tgt, nn_bucket);
    t_bucket = wtime() - t0;
    printf("  lat-lon buckets  build  %8.4f s   search %8.4f s   (%dx%d)\n",
           t_bucket_build, t_bucket, 2*nb, nb);

    kdtree kt;
    t0 = wtime();
    kdtree_build(&kt, &src);
    t_kd_build = wtime() - t0;
    t0 = wtime();
    kdtree_nearest(&kt, &tgt, nn_kd);
    t_kd = wtime() - t0;
    printf("  k-d tree         build  %8.4f s   search %8.4f s\n",
           t_kd_build, t_kd);

    printf("\n  --- correctness (must be EXACT, not approximate) ---\n");
    int ok = 1;
    if (nn_brute[0] < 0) {
        printf("     brute force is still a stub -- start with TODO 2.\n");
        ok = 0;
    } else {
        ok &= compare_results(nn_bucket, nn_brute, n_tgt, "buckets");
        ok &= compare_results(nn_kd,     nn_brute, n_tgt, "k-d tree");
    }

    /* Interpolation check: a constant field must survive IDW exactly. */
    double *f_src = malloc((size_t)n_src * sizeof(double));
    double *f_tgt = malloc((size_t)n_tgt * sizeof(double));
    for (int i = 0; i < n_src; i++) f_src[i] = 1.0;
    idw_interpolate(&src, &tgt, f_src, f_tgt, k);
    double maxerr = 0.0;
    for (int j = 0; j < n_tgt; j++) {
        double e = fabs(f_tgt[j] - 1.0);
        if (e > maxerr) maxerr = e;
    }
    printf("     IDW consistency  max |F-1| = %.3e  %s\n", maxerr,
           maxerr < 1e-14 ? "PASS" : "FAIL");
    if (maxerr >= 1e-14) ok = 0;

    printf("\n  %s\n", ok ? "All checks passed."
                          : "Work the TODOs in order: 1 -> 2 -> 3 -> 4 -> 5.");

    free(nn_brute); free(nn_bucket); free(nn_kd);
    free(f_src); free(f_tgt);
    bucket_free(&bg); kdtree_free(&kt);
    free_pointset(&src); free_pointset(&tgt);
    return ok ? 0 : 1;
}

/* ==========================================================================
 * TODO 6 — write your answers in notes/day3.md
 *
 * (a) Run `make scaling` (n = 10k, 40k, 160k, 640k). Fit the exponent for
 *     each method's search time. Do you measure O(n^2), O(n), O(n log n)?
 *     Where does the k-d tree BUILD time cross over and start to matter?
 *
 * (b) The bucket method uses a lat-lon grid, whose cells become slivers near
 *     the poles. Measure it: report the mean number of candidate points
 *     examined for targets with |lat| > 80 degrees versus |lat| < 10.
 *     Then name two fixes (search terms: equal-area binning, HEALPix).
 *
 * (c) Why can a plain 3D k-d tree be exactly correct for a great-circle
 *     nearest-neighbour query? State the property of the map from chord
 *     distance to arc distance that makes it work, and give a query where
 *     the trick would BREAK (hint: what if the points were not all on the
 *     same sphere?).
 *
 * (d) IDW is consistent but not conservative; Exercise 10's operator is
 *     both. For each of these coupled fields, say which you would use and
 *     why: sea-surface temperature; net surface heat flux; wind stress;
 *     river runoff into the ocean.
 *
 * (e) Parallel version: the source points live on the atmosphere's ranks
 *     and the targets on the ocean's, and no rank holds either grid whole.
 *     Sketch a distributed nearest-neighbour search. What is the minimum
 *     communication, and what do you replicate? (Search term: distributed
 *     spatial join. This is the hard part of coupler setup.)
 *
 * (f) DKRZ context: coupler setup for a high-resolution coupled run can
 *     take longer than a simulated day of model time. Given that the grids
 *     do not change, what is the obvious engineering answer, and what has
 *     to be true for it to be safe? (Think about what invalidates a cached
 *     weight file, and how you would detect that automatically.)
 * ========================================================================== */
