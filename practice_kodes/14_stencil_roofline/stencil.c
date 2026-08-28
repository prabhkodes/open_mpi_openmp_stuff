/* ==========================================================================
 * EXERCISE 14 — Roofline analysis of a memory-bound stencil
 * ==========================================================================
 *
 * GOAL
 *   Measure your machine's actual peak bandwidth and peak FLOP rate, compute
 *   a 3D stencil's arithmetic intensity, place it on the roofline, and then
 *   optimise it and watch the point move. Establish the ceiling FIRST, so
 *   you know when to stop optimising.
 *
 * WHY (DKRZ)
 *   You have done roofline before (low_level_optimisations/), so this is a
 *   refresher with a different purpose: the discipline of knowing the
 *   ceiling before you touch the code.
 *
 *   The dycore of an ESM is memory-bound. Almost every kernel is a
 *   low-intensity stencil over a huge array, and the single most useful
 *   thing you can tell a scientist is "this kernel is already at 85% of
 *   achievable bandwidth, so stop tuning it and go look at the halo
 *   exchange instead." That sentence saves weeks, and you can only say it
 *   if you measured the roof.
 *
 *   The related trap: comparing to the VENDOR's peak bandwidth. Nobody ever
 *   reaches it. Measure your own with a STREAM-style benchmark and compare
 *   to that — the gap between vendor peak and STREAM peak is 20-40% and
 *   pretending otherwise makes every kernel look bad.
 *
 * TASKS
 *   TODO 1  measure_bandwidth   — STREAM triad, your actual roof
 *   TODO 2  measure_peak_flops  — FMA-saturating loop, the other roof
 *   TODO 3  stencil_naive       — the baseline kernel
 *   TODO 4  compute the arithmetic intensity ON PAPER, then verify
 *   TODO 5  stencil_blocked     — cache blocking, and measure the move
 *   TODO 6  answer the questions at the bottom
 *
 * ACCEPTANCE
 *   - a measured bandwidth within 2x of your hardware's spec (if it is 10x
 *     off, your timing or your byte count is wrong -- find out which)
 *   - an arithmetic intensity you DERIVED before measuring, matching the
 *     measurement
 *   - the blocked version is measurably faster, and you can say which
 *     memory level it started hitting
 *   - a statement of what fraction of achievable bandwidth you reached
 *
 * HINTS
 *   - Arithmetic intensity = FLOPs / bytes moved from DRAM. The subtle part
 *     is "from DRAM": a value already in cache costs nothing. For a 7-point
 *     stencil sweeping a large array, count each element as loaded ONCE if
 *     your blocking works, and up to 3 times if it does not. That factor of
 *     3 IS the optimisation.
 *   - Use `volatile` or a printed checksum to stop the compiler deleting
 *     benchmark loops. A "1000 GB/s" result means it was optimised away.
 *   - Allocate with aligned_alloc for 64-byte alignment, and touch the
 *     arrays with the same thread that will use them (first-touch NUMA).
 *   - The blocked version should block in the two SLOW dimensions (j, k)
 *     and stream the fast one (i). Blocking the fast dimension destroys
 *     the vectorisation you already had.
 * ========================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static double wtime(void)
{
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
#endif
}

static double *alloc_aligned(size_t n)
{
    void *p = NULL;
    if (posix_memalign(&p, 64, n * sizeof(double)) != 0) {
        fprintf(stderr, "allocation of %zu doubles failed\n", n);
        exit(1);
    }
    return (double *)p;
}

/* --------------------------------------------------------------------------
 * TODO 1: STREAM triad -- a[i] = b[i] + s*c[i]
 *
 * Three arrays, each much larger than last-level cache. Bytes moved per
 * iteration: 8 (read b) + 8 (read c) + 8 (write a) = 24.
 *
 * Careful, and this is the classic subtlety: on a write-allocate machine
 * (i.e. essentially all of them) writing a[i] first READS the cache line
 * it is about to overwrite, so the true traffic is 32 bytes, not 24.
 * STREAM traditionally reports 24. Pick one, say which, and be consistent.
 *
 * Return GB/s.
 * ------------------------------------------------------------------------ */
double measure_bandwidth(size_t n, int reps)
{
    double *a = alloc_aligned(n), *b = alloc_aligned(n), *c = alloc_aligned(n);

    /* First touch in parallel: on a NUMA machine this decides which socket
     * each page lives on, and getting it wrong halves your bandwidth. */
    #pragma omp parallel for
    for (size_t i = 0; i < n; i++) { a[i] = 0.0; b[i] = 1.0; c[i] = 2.0; }

    double best = 0.0;
    for (int r = 0; r < reps; r++) {
        double t0 = wtime();
        /* TODO 1a: the triad loop, OpenMP-parallel */
        double t = wtime() - t0;
        /* TODO 1b: bytes = 24 * n (or 32 -- decide and document);
         *          gbs = bytes / t / 1e9; keep the best */
        (void)t;
    }

    /* Consume the result so the loop cannot be deleted. */
    volatile double sink = a[n / 2];
    (void)sink;

    free(a); free(b); free(c);
    return best;
}

/* --------------------------------------------------------------------------
 * TODO 2: peak FLOP rate.
 *
 * A loop of independent FMAs with everything in registers, so nothing
 * touches memory. You need enough independent accumulators to cover the
 * FMA latency -- typically 8 to 16 on a modern core. With too few, you
 * measure latency instead of throughput and get a number several times too
 * low.
 *
 * Count 2 FLOPs per FMA. Return GFLOP/s.
 * ------------------------------------------------------------------------ */
double measure_peak_flops(long iters)
{
    double best = 0.0;
    #pragma omp parallel reduction(max: best)
    {
        double acc[16];
        for (int k = 0; k < 16; k++) acc[k] = (double)k;
        const double m = 1.0000001, p = 0.0000001;

        double t0 = wtime();
        /* TODO 2: for (long it = 0; it < iters; it++)
         *             for (k = 0; k < 16; k++) acc[k] = acc[k]*m + p;      */
        double t = wtime() - t0;

        double sum = 0.0;
        for (int k = 0; k < 16; k++) sum += acc[k];
        if (sum == 12345.6789) printf(" ");   /* keep the loop alive */

        /* TODO 2b: flops = 2.0 * 16 * iters;  best = flops / t / 1e9;      */
        (void)t; (void)iters;
    }
    return best;
}

/* ==========================================================================
 * The kernel under study: a 7-point 3D stencil, the shape of every diffusion
 * and pressure-gradient term in a dycore.
 * ========================================================================== */

#define IDX(i, j, k, nx, ny) ((size_t)(k) * (ny) * (nx) + (size_t)(j) * (nx) + (i))

/* --------------------------------------------------------------------------
 * TODO 3: the naive version. Straight triple loop, i innermost.
 *
 *   out[i,j,k] = c0*in[i,j,k]
 *              + c1*(in[i-1,j,k] + in[i+1,j,k]
 *                  + in[i,j-1,k] + in[i,j+1,k]
 *                  + in[i,j,k-1] + in[i,j,k+1])
 *
 * Parallelise over k with OpenMP. Interior points only.
 * ------------------------------------------------------------------------ */
void stencil_naive(const double *in, double *out, int nx, int ny, int nz)
{
    const double c0 = 0.5, c1 = 1.0 / 12.0;
    /* TODO 3 */
    (void)in; (void)out; (void)nx; (void)ny; (void)nz; (void)c0; (void)c1;
}

/* --------------------------------------------------------------------------
 * TODO 5: cache-blocked version.
 *
 * The naive sweep touches three j-planes at a time. If one j-plane is
 * bigger than L2, the plane you loaded for j is evicted before j+1 needs
 * it, so every element gets fetched from DRAM three times instead of once.
 *
 * Block over j and k (NOT over i -- keep the innermost loop long and
 * stride-1 so it still vectorises):
 *
 *   for (kk = 1; kk < nz-1; kk += BK)
 *     for (jj = 1; jj < ny-1; jj += BJ)
 *       for (k = kk; k < min(kk+BK, nz-1); k++)
 *         for (j = jj; j < min(jj+BJ, ny-1); j++)
 *           for (i = 1; i < nx-1; i++)   ... same body ...
 *
 * Then sweep BJ and BK (make blocksweep) to find the pair that keeps the
 * working set in L2. Predict the best value from your cache size first.
 * ------------------------------------------------------------------------ */
void stencil_blocked(const double *in, double *out, int nx, int ny, int nz,
                     int bj, int bk)
{
    const double c0 = 0.5, c1 = 1.0 / 12.0;
    /* TODO 5 */
    (void)in; (void)out; (void)nx; (void)ny; (void)nz;
    (void)bj; (void)bk; (void)c0; (void)c1;
}

static double checksum(const double *a, size_t n)
{
    double s = 0.0;
    for (size_t i = 0; i < n; i++) s += a[i];
    return s;
}

int main(int argc, char **argv)
{
    int nx = 256, ny = 256, nz = 256;
    int bj = 16, bk = 16, reps = 5;
    if (argc > 1) nx = ny = nz = atoi(argv[1]);
    if (argc > 2) bj = atoi(argv[2]);
    if (argc > 3) bk = atoi(argv[3]);

    size_t n = (size_t)nx * ny * nz;
    printf("=== Exercise 14: roofline of a 7-point stencil ===\n");
    printf("grid %d^3 = %.2f Mcells, %.1f MiB per array\n",
           nx, n / 1e6, n * 8.0 / 1048576.0);
#ifdef _OPENMP
    printf("threads: %d\n", omp_get_max_threads());
#endif
    printf("\n");

    /* ---- establish the roofs BEFORE looking at the kernel --------------- */
    printf("  --- machine roofs (measure these first, always) ---\n");
    double bw = measure_bandwidth(1u << 24, 5);       /* 16M doubles = 128 MiB */
    double pf = measure_peak_flops(2000000L);
    printf("    STREAM triad bandwidth : %8.2f GB/s   %s\n", bw,
           bw <= 0.0 ? "<-- TODO 1" : "");
    printf("    peak FLOP rate         : %8.2f GFLOP/s %s\n", pf,
           pf <= 0.0 ? "<-- TODO 2" : "");
    if (bw > 0.0 && pf > 0.0)
        printf("    machine balance        : %8.2f FLOP/byte  (ridge point)\n",
               pf / bw);
    printf("\n");

    /* ---- the kernel ----------------------------------------------------- */
    double *in  = alloc_aligned(n);
    double *out = alloc_aligned(n);
    #pragma omp parallel for
    for (size_t i = 0; i < n; i++) { in[i] = (double)(i % 97) * 0.01; out[i] = 0.0; }

    double t0, t_naive, t_blocked;
    double chk_naive, chk_blocked;

    t0 = wtime();
    for (int r = 0; r < reps; r++) stencil_naive(in, out, nx, ny, nz);
    t_naive = (wtime() - t0) / reps;
    chk_naive = checksum(out, n);

    memset(out, 0, n * sizeof(double));
    t0 = wtime();
    for (int r = 0; r < reps; r++) stencil_blocked(in, out, nx, ny, nz, bj, bk);
    t_blocked = (wtime() - t0) / reps;
    chk_blocked = checksum(out, n);

    /* ---- TODO 4: the arithmetic intensity ------------------------------- */
    /* Interior points: (nx-2)*(ny-2)*(nz-2)
     * FLOPs per point:  6 adds + 1 add + 2 muls  = 9      (count them!)
     * Bytes per point:  BEST case  8 (read in) + 8 (write out) = 16
     *                   WORST case in is read ~3x -> 32
     * So AI is somewhere between 9/32 = 0.28 and 9/16 = 0.56 FLOP/byte.
     * Both are far left of the ridge point, which is WHY this kernel is
     * bandwidth-bound and why cache blocking (which moves you from the
     * worst case toward the best) is the optimisation that matters.        */
    double interior = (double)(nx - 2) * (ny - 2) * (nz - 2);
    double flops = 9.0 * interior;
    double bytes_best = 16.0 * interior;

    printf("  --- kernel ---\n");
    if (t_naive > 1e-9 && chk_naive != 0.0) {
        printf("    naive    %8.4f s   %7.2f GFLOP/s   %7.2f GB/s (best-case bytes)\n",
               t_naive, flops / t_naive / 1e9, bytes_best / t_naive / 1e9);
    } else {
        printf("    naive    -- stub (TODO 3)\n");
    }
    if (t_blocked > 1e-9 && chk_blocked != 0.0) {
        printf("    blocked  %8.4f s   %7.2f GFLOP/s   %7.2f GB/s   (bj=%d bk=%d)\n",
               t_blocked, flops / t_blocked / 1e9,
               bytes_best / t_blocked / 1e9, bj, bk);
        printf("    speedup over naive     : %.2fx\n", t_naive / t_blocked);
    } else {
        printf("    blocked  -- stub (TODO 5)\n");
    }

    if (chk_naive != 0.0 && chk_blocked != 0.0) {
        double rel = fabs(chk_naive - chk_blocked) / fabs(chk_naive);
        printf("    checksum agreement     : %.3e  %s\n", rel,
               rel < 1e-12 ? "PASS" : "FAIL -- blocking changed the answer");
    }

    if (bw > 0.0 && t_blocked > 1e-9 && chk_blocked != 0.0) {
        printf("\n    fraction of achievable bandwidth: %.1f %%\n",
               100.0 * (bytes_best / t_blocked / 1e9) / bw);
        printf("    (above ~80%% -- stop tuning this kernel and go look\n");
        printf("     somewhere else. That judgement IS the job.)\n");
    }

    free(in); free(out);
    return 0;
}

/* ==========================================================================
 * TODO 6 — write your answers in notes/day4.md
 *
 * (a) Derive the arithmetic intensity by hand before running anything.
 *     Count the FLOPs in the stencil body exactly (is c1*(sum of 6) 7 flops
 *     or 12?). Then give the AI for the best and worst caching cases, and
 *     say where each sits relative to your measured ridge point.
 *
 * (b) Run `make blocksweep`. Which (bj, bk) wins? Compute the working-set
 *     size for that block -- 3 j-planes of bj*nx doubles -- and compare it
 *     to your L2. Does the winner match your prediction?
 *
 * (c) You measured bandwidth with STREAM triad. Now measure it with a pure
 *     read loop (sum an array) and a pure write loop. Are the three the
 *     same? Explain any difference in terms of write-allocate and of read
 *     vs write buffers in the memory controller.
 *
 * (d) Your stencil reaches some fraction of STREAM. If it is below 50%,
 *     name three plausible causes and design a measurement that
 *     distinguishes them. (Candidates: unblocked traffic, TLB misses at
 *     this array size, NUMA first-touch, insufficient outstanding loads.)
 *
 * (e) Try huge pages if your platform has them. At 256^3 doubles the array
 *     is 128 MiB, which is ~32000 4 KiB pages -- far more than the TLB
 *     holds. How much does that cost, and how would you confirm it?
 *
 * (f) DKRZ context: an ICON dycore kernel on an UNSTRUCTURED grid cannot do
 *     this blocking, because neighbours are reached through an index array
 *     rather than i+1. What is the equivalent optimisation there? (You have
 *     already met the answer, in Exercises 07 and 08 -- connect them.)
 * ========================================================================== */
