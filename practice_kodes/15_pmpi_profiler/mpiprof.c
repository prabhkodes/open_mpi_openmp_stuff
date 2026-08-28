/* ==========================================================================
 * EXERCISE 15 — A PMPI interposition profiler
 * ==========================================================================
 *
 * GOAL
 *   Write a profiling library that intercepts every MPI call an application
 *   makes, with ZERO changes to that application's source, and reports call
 *   counts, time, message volume, a message-size histogram, and a
 *   load-imbalance estimate.
 *
 * WHY (DKRZ)
 *   This is the highest-leverage item on the whole week's list.
 *
 *   The job is optimising code you did not write. A scientist hands you a
 *   500k-line Fortran model and says it is slow at 512 nodes. You cannot
 *   read it all, you cannot instrument it by hand, and rebuilding it with a
 *   heavyweight tool may not even be possible on that machine that week.
 *   What you CAN do is link one object file in and get a complete picture of
 *   its communication behaviour in a single run.
 *
 *   This is also exactly how Score-P, mpiP, IPM, and Vampir work
 *   underneath. Building a small one means that when you use the big ones
 *   you know what they can and cannot see -- and in an interview, "I wrote
 *   my own PMPI profiler" lands very differently from "I have used Score-P".
 *
 * HOW IT WORKS
 *   The MPI standard requires every MPI_Xxx to also be callable as PMPI_Xxx.
 *   So you define your OWN MPI_Send, do your bookkeeping, and call
 *   PMPI_Send to do the real work. At link time your definition wins over
 *   the library's, and the application -- which called MPI_Send -- lands in
 *   yours without knowing.
 *
 *   Link-time interposition is used here because it is portable. macOS has
 *   no LD_PRELOAD (it has DYLD_INSERT_LIBRARIES, restricted by SIP), so
 *   linking the object file directly is the approach that works everywhere.
 *
 * TASKS
 *   TODO 1  prof_record       — the accounting helper
 *   TODO 2  wrap Send/Recv/Isend/Irecv/Sendrecv
 *   TODO 3  wrap Barrier/Allreduce/Bcast/Alltoallv
 *   TODO 4  size histogram in log2 bins
 *   TODO 5  report at Finalize, aggregated across ranks
 *   TODO 6  estimate load imbalance from time spent in synchronising calls
 *   TODO 7  answer the questions at the bottom
 *
 * ACCEPTANCE
 *   - `make run` profiles victim.c and the counts match what victim.c says
 *     it did (it prints its own expected totals -- they must agree exactly)
 *   - the histogram shows the two distinct message sizes victim.c uses
 *   - the imbalance estimate correctly fingers the rank victim.c delays
 *   - you then run it on Exercise 06 and Exercise 09 unchanged
 *
 * HINTS
 *   - MPI_Type_size gives bytes per element; multiply by count.
 *   - Time only the PMPI call itself. Your bookkeeping must be outside the
 *     timed region or you profile your own profiler.
 *   - Time in MPI_Barrier and MPI_Allreduce is mostly WAITING for other
 *     ranks, not communicating. That is the signal for TODO 6: a rank that
 *     spends little time in collectives is the slow one everyone waits for.
 *   - Do the cross-rank reduction BEFORE calling PMPI_Finalize. After it,
 *     MPI is gone.
 *   - Non-blocking calls are the subtle case: MPI_Isend returns immediately,
 *     so its "time" is meaningless. The cost is in MPI_Wait. Record the
 *     BYTES at Isend and the TIME at Wait, and say so in your report.
 * ========================================================================== */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- the call table ---------------------------------------------------- */
enum {
    P_SEND, P_RECV, P_ISEND, P_IRECV, P_SENDRECV,
    P_WAIT, P_WAITALL,
    P_BARRIER, P_ALLREDUCE, P_BCAST, P_ALLTOALL, P_ALLTOALLV,
    P_NCALLS
};

static const char *call_name[P_NCALLS] = {
    "MPI_Send", "MPI_Recv", "MPI_Isend", "MPI_Irecv", "MPI_Sendrecv",
    "MPI_Wait", "MPI_Waitall",
    "MPI_Barrier", "MPI_Allreduce", "MPI_Bcast", "MPI_Alltoall", "MPI_Alltoallv"
};

/* Which calls are SYNCHRONISING -- time here is mostly waiting for peers,
 * and is the basis of the load-imbalance estimate in TODO 6. */
static const int is_sync[P_NCALLS] = {
    0, 1, 0, 0, 1,
    1, 1,
    1, 1, 1, 1, 1
};

#define NBINS 24        /* log2 bins: <1B, <2B, <4B ... up to 8 MiB+ */

typedef struct {
    long   count;
    double time;
    double bytes;
} callstat;

static callstat stats[P_NCALLS];
static long     hist[NBINS];
static double   t_start_run;
static int      prof_enabled = 1;

/* --------------------------------------------------------------------------
 * TODO 1: the accounting helper.
 *
 *   stats[id].count += 1
 *   stats[id].time  += dt
 *   stats[id].bytes += bytes
 *   and bin `bytes` into the histogram (TODO 4)
 *
 * Keep it branch-light: this runs on every single MPI call in the
 * application, and a profiler that measurably slows the program down is
 * reporting a program that no longer exists.
 * ------------------------------------------------------------------------ */
static void prof_record(int id, double dt, double bytes)
{
    if (!prof_enabled) return;
    /* TODO 1 */
    (void)id; (void)dt; (void)bytes;
}

/* --------------------------------------------------------------------------
 * TODO 4: bin a message size into log2 buckets.
 *
 * bin 0 = 0 bytes, bin k = sizes in [2^(k-1), 2^k). Clamp at NBINS-1.
 * The histogram is the single most useful plot a communication profiler
 * produces: it instantly separates "latency-bound, millions of tiny
 * messages" from "bandwidth-bound, a few big ones", and those two problems
 * have completely different fixes.
 * ------------------------------------------------------------------------ */
static int size_bin(double bytes)
{
    /* TODO 4 */
    (void)bytes;
    return 0;
}

static double bytes_of(int count, MPI_Datatype dt)
{
    int sz = 0;
    PMPI_Type_size(dt, &sz);
    return (double)count * (double)sz;
}

/* ==========================================================================
 * TODO 2 — point-to-point wrappers.
 *
 * The pattern for every one of them:
 *
 *     int MPI_Send(...) {
 *         double t = PMPI_Wtime();
 *         int rc = PMPI_Send(...same args...);
 *         prof_record(P_SEND, PMPI_Wtime() - t, bytes_of(count, datatype));
 *         return rc;
 *     }
 *
 * MPI_Send is written out in full as the worked example. Do the rest.
 * ========================================================================== */

int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest,
             int tag, MPI_Comm comm)
{
    double t = PMPI_Wtime();
    int rc = PMPI_Send(buf, count, datatype, dest, tag, comm);
    prof_record(P_SEND, PMPI_Wtime() - t, bytes_of(count, datatype));
    return rc;
}

/* TODO 2a: MPI_Recv
 *   int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source,
 *                int tag, MPI_Comm comm, MPI_Status *status)               */

/* TODO 2b: MPI_Isend and MPI_Irecv.
 *   Record the BYTES here but expect the TIME to be ~0 -- these return
 *   immediately. The real cost shows up in MPI_Wait/MPI_Waitall.           */

/* TODO 2c: MPI_Sendrecv
 *   Count the bytes of BOTH halves. Decide whether to record it as one call
 *   or two, and be consistent with what you print.                          */

/* TODO 2d: MPI_Wait and MPI_Waitall
 *   Time only -- the bytes were already counted at Isend/Irecv time.
 *   Double-counting them here would inflate your volume by 2x, which is a
 *   very easy mistake to make and a very embarrassing one to present.      */

/* ==========================================================================
 * TODO 3 — collective wrappers.
 *
 * MPI_Barrier is the interesting one: it moves no data at all, so 100% of
 * its time is waiting. A rank with a SMALL barrier time arrived LAST, which
 * means it is the slow rank. That inversion trips people up constantly --
 * say it out loud until it sticks.
 * ========================================================================== */

/* TODO 3a: MPI_Barrier(MPI_Comm comm)                                      */
/* TODO 3b: MPI_Allreduce(const void*, void*, int, MPI_Datatype, MPI_Op,
 *                        MPI_Comm)                                          */
/* TODO 3c: MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm)              */
/* TODO 3d: MPI_Alltoall and MPI_Alltoallv                                  */

/* ==========================================================================
 * Init / Finalize: where the report happens.
 * ========================================================================== */

int MPI_Init(int *argc, char ***argv)
{
    int rc = PMPI_Init(argc, argv);
    memset(stats, 0, sizeof(stats));
    memset(hist,  0, sizeof(hist));
    t_start_run = PMPI_Wtime();
    return rc;
}

int MPI_Init_thread(int *argc, char ***argv, int required, int *provided)
{
    int rc = PMPI_Init_thread(argc, argv, required, provided);
    memset(stats, 0, sizeof(stats));
    memset(hist,  0, sizeof(hist));
    t_start_run = PMPI_Wtime();
    return rc;
}

/* --------------------------------------------------------------------------
 * TODO 5 + TODO 6: the report.
 *
 * Everything here must happen BEFORE PMPI_Finalize -- afterwards there is no
 * MPI left to reduce with.
 *
 * TODO 5: reduce across ranks and print on rank 0:
 *           - per call: total count, total time, total bytes, mean size
 *           - the message-size histogram, summed over ranks
 *           - MPI time as a percentage of total run time
 *
 * TODO 6: the imbalance estimate.
 *           - each rank's total time in SYNCHRONISING calls (is_sync[])
 *           - MPI_Reduce with MPI_MIN and MPI_MAX to find the extremes
 *           - the rank with the LOWEST sync time is the slowest rank:
 *             it arrives last, so it never waits
 *           - report  (max_sync - min_sync)  as the imbalance, and use
 *             MPI_MINLOC to name the guilty rank
 * ------------------------------------------------------------------------ */
int MPI_Finalize(void)
{
    int rank = 0, nprocs = 1, i;
    double t_run = PMPI_Wtime() - t_start_run;

    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    prof_enabled = 0;   /* stop profiling our own reduction calls */

    double tot_time[P_NCALLS], tot_bytes[P_NCALLS];
    long   tot_count[P_NCALLS];
    long   tot_hist[NBINS];
    double loc_time[P_NCALLS], loc_bytes[P_NCALLS];
    long   loc_count[P_NCALLS];

    for (i = 0; i < P_NCALLS; i++) {
        loc_time[i]  = stats[i].time;
        loc_bytes[i] = stats[i].bytes;
        loc_count[i] = stats[i].count;
    }

    PMPI_Reduce(loc_time,  tot_time,  P_NCALLS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    PMPI_Reduce(loc_bytes, tot_bytes, P_NCALLS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    PMPI_Reduce(loc_count, tot_count, P_NCALLS, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
    PMPI_Reduce(hist,      tot_hist,  NBINS,    MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);

    /* ---- TODO 6: imbalance ---------------------------------------------- */
    double my_sync = 0.0;
    for (i = 0; i < P_NCALLS; i++) if (is_sync[i]) my_sync += stats[i].time;

    struct { double val; int rank; } smin, smax, in;
    in.val = my_sync; in.rank = rank;
    PMPI_Reduce(&in, &smin, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
    PMPI_Reduce(&in, &smax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double mpi_time = 0.0, mpi_bytes = 0.0;
        long   mpi_calls = 0;
        for (i = 0; i < P_NCALLS; i++) {
            mpi_time  += tot_time[i];
            mpi_bytes += tot_bytes[i];
            mpi_calls += tot_count[i];
        }

        printf("\n");
        printf("================ MPI PROFILE (%d ranks) ================\n", nprocs);
        printf("  wall time (rank 0)   : %10.4f s\n", t_run);
        printf("  total MPI time       : %10.4f s  (%.1f %% of %d x wall)\n",
               mpi_time, 100.0 * mpi_time / (t_run * nprocs), nprocs);
        printf("  total MPI calls      : %10ld\n", mpi_calls);
        printf("  total bytes moved    : %10.3f MiB\n", mpi_bytes / 1048576.0);
        if (mpi_calls == 0)
            printf("  (zero calls recorded -- prof_record is still a stub, TODO 1)\n");
        printf("\n");
        printf("  %-16s %10s %12s %14s %12s\n",
               "call", "count", "time (s)", "bytes", "mean size");
        printf("  ---------------------------------------------------------------------\n");
        for (i = 0; i < P_NCALLS; i++) {
            if (tot_count[i] == 0) continue;
            printf("  %-16s %10ld %12.5f %14.0f %12.1f\n",
                   call_name[i], tot_count[i], tot_time[i], tot_bytes[i],
                   tot_bytes[i] / (double)tot_count[i]);
        }

        printf("\n  message size histogram (all ranks)\n");
        printf("  ----------------------------------\n");
        long hmax = 1;
        for (i = 0; i < NBINS; i++) if (tot_hist[i] > hmax) hmax = tot_hist[i];
        for (i = 0; i < NBINS; i++) {
            if (tot_hist[i] == 0) continue;
            int bar = (int)(40.0 * tot_hist[i] / hmax);
            char label[32];
            if (i == 0) snprintf(label, sizeof label, "0 B");
            else snprintf(label, sizeof label, "< %ld B", 1L << i);
            printf("  %-12s %8ld  ", label, tot_hist[i]);
            for (int b = 0; b < bar; b++) putchar('#');
            putchar('\n');
        }

        printf("\n  load imbalance (from time in synchronising calls)\n");
        printf("  ------------------------------------------------\n");
        printf("    most  time waiting : rank %d  (%.4f s)  <- fastest rank\n",
               smax.rank, smax.val);
        printf("    least time waiting : rank %d  (%.4f s)  <- SLOWEST rank\n",
               smin.rank, smin.val);
        printf("    imbalance          : %.4f s\n", smax.val - smin.val);
        if (smax.val - smin.val <= 0.0)
            printf("    (TODO 6: zero -- are the collective wrappers written?)\n");
        printf("=======================================================\n\n");
    }

    return PMPI_Finalize();
}

/* ==========================================================================
 * TODO 7 — write your answers in notes/day4.md
 *
 * (a) Run your profiler on Exercise 06 and Exercise 09 with NO source
 *     changes to either. Paste the profiles into your notes. For Exercise 09,
 *     does the imbalance report correctly identify the atmosphere as the
 *     slow component?
 *
 * (b) Why is a LOW barrier time the signature of a SLOW rank? Draw the
 *     timeline for 3 ranks where rank 1 takes twice as long, and mark each
 *     rank's time-in-barrier.
 *
 * (c) Non-blocking calls: MPI_Isend records ~0 time. So where does the cost
 *     of a non-blocking exchange actually appear in your profile, and what
 *     would you have to add to distinguish "the message took a long time"
 *     from "we waited a long time before needing it"? (This is precisely
 *     what a real tool's late-sender analysis does.)
 *
 * (d) Your profiler sees calls, not phases. A model has a dynamics phase and
 *     a physics phase with very different communication. Design a minimal
 *     API the application could call to mark phases, and say what it costs
 *     you in the "zero source changes" property you started with.
 *
 * (e) Overhead: time a run with and without the profiler linked in. What is
 *     the percentage cost? At what call rate would this profiler start to
 *     distort what it measures?
 *
 * (f) Compare with Score-P/mpiP: name two things they do that this cannot,
 *     and one thing this does better. (Hint for the last one: think about
 *     what it takes to deploy either on a cluster you just got an account on.)
 * ========================================================================== */
