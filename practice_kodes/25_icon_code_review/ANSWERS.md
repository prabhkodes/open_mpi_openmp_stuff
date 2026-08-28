# Exercise 25 — answer key

**Write your own answer first.** The gap between yours and this is your study list.

Where a number depends on hardware, the arithmetic is shown so you can redo it
for whatever machine you are asked about. Getting the *method* right matters
more than memorising a factor.

---

# Part A — memory coalescing and data layout

## Q1. The loop nest

**1. What is wrong.**

The loop order is inverted for the GPU. With `f(nproma, nlev, nblks)`, the
address of element `f(jc,jk,jb)` is

```
base + [ (jc-1) + (jk-1)*nproma + (jb-1)*nproma*nlev ] * 8
```

so `jc` is the **fastest-varying** index — consecutive `jc` are consecutive in
memory. `COLLAPSE(2)` linearises the two loops with the **inner** loop varying
fastest, and here the inner loop is `jk`. So consecutive threads get consecutive
`jk` at fixed `jc`, and stride through memory by `nproma` elements.

**2. Address stride.**

Consecutive threads differ by `nproma * 8 = 20000 * 8 = 160,000 bytes`.

**3. Cache lines per warp.**

- Ideal (`jc` innermost): 32 threads × 8 B = 256 B contiguous → **8 sectors** of
  32 B, or **2** lines of 128 B.
- Actual: 32 threads, each 160 kB from the next → **32 distinct sectors**
  (1024 B fetched to use 256 B) and **32 distinct 128 B lines**.

So **4× more memory transactions** at 32-byte sector granularity, and a **16×**
larger cache-line footprint. Quote it as "4–16× depending on access granularity"
and say which granularity you mean — that is the answer that shows you know the
hardware rather than the slogan.

**4. Fixed version.**

```fortran
!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
DO jk = 1, nlev
  DO jc = i_startidx, i_endidx
    z_theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) &
      &                 * ( 1._wp + vwp1*p_prog%tracer(jc,jk,jb,iqv) &
      &                           - p_prog%tracer(jc,jk,jb,iqc) )
  END DO
END DO
!$ACC END PARALLEL
```

Only the two `DO` lines swapped. This is the single most common ICON GPU
review finding.

**5. Is it wrong on CPU too?**

Yes, for a *different* mechanism. At `nproma = 256` the inner `jk` loop strides
by `256 * 8 = 2048 B` — a new cache line every iteration, using 8 of every 64
bytes, and **not stride-1, so it will not vectorise**.

The important observation: *the same fix helps both machines*. That is what
makes one source viable across CPU and GPU, and it is why ICON's convention is
always `jk` outer, `jc` inner. The thing that differs between machines is
`nproma` (a namelist value), not the loop order.

> Interview line: "The loop order is a correctness-of-performance property of
> the layout, not of the machine. `nproma` is the machine-specific knob."

---

## Q2. The neighbour gather

**1. Why it may still not coalesce.**

`jc` is innermost, so the *thread-to-index* mapping is right. But the address
actually read is `flux(iidx(jc,jb,1), ...)` — it depends on the **value** of
`iidx`, not on `jc`. Coalescing therefore depends entirely on whether cells
with adjacent `jc` have neighbours with adjacent indices, which is a property of
the **cell numbering**, decided at setup.

**2. Loads per output element.**

12 loads:

| Load | Coalesces? |
|---|---|
| `iidx(jc,jb,1..3)` — 3 | yes, consecutive `jc` |
| `iblk(jc,jb,1..3)` — 3 | yes |
| `div_coeff(jc,jb,1..3)` — 3 | yes |
| `flux(iidx(...), jk, iblk(...))` — 3 | **depends on the numbering** |

Nine of twelve are perfect; the three that matter are the gather.

**3. The setup-time decision.**

The **cell renumbering / partitioning**. A space-filling-curve ordering keeps
geometric neighbours close in index space, so the gather is near-contiguous; a
random ordering makes every gather a separate sector. This is exactly what
[Exercise 07](../07_partitioning) measures (edge cut, locality) and
[Exercise 08](../08_nproma_blocking) TODO 5 quantifies (`make permute` —
sequential vs randomised numbering, same arithmetic, very different time).

**4. Caching the neighbour values in a local array.**

As stated, no — there is no reuse *within* one output element; each of the three
values is read once. A per-thread local array of any size risks spilling to
local memory (which is global memory), making it worse.

But there **is** real reuse, and it is in the other direction: `iidx(jc,jb,n)`
and `div_coeff(jc,jb,n)` do **not depend on `jk`**. With `COLLAPSE(2)` you
re-read them for all 90 levels. Restructure so `jc` is the vector loop and `jk`
is a sequential inner loop per thread, and each thread loads the six index
values and three coefficients **once** and reuses them across all levels:

```fortran
!$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
DO jc = i_startidx, i_endidx
  i1 = iidx(jc,jb,1);  b1 = iblk(jc,jb,1);  c1 = div_coeff(jc,jb,1)
  i2 = iidx(jc,jb,2);  b2 = iblk(jc,jb,2);  c2 = div_coeff(jc,jb,2)
  i3 = iidx(jc,jb,3);  b3 = iblk(jc,jb,3);  c3 = div_coeff(jc,jb,3)
  !$ACC LOOP SEQ
  DO jk = 1, nlev
    p_diag%div(jc,jk,jb) = flux(i1,jk,b1)*c1 + flux(i2,jk,b2)*c2 + flux(i3,jk,b3)*c3
  END DO
END DO
```

This cuts index traffic by `nlev` (90×) and **stays coalesced**, because
consecutive threads still hold consecutive `jc`. The cost is register pressure
(9 extra live values per thread), which can reduce occupancy — so measure, do
not assume. On CPU the same hoisting is what the compiler would do anyway if it
can prove the indices are loop-invariant.

**5. How to measure coalescing specifically.**

- `ncu --metrics l1tex__t_sectors_pipe_lsu_mem_global_op_ld.sum,l1tex__t_requests_pipe_lsu_mem_global_op_ld.sum` — sectors per request. 4 is ideal for 8-byte loads; 32 is fully scattered.
- Or the direct experiment: run the same kernel with an identity neighbour map
  (`iidx(jc,jb,n) = jc`) versus the real one. The difference is *entirely*
  coalescing, with all other effects held constant. That is the cleanest
  possible attribution and it needs no profiler.

---

## Q3. The derived type

**1. Why coalescing is broken.**

`t_cell_state` is 4 × 8 = **32 bytes**. `cells(jc,jk,jb)%rho` for consecutive
`jc` is 32 bytes apart, not 8. This is **array of structures (AoS)**: the
`rho` values a warp wants are interleaved with `theta_v`, `exner`, `w`, which
this kernel never touches.

**2. Fraction of each line used.**

8 bytes of every 32 → **25%**. A warp fetches 1024 B to use 256 B.

**3. The two layouts.**

- **AoS** — `cells(:)%rho`: one array of structs.
- **SoA** — `rho(:,:,:)`, `theta_v(:,:,:)`: separate arrays, usually grouped in a
  derived type of *arrays* (`p_prog%rho`, `p_prog%theta_v`).

ICON uses **SoA**: `TYPE t_nh_prog` holds `REAL(wp), POINTER :: rho(:,:,:)`,
`theta_v(:,:,:)`, etc. Each field is separately contiguous, so any kernel
touching one field gets full bandwidth utilisation, and each field can be
independently placed on the device with `ENTER DATA`.

**4. When AoS is better.**

When a kernel touches **all** members for the same cell — then the 32-byte
struct is exactly one useful unit and the line is fully used. A column physics
routine that reads all prognostic variables for a cell is the classic case. Also
when access is scattered anyway (a gather with poor locality), AoS gives you all
four values for one sector fetch instead of four.

The general rule: **layout should follow the access pattern of the dominant
kernel.** ESMs have far more single-field kernels than all-field ones, so SoA
wins overall.

---

# Part B — exploiting data parallelism

## Q4. The block loop

**1. What it launches.**

On the GPU build `nproma` = the whole patch, so `nblks = 1`. The only loop
carrying `GANG` is `jb`, with **one iteration**. So the kernel launches **1 gang
of 1 thread**, and the `jk`/`jc` loops run sequentially inside it.

On a device with 108 SMs and ~2048 threads per SM (~221,000 threads of
capacity), you are using **one**. Utilisation ≈ `1 / 221000` ≈ **0.0005%**.
In practice this will be hundreds to thousands of times slower than the CPU
version — a "GPU port" that is a catastrophic regression.

**2. The second, independent problem.**

`CALL get_indices_c(...)` is inside the `!$ACC PARALLEL` region. It is a host
routine and is not marked `!$ACC ROUTINE`, so it **cannot be called from device
code**. Depending on the compiler this is a compile error ("procedure not
available on device") or, if someone "fixed" it by slapping `!$ACC ROUTINE SEQ`
on it, it becomes a fully serial device call — correct but disastrous.

The right structure is: keep the `jb` loop and the index computation **on the
host**, and launch a kernel per block.

**3. Corrected.**

```fortran
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1)
  DO jk = 1, nlev
    DO jc = i_startidx, i_endidx
      p_diag%temp(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) * p_diag%exner(jc,jk,jb)
    END DO
  END DO
  !$ACC END PARALLEL
END DO
```

With `nblks = 1` this is one launch over `nlev × nproma` = 90 × 20000 = 1.8M
threads. That fills the device.

**4. Why the same loop is right on CPU and wrong on GPU.**

| | CPU | GPU |
|---|---|---|
| `nproma` | ~256 (cache-block sized) | = whole patch |
| `nblks` | thousands | 1 |
| parallelism needed | ~128 threads | ~10⁵–10⁶ threads |
| where it comes from | `jb` — thousands of blocks | `jc` — millions of cells |
| why | each block's working set fits L2; blocks are independent → no false sharing | only `jc` has enough width |

Same source, same loop nest, opposite parallel dimension — reconciled by making
`nproma` a namelist parameter and putting the directives on the right loops for
each. That is the whole design of ICON's GPU port in one table, and it is a
very good thing to be able to draw.

---

## Q5. The vertical solve

**1. What is wrong.**

`w(jc,jk,jb)` reads `w(jc,jk-1,jb)` — a **loop-carried dependence in `jk`**.
`COLLAPSE(2)` tells the compiler both loops are independent, which is false.
Threads computing level `jk` may read `w(jc,jk-1)` before the thread computing
level `jk-1` has written it, so they read the *old* value. The result is
garbage: instead of a cumulative sum down the column you get something close to
`w_old(jc,jk-1) + dz*src`.

**2. Deterministic?**

Formally no — it depends on warp scheduling. In practice it is often
*repeatably* wrong on a given device and problem size, because scheduling is
deterministic enough. That is worse than a flaky bug: it looks like a systematic
physics error, and it will change when you move to a different GPU or change
`nproma`, which is when someone finally notices.

**3. Correct directive.**

```fortran
!$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
DO jc = i_startidx, i_endidx
  !$ACC LOOP SEQ
  DO jk = 2, nlev
    w(jc,jk,jb) = w(jc,jk-1,jb) + dz(jc,jk,jb) * src(jc,jk,jb)
  END DO
END DO
!$ACC END PARALLEL
```

**`jc` carries the parallelism; `jk` must be sequential.** One thread per
column, walking down it.

Note the happy accident that is not an accident: because `jc` is the
fastest-varying dimension, consecutive threads at a given `jk` still access
consecutive memory, so **this is fully coalesced**. The `(nproma, nlev, nblks)`
layout was chosen precisely so that "parallel over cells, sequential over
levels" — the natural structure of atmospheric physics — is also the
coalescing-optimal structure. Say this out loud in an interview; it shows you
understand *why* the layout is what it is.

**4. Available parallelism.**

`i_endidx - i_startidx + 1` threads — the cells in this block. On a GPU build
with `nproma` = whole patch, that is O(10⁵–10⁶) per rank: plenty. It becomes
thin if you strong-scale to many ranks (few cells each) or on a small nested
domain. The check: threads ≥ ~4× (SMs × max warps per SM) to hide latency.

**5. Exposing more parallelism.**

A **parallel scan** (Blelloch prefix sum) over `jk`: depth `log₂(90) ≈ 7`
instead of 90, at the cost of ~2× the total work, more memory traffic, and much
more code. For a genuine tridiagonal solve, **parallel cyclic reduction** or
**thomas-per-thread with PCR fallback** is the standard answer.

Worth it only if `nlev` is large and cell-parallelism is exhausted — which for
90 levels and millions of cells it is not. The right answer here is "no, the
cell dimension already has more parallelism than the device can use." Knowing
when *not* to apply a clever algorithm is the more valuable judgement.

**6. Would `-Minfo=accel` have caught it?**

**No — and this is the important part.**

- `!$ACC KERNELS` asks the compiler to *analyse* the loop nest. It would detect
  the dependence and report `loop carried dependence of w prevents
  parallelization`, then run that loop sequentially.
- `!$ACC PARALLEL LOOP` is an **assertion by the programmer** that the
  iterations are independent. The compiler trusts you. `COLLAPSE(2)` extends
  that assertion to both loops.

So `parallel loop` + `collapse` silently accepts a wrong program.
`-Minfo=accel` will happily report `Generating Tesla code / 90, !$acc loop gang,
vector` and say nothing about the race.

> Interview line: "`kernels` is a question; `parallel loop` is a promise.
> `collapse` is a bigger promise. The compiler only checks the question."

---

# Part C — race conditions

## Q6. Edge to cell

**1. The defect.**

A **write–write / read-modify-write race** on `p_diag%div`. Multiple concurrent
edge iterations perform `div = div + flux` on the same cell element. The
read-modify-write is not atomic, so updates are lost.

**2. Why it is guaranteed, not merely possible.**

Every interior edge is shared by exactly **two** cells, and every triangular
cell has exactly **three** edges. So each cell's `div` is written by three
different edge iterations. In a vectorised loop over edges, iterations are
executed in warps of 32 and blocks of hundreds — the three edges of a cell are
essentially certain to be in flight simultaneously. This is not a rare
interleaving; it is the normal case.

**3. Three fixes.**

| Fix | Cost | Bit-reproducible? |
|---|---|---|
| `!$ACC ATOMIC UPDATE` on each accumulation | Serialisation under contention; atomics on doubles are slow on some architectures | **No** — FP addition order is non-deterministic |
| **Colour the edges** so no two edges in a colour share a cell; one kernel per colour | Multiple launches; less parallelism per launch; colouring must be computed and stored at setup | **Yes** — fixed order within a fixed colouring |
| **Invert the loop**: iterate over *cells*, gather from their 3 edges | Needs the flux stored in an edge array first (extra array + extra pass), or recompute | **Yes** — each output written by exactly one thread |

**4. Why ICON uses the cell-gather form despite the extra arithmetic.**

Four reasons, in order of importance:

1. **Determinism.** Each output element is written by exactly one thread, in a
   fixed order. That is a precondition for the bit-identical restart test
   ([Exercise 22](../22_regression_harness)). Atomics would forfeit it
   permanently.
2. **The kernel is memory-bound anyway.** Doubling flux flops costs almost
   nothing when you are at bandwidth. In practice ICON computes the edge flux
   once into an edge array, so there is no recompute at all — just one more
   array.
3. **Coalesced writes.** The cell loop writes `div(jc,jk,jb)` for consecutive
   `jc` — perfectly coalesced. The edge scatter writes to scattered addresses.
4. **No atomics, no colouring machinery** to maintain across grid changes.

> The general principle worth stating: **prefer gather to scatter.** Gather is
> race-free, deterministic, and has coalesced writes. Almost every scatter in a
> GPU code can be rewritten as a gather by inverting the connectivity at setup,
> and the connectivity inversion is a setup cost you pay once.

**5. On CPU / OpenMP.**

Yes, the same race exists if threads process edges that share cells. It is
**hidden** when:

- the run is single-threaded;
- the decomposition happens to give each thread a set of edges whose cells are
  disjoint (true in the interior of a block, false at block boundaries);
- there are few threads, so the collision window is small.

"It works on CPU" is not evidence of correctness because a race is a property of
the *program*, not of the run. Fewer threads means lower probability, not
absence. The right tools are a data-race detector (Intel Inspector, Archer/
ThreadSanitizer for OpenMP), not a test that passed.

---

## Q7. The OpenMP block loop

**1. The two defects.**

**(a) `i_startidx` and `i_endidx` are not `PRIVATE`.** They are assigned by
`get_indices_c` inside the parallel loop, so all threads write and read the same
two shared variables. Thread A can overwrite the bounds thread B is about to use.

**(b) `z_tmp` is a shared scratch array.** Every thread writes `z_tmp(jc,jk)`
for its own `jb` into the *same* array, then reads it back in the second loop.
Whatever the last writer put there is what everyone reads.

(`jb` in the `PRIVATE` list is harmless but redundant — the loop variable of an
`OMP DO` is private by definition.)

**2. Symptoms.**

- (a) Wrong loop bounds. With `-fcheck=bounds` you may get a clean abort; without
  it, silent out-of-range reads/writes, or simply computing the wrong subset of
  cells. Occasionally a segfault.
- (b) **Intermittent wrong answers.** No crash, no warning, plausible-looking
  fields, different results run to run.

**(b) is far worse.** A crash tells you where to look on the first run. A silent
non-deterministic numerical error in a climate model can survive for months,
gets attributed to "model variability", and poisons every result produced in
the meantime. It also breaks the run-to-run reproducibility test, which is often
how it is finally caught — which is the argument for having that test at all.

**3. Would a 2-thread or 128-thread test catch it?**

Neither reliably. 2 threads: the window is small and you may get lucky for
thousands of runs. 128 threads: much more likely, still not guaranteed, and now
you cannot tell *which* defect fired. Race detection needs a **tool**
(ThreadSanitizer/Archer, Intel Inspector), not more threads. This is worth
saying explicitly — "I would not try to test my way to confidence here" is a
strong answer.

**4. Corrected.**

```fortran
!$OMP PARALLEL DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_tmp) ICON_OMP_DEFAULT_SCHEDULE
```

...but privatising a large array via a clause allocates it per thread on the
stack every iteration. The idiomatic Fortran fix, and what ICON does, is to make
it an **automatic array local to the loop body**, so it is thread-local by
construction:

```fortran
!$OMP PARALLEL DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_tmp) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk
  BLOCK
    REAL(wp) :: z_tmp(nproma, nlev)     ! automatic -> per-thread
    ...
  END BLOCK
END DO
```

Watch the stack: `nproma × nlev × 8` per thread. At `nproma = 256`, `nlev = 90`
that is 184 kB — fine. At GPU-sized `nproma` it would blow the stack, which is
another reason the two builds use different `nproma`.

**5. `ICON_OMP_DEFAULT_SCHEDULE`.**

Work per block here is uniform — same `nlev`, and `nproma` cells in every block
except the last. So `SCHEDULE(STATIC)` is right: lowest overhead, best locality,
and each thread touches the same blocks every timestep (good for NUMA
first-touch).

What would change my mind, from [Exercise 02](../02_fortran_omp_reduction) and
[Exercise 16](../16_load_imbalance):

- **boundary blocks** have a smaller `i_endidx - i_startidx`, so the last block
  and halo blocks are cheaper — a small, fixed imbalance that `STATIC` handles
  badly only if `nblks` is close to the thread count;
- **data-dependent physics** in the loop body (convection firing, radiation on
  the day side) makes cost vary per block by 2–10× → `GUIDED` or
  `DYNAMIC` with a chunk large enough to amortise the scheduling overhead;
- if it is data-dependent *and* spatially clustered, the better answer is not a
  schedule at all but a **different decomposition** (Exercise 16 variant B).

---

## Q8. The reduction

**1. The three defects.**

**(a) `vcfl` is a shared scalar.** It is assigned inside the loop and read in
the next line. OpenACC will usually privatise a scalar assigned in a loop body,
but this is exactly the inference you should not rely on across compilers —
write `PRIVATE(vcfl)`.

**(b) No `REDUCTION(MAX:max_vcfl)`.** The `IF (vcfl > max_vcfl) max_vcfl = vcfl`
is a read-modify-write on a shared variable from every thread. Updates are lost;
the reported maximum is some arbitrary thread's value, essentially always an
under-estimate.

**(c) The result is read on the host with no `WAIT`.** The kernel is `ASYNC(1)`,
so `IF (max_vcfl > 1._wp)` executes on the host while the kernel may still be
running. Even with the reduction clause present, the value is not valid until
the async region completes.

**2. Corrected.**

```fortran
max_vcfl = 0._wp

!$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) &
!$ACC   PRIVATE(vcfl) REDUCTION(MAX:max_vcfl) ASYNC(1)
DO jk = 1, nlev
  DO jc = i_startidx, i_endidx
    vcfl = ABS(w(jc,jk,jb)) * dtime / dz(jc,jk,jb)
    max_vcfl = MAX(max_vcfl, vcfl)
  END DO
END DO
!$ACC END PARALLEL

!$ACC WAIT(1)                                   ! reduction result now valid
CALL MPI_Allreduce(MPI_IN_PLACE, max_vcfl, 1, MPI_DOUBLE_PRECISION, &
                   MPI_MAX, p_comm_work, ierr)  ! global, not just this rank
IF (max_vcfl > 1._wp) CALL finish('velocity_tendencies', 'CFL violated')
```

(Using `MAX(...)` rather than the `IF` also lets the compiler emit the intrinsic
reduction directly.)

**3. Why a race here is worse than in a printed diagnostic.**

This gates an **abort**. A lost update means the reported maximum is too small,
so a genuine CFL violation goes undetected. The model continues, goes
numerically unstable some number of steps later, and produces either a crash far
from the cause or — worse — plausible garbage. The safety check silently stops
being a safety check, and you only find out when someone questions a result.

A racy printed diagnostic gives you a wrong number on a log line. Bad, but
self-evidently a number, and it does not change the model's behaviour.

**4. Where the `MPI_Allreduce` goes, and what it costs.**

After `WAIT(1)`, on the host, over the work communicator — as above. The CFL
condition is global; a violation on any rank invalidates the timestep.

What it does to the pipeline: it is a **hard global synchronisation every
timestep**. It drains the async queue (you must `WAIT` before you have the
value) and then blocks on a collective. That is exactly the barrier you spent
Part E removing.

Mitigations, in increasing order of cleverness:

- check every N steps rather than every step (accepting that you detect a
  violation up to N steps late — usually fine, since the instability takes many
  steps to grow);
- `MPI_Iallreduce` immediately after the kernel, and `MPI_Test` it at the *end*
  of the timestep, so the collective overlaps the rest of the step;
- compute the check on the *previous* step's data, so it never sits on the
  critical path at all — the same one-window-lag argument as the ice sheet in
  [Exercise 23](../23_add_new_component).

**5. Reproducibility of `MAX` vs `SUM`.**

**`MAX` is safe.** It *selects* an element rather than combining values, so no
rounding occurs and the result is independent of evaluation order. A
non-deterministic reduction order gives a bit-identical answer.

**`SUM` is not.** Floating-point addition is not associative, so a different
combination order gives a different last few bits. This is the
[Exercise 02](../02_fortran_omp_reduction) TODO 4 problem and the
[Exercise 22](../22_regression_harness) rank-invariance problem: a global mass
sum will differ between rank counts and between GPU runs unless you use a
fixed-order or compensated reduction.

Two footnotes that show depth: `MAX` order-independence breaks with **NaN**
(propagation depends on order and on whether the intrinsic is `MAX` or a
comparison), and with **signed zeros**. In a CFL check on `ABS()` values neither
can arise — but say the caveat, because it demonstrates you know why the general
claim is not universal.

---

# Part D — data regions and host↔device transfers

## Q9. The timeloop

**1. The biggest problem.**

The `!$ACC DATA` region is **inside** the timestep loop. Every field named in
`COPYIN` is transferred host→device at the top of *every* timestep, and every
`COPYOUT` field is transferred back at the bottom. The arrays never stay
resident on the device.

**2. The arithmetic.**

One field: `5e6 cells × 90 levels × 8 B` = **3.6 GB**.

| Direction | Fields | Bytes |
|---|---|---|
| `COPYIN` | `rho`, `theta_v`, `ddqz_z_full` | 3 × 3.6 = 10.8 GB |
| `COPYOUT` | `div`, `temp` | 2 × 3.6 = 7.2 GB |
| **Total per timestep** | | **18 GB** |

At 25 GB/s effective (PCIe gen4 x16): **0.72 s per timestep of pure transfer.**

A real ICON dynamics step at this size is O(0.1 s) of compute. So transfers are
roughly **7× the compute** — the "GPU port" is a large net *slowdown*, and it
would look like the GPU is slow rather than like the data management is wrong.
This is the number that ends arguments.

**3. The rewrite.**

```fortran
! --- at model init, once ---
!$ACC ENTER DATA COPYIN(p_prog%rho, p_prog%theta_v, p_metrics%ddqz_z_full) &
!$ACC            CREATE(p_diag%div, p_diag%temp)

DO jstep = 1, nsteps
  CALL compute_divergence(p_patch, p_prog, p_diag)     ! DEFAULT(PRESENT) inside
  CALL compute_temperature(p_patch, p_prog, p_diag)
  CALL vertical_diffusion(p_patch, p_prog, p_diag, p_metrics)
  CALL nh_solve(p_patch, p_prog, p_diag)
END DO

! --- at finalize ---
!$ACC EXIT DATA COPYOUT(p_diag%div, p_diag%temp) DELETE(p_prog%rho, ...)
```

Use **`ENTER DATA` / `EXIT DATA`** (unstructured), not `DATA` (structured),
because the lifetime must span subroutine boundaries and the whole run — a
structured `DATA` region is lexically scoped and cannot.

Note `CREATE` rather than `COPYIN` for `div` and `temp`: they are outputs, and
their initial host contents are meaningless. Copying them in is pure waste, and
it is a very common oversight.

**4. `nh_solve` runs on the host.**

If it stays on the host it needs the data there, which forces `UPDATE HOST`
before and `UPDATE DEVICE` after — potentially reintroducing the whole problem.

Before assuming you can move it inside, check:

- **What fraction of runtime is it?** If it is 40% of the step, a partial port
  gives you Amdahl's ceiling regardless of how good the rest is.
- **Which fields does it actually touch?** Usually a subset, not the whole
  state. Update only those. Getting from "18 GB/step" to "one 3.6 GB field
  round-tripped" is already a 5× improvement and may be enough to be net
  positive while the port continues.
- **Is it a solver with a sequential structure** (tridiagonal, implicit)? Then
  it is Q5's problem and needs real work, not directives.

The honest engineering answer is usually: port it, but stage it — measure, port
the biggest piece, re-measure. A half-ported timeloop is often *slower* than a
CPU one, which is why GPU ports look bad in the middle and you have to warn
people in advance.

**5. Why per-module `ENTER DATA` rather than one big region.**

1. **Maintainability.** One region around the timeloop would have to name every
   array in a model with hundreds of modules, in one place, kept in sync with
   every change. It would be permanently wrong.
2. **Runtime-dependent field sets.** Which modules are active (which physics
   packages, which tracers) is a namelist decision. The set of live arrays is
   not known at compile time.
3. **Ownership matches allocation.** A module allocates its own arrays; it
   should manage its own device residency in the same place, so the two cannot
   drift apart.
4. **Memory.** Not everything needs to be resident simultaneously; per-module
   control lets you free what an inactive component does not need.

---

## Q10. The missing clause

**1. Which clause.**

`DEFAULT(PRESENT)` on the `!$ACC PARALLEL`.

Without it, OpenACC's implicit data handling for arrays referenced in a compute
construct is effectively **"present or copy"**: if the array is already on the
device, use it; **if it is not, silently allocate device memory and copy it in
(and out) around every kernel launch.**

**2. Does it cost anything at runtime if the data *is* present?**

**No.** If the data is present, the behaviour is identical — a present-table
lookup either way. The clause is free.

**Why ICON insists on it anyway:** the clause exists to change what happens in
the case where the data is **not** present. With `DEFAULT(PRESENT)`, a missing
`ENTER DATA` becomes a **runtime error** ("Present table lookup failed for ...")
that stops the run at the exact kernel and array. Without it, the same mistake
becomes a silent 10–100× slowdown that still produces correct results.

**3. The failure mode it converts.**

Someone adds a new prognostic field, wires it into a kernel, and forgets the
`ENTER DATA`. The model runs. The answers are right. It is just much slower —
and the slowdown is smeared across every kernel that touches the field, so no
single hotspot stands out.

It is hard to find because every normal signal says the code is fine: no error,
no warning, correct output, tests pass. You would only catch it by profiling
*and* by knowing what the timeline should look like. In a codebase where dozens
of people add fields, this happens constantly — hence the convention that every
single kernel carries `DEFAULT(PRESENT)`.

**4. What `nsys` shows.**

- **Data present:** one `cuLaunchKernel` per iteration. Clean timeline.
- **Data absent:** `cuMemAlloc` → `cuMemcpyHtoDAsync` → `cuLaunchKernel` →
  `cuMemcpyDtoHAsync` → `cuMemFree` **around every launch**. The timeline is a
  dense picket fence of tiny memcpys with slivers of kernel between them, and
  the memcpy row dominates. It is unmistakable once you have seen it once —
  which is a good reason to deliberately induce it once.

**5. The general principle.**

**Make the slow path a loud failure rather than a silent fallback.**

Equivalently: prefer an error to a degradation. A system that silently does the
wrong-but-correct thing will accumulate those cases faster than anyone can find
them. This applies well beyond OpenACC — it is the same argument as
`-Werror`, as failing a build with no `CMAKE_BUILD_TYPE`
([Exercise 21](../21_cmake_mixed_build)), and as making an unimplemented stub
abort rather than return zero.

---

## Q11. The diagnostic

**1. The two bugs.**

**Correctness:** `!$ACC UPDATE HOST(...) ASYNC(1)` is asynchronous, and the very
next line reads `p_prog%rho` on the host with no `!$ACC WAIT(1)`. The host races
the transfer.

**Performance:** the transfer happens **every step**, but the result is only used
every 100 steps. 99% of the traffic is discarded.

(Third, minor: `max_w` is computed every step and never used at all in the code
shown — dead work over a 3.6 GB array.)

**2. Mechanism and symptom.**

The `UPDATE` is queued on stream 1 and returns immediately. The host proceeds to
`SUM(p_prog%rho)` over the host array, which at that instant holds some mixture
of the previous step's values and however much of the current step's transfer
has landed.

If you printed `total_mass` every step you would see values that are
**plausible** — right order of magnitude, smoothly varying — but lagging by
roughly a step and jittering. Crucially, **they would differ between runs of the
same executable**, because the amount of transfer completed depends on timing.
That run-to-run difference is the tell, and it is why the reproducibility test
in [Exercise 22](../22_regression_harness) is worth having: it catches exactly
this class of bug, which no amount of eyeballing the output will.

**3. Both fixed.**

```fortran
DO jstep = 1, nsteps
  CALL dynamics_step(p_patch, p_prog, p_diag)

  IF (MOD(jstep, 100) == 0) THEN
    !$ACC UPDATE HOST(p_prog%rho, p_prog%w) ASYNC(1)
    !$ACC WAIT(1)
    total_mass = SUM(p_prog%rho) * cell_volume
    max_w      = MAXVAL(ABS(p_prog%w))
    WRITE(message_text,'(a,e18.10)') 'total mass: ', total_mass
    CALL message('dynamics', message_text)
  END IF
END DO
```

**Ideal frequency: exactly the diagnostic output frequency, and no more.** The
general rule — *transfer on demand, at the cadence of the consumer, never "just
in case"*.

**4. A design that transfers nothing.**

Do the reduction **on the device**:

```fortran
local_mass = 0._wp
!$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) &
!$ACC   REDUCTION(+:local_mass) ASYNC(1)
DO jb = 1, nblks
  DO jk = 1, nlev
    DO jc = 1, nproma
      local_mass = local_mass + p_prog%rho(jc,jk,jb)
    END DO
  END DO
END DO
!$ACC END PARALLEL
!$ACC WAIT(1)
CALL MPI_Allreduce(MPI_IN_PLACE, local_mass, 1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, p_comm_work, ierr)
total_mass = local_mass * cell_volume
```

This moves **8 bytes** instead of 3.6 GB — a factor of ~4.5 × 10⁸. It is also
faster in absolute terms even ignoring the transfer, because the GPU has far
more bandwidth to its own memory than the host has to the array.

**5. Reproducibility versus the device reduction.**

The two facts are reconciled by asking **what the number is used for**:

- **If it only feeds a printed diagnostic**, non-determinism in the last bits is
  acceptable. Document it, and make sure your regression test compares model
  *state*, not log lines — otherwise a harmless diagnostic will fail your CI
  and people will start ignoring CI.
- **If it feeds a decision** — a mass fixer, a conservation correction, an abort
  threshold — it is model state and must be deterministic. Then you need a
  fixed-order or compensated reduction (the `NCHUNK` approach from
  [Exercise 02](../02_fortran_omp_reduction) TODO 4), which costs perhaps 2–3×
  on the reduction but is still negligible against a 3.6 GB transfer.

> The reusable principle: **reproducibility requirements follow from the
> consumer, not from the quantity.** The same global sum can be allowed to
> wobble in one place and required to be exact in another.

---

# Part E — asynchronous execution and overlap

## Q12. The halo exchange

**1. Synchronisation points, and which are needed.**

Per timestep:

| # | Sync | Needed? |
|---|---|---|
| 1 | `!$ACC WAIT` after the `vt` kernel | **No** — the pack kernel is on the same queue and is ordered behind it automatically |
| 2 | `!$ACC WAIT` after the pack kernel | **Yes** (in the staged version) — the host must not read `sendbuf` before it is written |
| 3 | `!$ACC UPDATE HOST(sendbuf)` | implicit sync; **removable** with GPU-aware MPI |
| 4 | `MPI_Waitall` | **Yes**, but in the wrong place — see (3) below |
| 5 | `!$ACC UPDATE DEVICE(recvbuf)` | implicit sync; **removable** |
| 6 | `!$ACC WAIT` after the unpack kernel | **No** — subsequent kernels on the same queue are already ordered |

So of six, **one** is genuinely required, and even that disappears with
GPU-aware MPI.

**2. Cost of the waits.**

Each `!$ACC WAIT` drains the device and forces the host to round-trip to the
driver: typically **5–10 µs**, and it destroys the ability to overlap the
*launch* of the next kernel with the *execution* of the current one.

Order-of-magnitude for a real dynamics step: 30 exchanged fields × 3 waits ×
~10 µs ≈ **0.9 ms per timestep** of pure synchronisation. Whether that matters
depends entirely on the step time:

- step = 50 ms → ~2%. Annoying, not urgent.
- step = 5 ms (strong-scaled to many GPUs, which is the interesting regime) →
  **~18%**. Now it is the thing to fix.

Note the direction: **synchronisation overhead gets relatively worse as you
strong-scale**, exactly when you are trying to scale. That is why it is worth
fixing before you need it.

**3. The correctness bug in the MPI sequence.**

Look at the order:

```fortran
!$ACC UPDATE HOST(sendbuf)                              ! overwrites sendbuf
CALL MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)  ! waits on PREVIOUS step
CALL start_halo_exchange(sendbuf, recvbuf, req)         ! posts new Isends
```

Two problems, one fatal:

- **`sendbuf` is overwritten while the previous step's `MPI_Isend` from it may
  still be in flight.** A non-blocking send does not permit you to touch the
  buffer until the corresponding wait completes. Here the buffer is modified
  *before* `MPI_Waitall`. This is a classic buffer-reuse-before-completion bug:
  it corrupts the message the neighbour receives, and it is timing-dependent, so
  it appears as intermittent wrong halo values.
- **On the first iteration `req` is uninitialised**, so `MPI_Waitall` is
  undefined behaviour.

The fix is simply to wait *before* touching the buffer:

```fortran
CALL MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)  ! previous step done
!$ACC UPDATE HOST(sendbuf)                              ! now safe to overwrite
CALL start_halo_exchange(sendbuf, recvbuf, req)
```

...and initialise `req = MPI_REQUEST_NULL` before the loop, so the first
`Waitall` is a no-op.

> This is the kind of defect that is worth *leading with* in an interview. It is
> a correctness bug hiding inside a performance question, and noticing it
> signals that you read code rather than pattern-match on directives.

**4. GPU-aware MPI.**

```fortran
!$ACC HOST_DATA USE_DEVICE(sendbuf, recvbuf)
CALL start_halo_exchange(sendbuf, recvbuf, req)
!$ACC END HOST_DATA
```

Inside `HOST_DATA USE_DEVICE`, the Fortran array names resolve to **device**
addresses, so MPI receives device pointers and (with GPUDirect RDMA) moves bytes
GPU→NIC→GPU without touching host memory.

Requirements:

- MPI built with CUDA/ROCm support — verify, do not assume:
  `ompi_info --parsable --all | grep mpi_built_with_cuda_support`, or
  `ucx_info -d | grep -i cuda`;
- often `export OMPI_MCA_opal_cuda_support=1`;
- if the MPI is **not** GPU-aware, passing a device pointer typically
  **segfaults** rather than falling back gracefully. Test with a tiny message
  first.

This is [Exercise 19](../19_gpu_aware_mpi_halo) exactly.

**5. Restructured with overlap.**

```fortran
req = MPI_REQUEST_NULL
DO jstep = 1, nsteps

  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1)
  DO jk = 1, nlev
    DO jc = 1, nproma
      p_diag%vt(jc,jk,jb) = compute_tangential(jc,jk,jb)
    END DO
  END DO

  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
  DO i = 1, n_send
    sendbuf(i) = p_diag%vt(send_idx(i), send_lev(i), send_blk(i))
  END DO

  CALL MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)  ! previous step
  !$ACC WAIT(1)                        ! the ONE necessary sync: sendbuf ready

  !$ACC HOST_DATA USE_DEVICE(sendbuf, recvbuf)
  CALL start_halo_exchange(sendbuf, recvbuf, req)         ! non-blocking
  !$ACC END HOST_DATA

  ! ---- interior work, independent of the halo, on a DIFFERENT queue ----
  CALL update_interior(p_patch, p_diag)                   ! internally ASYNC(2)

  CALL MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
  req = MPI_REQUEST_NULL

  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
  DO i = 1, n_recv
    p_diag%vt(recv_idx(i), recv_lev(i), recv_blk(i)) = recvbuf(i)
  END DO

  !$ACC WAIT                            ! join queues 1 and 2 before next step
END DO
```

Six syncs became two, the host staging disappeared, and `update_interior` now
runs on queue 2 *while* the MPI transfer is in flight.

**6. The timeline.**

```
BEFORE
  device  [vt][W][pack][W]..........idle..........[unpack][W][interior]
  host                     [D2H][MPI              ][H2D]
                            ^--------- device idle ---------^

AFTER
  queue1  [vt][pack]                                 [unpack]
  queue2                 [========= interior ========]
  network            [========= MPI (device ptrs) ===]
                      ^-- device stays busy throughout --^
```

Practise drawing this. Together with the coupled timeline from
[Exercise 23](../23_add_new_component), it is the diagram most likely to earn
you the offer.

**7. The one measurement that proves overlap happened.**

The rigorous version is the `nsys` timeline showing `update_interior`'s kernels
executing *within* the MPI window. But the cheap, quantitative version needs no
profiler and is the better answer:

Measure three things — `t_before`, `t_after`, and `t_comm_alone` (the exchange
with no computation at all, as in
[Exercise 03](../03_hybrid_mpi_omp)) — then compute

```
overlap efficiency  η = (t_before - t_after) / t_comm_alone
```

η ≈ 1 means the communication is fully hidden. **η ≈ 0 means you moved the wait
but overlapped nothing**, which is the common outcome and the reason this
measurement matters.

The usual cause of η ≈ 0 is that **MPI makes no asynchronous progress**: a
non-blocking call may not actually transfer anything until you next call into
the MPI library. "Non-blocking" is not the same as "overlapped". Check for a
progress thread (`MPICH_ASYNC_PROGRESS`, UCX progress thread) or insert periodic
`MPI_Test` calls. This is Exercise 03 TODO 6(b), and it catches people out
constantly.

---

# Part F — the open question

## Q13. "What would you optimise?"

**1. The first three questions.**

1. *What is the run configuration?* Resolution, levels, node count, ranks per
   node, threads per rank, GPU or CPU partition. A kernel that matters at R2B4
   may be irrelevant at R2B9.
2. *What does the profile say?* I do not want to guess; I want the top five by
   inclusive time, and the MPI fraction. If there is no profile, that is my
   first task, not my first opinion.
3. *What is the goal?* Throughput (simulated years per wall day), or capability
   (fit a bigger problem), or cost (core-hours per simulated year)? They lead to
   different work — throughput may say "use fewer nodes more efficiently", cost
   may say the opposite.

**2. What to know about the machine.**

Memory bandwidth per node (measured, not vendor), cache sizes, cores/NUMA
domains per node, interconnect bandwidth *and* message rate, GPUs per node and
their interconnect (PCIe vs NVLink), whether MPI is GPU-aware, filesystem type
and striping. Roughly: enough to build the roofline and to know what a message
costs.

**3. The Amdahl arithmetic.**

Normalise total runtime to 1.

- **2× on a 5% kernel:** `0.95 + 0.05/2 = 0.975` → **2.5% faster**.
- **10% on a 60% kernel:** `0.40 + 0.60/1.10 = 0.40 + 0.545 = 0.945` → **5.5%
  faster**.

Take the second — it is more than twice the win. And note the ceiling: even
making the 5% kernel *infinitely* fast only buys 5%. Always compute the ceiling
before starting work; it is the cheapest way to avoid wasting a month.

**4. 35% of time in MPI — four causes and how to tell them apart.**

| Cause | Distinguishing measurement |
|---|---|
| **Load imbalance** (the time is *waiting*, not communicating) | Per-rank time in synchronising calls; the rank with the **least** collective time is the slow one ([Exercise 15](../15_pmpi_profiler) TODO 6). Confirm by timing compute alone per rank. |
| **Bandwidth-limited** | Message-size histogram is dominated by large messages; bytes/s approaches link bandwidth. |
| **Message-rate / latency-limited** | Histogram dominated by tiny messages; messages/s approaches the NIC limit (~10⁶–10⁷/node). Fix: aggregate fields into one exchange. |
| **Bad decomposition** (too much halo) | Halo bytes vs interior cells; compare partitioners ([Exercise 07](../07_partitioning)). Surface-to-volume ratio rising faster than expected with rank count. |

A fifth worth naming: **no overlap** — lots of time in `MPI_Wait` while the
device or CPU sits idle. Distinguished by comparing time in `Wait` against
`t_comm_alone`.

The important framing: **"35% in MPI" is not a diagnosis.** Most of that time is
usually not the network at all — it is imbalance showing up at the
synchronisation point. Saying that is a strong answer on its own.

**5. "It got slower on the new machine" — the first hour.**

1. **Establish the comparison is real.** Same resolution, same rank count, same
   number of steps? Get both logs side by side. Surprisingly often it is not the
   same run.
2. **Check the build.** Did `-O2/-O3` survive? Different compiler or version?
   A missing build type is a classic (see
   [Exercise 21](../21_cmake_mixed_build)). Diff the compile lines.
3. **Check placement.** Ranks per node, threads per rank, `OMP_PROC_BIND`,
   NUMA binding, GPU affinity. A wrong pinning routinely costs 2×.
4. **Split compute from communication.** Run the scaling harness
   ([Exercise 13](../13_scaling_harness)) at two node counts on both machines.
   If single-node compute time is the same and only multi-node differs, it is
   the network or the decomposition, not the code.
5. **Profile both** with the same tool ([Exercise 15](../15_pmpi_profiler)) and
   compare the *shape*, not just the total.
6. **Check I/O.** Different filesystem, different striping, different number of
   OSTs ([Exercise 20](../20_parallel_io)). I/O regressions on a new machine are
   extremely common and easy to miss because they do not show up in an MPI
   profile.

The meta-point: **spend the first hour narrowing, not fixing.** Almost everyone
starts optimising before they know what changed.

**6. A 3× win that changes the last three bits.**

- **Do not ship it silently.** Bit-identity is the model's regression contract
  and the basis of its entire test strategy
  ([Exercise 22](../22_regression_harness)). Breaking it quietly destroys the
  ability to distinguish a bug from a change, for everyone, forever.
- **Characterise the change.** Is it a rounding-order difference (reduction
  order, FMA contraction, a different transcendental library) or an algorithmic
  approximation? The first is usually acceptable; the second needs scientific
  review.
- **Produce evidence, not reassurance.** Run an ensemble of both versions with
  perturbed initial conditions and show that the new results are statistically
  indistinguishable from the old — the standard tools are a climate-statistics
  tolerance test or a KS/Wilcoxon test on the ensemble distributions. "It looks
  the same" is not evidence.
- **Who decides: not you.** The model's scientific owner or the working group.
  Your job is to bring the speedup, the mechanism, and the ensemble evidence,
  and to be clear about what is and is not established.
- **Practical resolution:** gate it behind a namelist switch, default off, so
  the reference configuration stays bit-reproducible and experiments can opt in.
  Then it can be validated over time and promoted once there is confidence.

> This question is really testing whether you understand that in climate
> modelling, **reproducibility is a scientific requirement, not an engineering
> preference** — and whether you know where your authority ends.

---

## After you have worked through these

Put the following in `notes/INTERVIEW.md`:

- **Three diagrams** you can draw in 30 seconds: Q12's before/after timeline,
  Q4's CPU-vs-GPU parallelism table, Q1's stride arithmetic.
- **Five one-liners:**
  - "`kernels` is a question, `parallel loop` is a promise."
  - "Prefer gather to scatter — race-free, deterministic, coalesced writes."
  - "Make the slow path a loud failure, not a silent fallback."
  - "Time in a collective is time spent waiting; the rank with the least of it is
    the slow one."
  - "Non-blocking is not the same as overlapped."
- **One number:** Q9's 18 GB/step → 0.72 s. It is the most vivid illustration of
  why data residency dominates everything else in a GPU port.
