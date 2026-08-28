# DKRZ interview — phone revision

Self-contained. No code to run, no files to open. Answers are collapsed —
**say your answer out loud before tapping.**

Target role: Research Software Engineer, HPC Code Performance Optimisation, DKRZ Hamburg.

Reported interview format: *shown a block of ICON code, asked what needs optimising.*
Themes: coalescing, data parallelism, races, data regions/transfers, async overlap.

**Sessions:** each Part is 5–15 min. Part 2 is the core — do it daily.

---

## Part 0 — the 60-second pitch

Have this ready. Rehearse it timed.

> I'm a computational scientist with an HPC background — CUDA, MPI/OpenMP,
> Fortran, PETSc. Most of my work has been performance: blocked matmul with a
> roofline analysis, a distributed Jacobi solver, an MPI+CUDA FFT, Cannon's
> algorithm on a 2D process grid.
>
> What drew me to this role is the coupling side. I've been working through
> conservative remapping between mismatched grids, multi-rate coupling where
> a slow component consumes a time-averaged flux, and asynchronous component
> coupling so an expensive infrequent model doesn't sit on the critical path.
>
> The gap I've been closing deliberately is unstructured grids and
> directive-based GPU — I knew CUDA, not OpenACC.

Adjust to be true. **Do not claim ICON experience you don't have.** "I read the
ICON grid documentation and built a small icosahedral connectivity code" is
strong and honest. "I've worked on ICON" is neither.

---

## Part 1 — read this idiom fluently

Every snippet you're shown will look like this.

```fortran
! f(nproma, nlev, nblks)
!  jc = cell in block   <- FASTEST-VARYING
!  jk = vertical level
!  jb = block

DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, ...,     &
       i_startidx, i_endidx, ...)

  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
  !$ACC   DEFAULT(PRESENT) ASYNC(1)
  DO jk = 1, nlev
    DO jc = i_startidx, i_endidx
      ...
    END DO
  END DO
END DO
```

Indirect neighbours:

```fortran
iidx => p_patch%cells%neighbor_idx  ! (nproma,nblks,3)
f(iidx(jc,jb,1), jk, iblk(jc,jb,1))
```

**The single fact behind half of all defects:**

- CPU: `nproma` ≈ 256, `nblks` large → thread over `jb`
- GPU: `nproma` = whole patch, `nblks` = 1 → parallel over `jc`

Same source. Opposite parallel dimension. Reconciled by `nproma` being a
namelist value.

---

## Part 2 — flashcards

### 2.1 Coalescing

**Q. Which index is fastest-varying, and what does `COLLAPSE(2)` do with it?**

<details><summary>Answer</summary>

`jc` is fastest-varying in `f(nproma,nlev,nblks)`.

`COLLAPSE(2)` linearises with the **inner** loop varying fastest. So the inner
loop must be `jc`, i.e. `jk` outer.

Address: `base + [(jc-1) + (jk-1)*nproma + (jb-1)*nproma*nlev] * 8`

</details>

**Q. `DO jc` outer, `DO jk` inner, `COLLAPSE(2)`. What's the stride between
consecutive threads at `nproma=20000`, real64?**

<details><summary>Answer</summary>

`nproma * 8` = **160,000 bytes**.

Warp of 32: touches 32 distinct sectors (1024 B fetched, 256 B used) vs the
ideal 8 sectors → **4× at 32-byte sector granularity**, **16×** at 128-byte
line granularity.

Say which granularity you mean. That's the answer that shows you know the
hardware rather than the slogan.

</details>

**Q. Is the wrong loop order also wrong on CPU?**

<details><summary>Answer</summary>

Yes — different mechanism. At `nproma=256` the inner `jk` loop strides 2048 B:
new cache line every iteration, uses 8 of 64 bytes, **and it's not stride-1 so
it won't vectorise**.

Key point: *the same fix helps both machines.* That's what makes one source
viable. The machine-specific knob is `nproma`, not the loop order.

</details>

**Q. Loop order is right, but you gather `flux(iidx(jc,jb,1), jk, ...)`.
Does it coalesce?**

<details><summary>Answer</summary>

Depends on the **cell numbering**, decided at setup.

The *index* loads (`iidx`, `iblk`, `div_coeff`) all coalesce — consecutive `jc`.
The gathered `flux` address depends on the *value* of `iidx`. Space-filling-curve
numbering → neighbours close in index space → near-contiguous. Random numbering
→ every gather a separate sector.

9 of the 12 loads per output element are perfect; the 3 that matter are the gather.

</details>

**Q. Real optimisation in that gather kernel?**

<details><summary>Answer</summary>

`iidx`/`iblk`/`div_coeff` **don't depend on `jk`**. With `COLLAPSE(2)` you
re-read them for all 90 levels.

Restructure: `GANG VECTOR` over `jc`, hoist the 9 index/coeff loads into
registers, `LOOP SEQ` over `jk`. Cuts index traffic 90×, stays coalesced.

Cost: register pressure → possibly lower occupancy. Measure.

</details>

**Q. AoS vs SoA — a 4-double struct, kernel touches one member. Efficiency?**

<details><summary>Answer</summary>

**25%** — 8 useful bytes of every 32.

ICON uses **SoA**: `p_prog%rho(:,:,:)`, `p_prog%theta_v(:,:,:)` — a derived type
*of arrays*, not an array of structs.

AoS wins only when a kernel touches **all** members for the same cell (a column
physics routine), or when access is scattered anyway.

Rule: layout follows the access pattern of the dominant kernel. ESMs have far
more single-field kernels.

</details>

### 2.2 Data parallelism

**Q. `!$ACC PARALLEL LOOP GANG` over `jb`, GPU build. What launches?**

<details><summary>Answer</summary>

`nblks = 1` → **one gang, one thread**. On 108 SMs (~221,000 threads capacity)
that's ~0.0005% utilisation. Hundreds to thousands of times slower than the CPU.

Second, independent bug: `get_indices_c` is a **host** routine called inside a
`PARALLEL` region. Not valid on device without `!$ACC ROUTINE` — and marking it
`ROUTINE SEQ` would make it correct but serial.

Fix: keep `jb` and the index call on the host, launch a kernel per block over
`jk`/`jc`.

</details>

**Q. Why is `jb` the right parallel dimension on CPU and the wrong one on GPU?**

<details><summary>Answer</summary>

- CPU needs ~128-way parallelism; `nblks` is thousands. Each block's working set
  fits L2. Blocks are independent → no false sharing.
- GPU needs 10⁵–10⁶-way; only `jc` is that wide, so `nproma` is set to the whole
  patch and `nblks` collapses to 1.

Draw this as a two-column table. It's the design of ICON's GPU port in one image.

</details>

### 2.3 Races

**Q. `w(jc,jk,jb) = w(jc,jk-1,jb) + ...` under `COLLAPSE(2)`. What happens?**

<details><summary>Answer</summary>

Loop-carried dependence in `jk`. `COLLAPSE(2)` **asserts** independence — it's
false. Threads read `w(jc,jk-1)` before it's written → garbage.

Fix:
```fortran
!$ACC PARALLEL LOOP GANG VECTOR ...
DO jc = i_startidx, i_endidx
  !$ACC LOOP SEQ
  DO jk = 2, nlev
    w(jc,jk,jb) = w(jc,jk-1,jb) + ...
```

**`jc` carries the parallelism. `jk` must be sequential.**

</details>

**Q. Why is that fix also the coalescing-optimal form?**

<details><summary>Answer</summary>

Because `jc` is the fastest-varying dimension — consecutive threads at a given
`jk` still access consecutive memory.

**The `(nproma,nlev,nblks)` layout was chosen so that "parallel over cells,
sequential over levels" — the natural structure of atmospheric physics — is
simultaneously the coalescing-optimal structure.**

Say this out loud. It shows you understand *why* the layout exists.

</details>

**Q. Would `-Minfo=accel` have caught the dependence?**

<details><summary>Answer</summary>

**No.**

- `!$ACC KERNELS` = you ask the compiler to analyse. It *would* report
  "loop carried dependence prevents parallelization".
- `!$ACC PARALLEL LOOP` = **you assert** independence. `COLLAPSE` extends the
  assertion. The compiler trusts you and says nothing.

> "`kernels` is a question, `parallel loop` is a promise. The compiler only
> checks the question."

</details>

**Q. Loop over edges, `div(cell) = div(cell) + flux`. Why is the race
*guaranteed*, not just possible?**

<details><summary>Answer</summary>

Every interior edge is shared by exactly **2** cells; every triangle has **3**
edges. So each cell's `div` is written by 3 different edge iterations — and in a
warp of 32 they're certainly in flight together.

Not a rare interleaving. The normal case.

</details>

**Q. Three fixes for that scatter, and their costs?**

<details><summary>Answer</summary>

1. **`!$ACC ATOMIC UPDATE`** — serialises under contention. **Not
   bit-reproducible** (FP addition order).
2. **Colour the edges**, one kernel per colour — more launches, less parallelism
   each, setup cost. **Reproducible.**
3. **Invert to a cell-gather** — loop over cells, gather from their 3 edges.
   **Reproducible**, coalesced writes, no atomics. **ICON's choice.**

Why gather wins despite more arithmetic: the kernel is memory-bound anyway, so
extra flops are ~free; and determinism is a precondition for the bit-identical
restart test.

> **Prefer gather to scatter.** Race-free, deterministic, coalesced writes.

</details>

**Q. `!$OMP PARALLEL DO PRIVATE(jb,jc,jk)` over blocks, with
`get_indices_c` inside and a shared `z_tmp` scratch array. What's wrong?**

<details><summary>Answer</summary>

Two defects:

1. **`i_startidx`/`i_endidx` not `PRIVATE`** — assigned inside the loop, shared
   by all threads → wrong loop bounds.
2. **`z_tmp` is shared** — every thread writes it for its own `jb`, then reads
   it back → intermittent wrong answers.

(2) is worse: no crash, plausible fields, different every run.

Fix: add both to `PRIVATE`; better, declare `z_tmp` as an **automatic array
inside the loop body** so it's thread-local by construction. Watch the stack.

`jb` in `PRIVATE` is redundant — an `OMP DO` loop variable is private by
definition.

</details>

**Q. Would more threads reliably expose a race?**

<details><summary>Answer</summary>

No. More threads = higher probability, not certainty. A race is a property of
the *program*, not the run.

Use a detector: ThreadSanitizer/Archer, Intel Inspector. "I wouldn't try to test
my way to confidence here" is a strong answer.

</details>

**Q. A CFL check accumulating `max_vcfl` in an `ASYNC(1)` kernel, read on the
host straight after. Three defects?**

<details><summary>Answer</summary>

1. `vcfl` scalar not explicitly `PRIVATE` (usually inferred, don't rely on it)
2. **No `REDUCTION(MAX:max_vcfl)`** → lost updates → under-estimate
3. **No `!$ACC WAIT(1)`** before the host reads it → invalid value

Plus: it needs an `MPI_Allreduce(MPI_MAX)` — CFL is global.

Why this race is worse than in a diagnostic: it gates an **abort**. An
under-reported max means a real CFL violation goes undetected, the run goes
unstable later, and the failure appears far from the cause.

</details>

**Q. Is a GPU `MAX` reduction bit-reproducible? A `SUM`?**

<details><summary>Answer</summary>

- **`MAX` — yes.** It *selects* an element, no rounding, order-independent.
- **`SUM` — no.** FP addition isn't associative; different combination order →
  different last bits.

So a global max is safe with any reduction order; a global mass sum needs a
fixed-order or compensated reduction to survive a rank-count change.

Depth footnote: `MAX` order-independence breaks with **NaN** and **signed
zeros**. Irrelevant for `ABS()` values — but saying the caveat shows you know
why the general claim isn't universal.

</details>

### 2.4 Data regions

**Q. `!$ACC DATA COPYIN(...)` *inside* the timestep loop. Quantify.**

<details><summary>Answer</summary>

One field, 5M cells × 90 levels × 8 B = **3.6 GB**.

3 `COPYIN` + 2 `COPYOUT` = **18 GB per timestep**.
At 25 GB/s → **0.72 s/step of pure transfer**.

A real step is ~0.1 s of compute. **Transfers are ~7× the compute** — the "GPU
port" is a large net slowdown.

This number ends arguments. Memorise it.

</details>

**Q. Correct data management?**

<details><summary>Answer</summary>

`!$ACC ENTER DATA` at init, `EXIT DATA` at finalize. **Unstructured**
(`ENTER`/`EXIT`), not structured `DATA`, because the lifetime spans subroutine
boundaries.

Use **`CREATE`** not `COPYIN` for outputs — their initial host contents are
meaningless. Very common oversight.

Why ICON does it **per module** rather than one big region: hundreds of modules;
which arrays are live is a runtime (namelist) property; ownership should match
allocation; and not everything needs to be resident at once.

</details>

**Q. `DEFAULT(PRESENT)` — what does it cost, and why insist on it?**

<details><summary>Answer</summary>

**Costs nothing** if the data is present — same present-table lookup.

Its purpose is the case where data is **absent**. Without it: OpenACC silently
allocates and copies in/out **around every kernel launch**. Correct answers,
10–100× slower, smeared across every kernel — nearly unfindable.

With it: a **runtime error** naming the exact kernel and array.

> **Make the slow path a loud failure, not a silent fallback.**

In `nsys`: data present = clean kernel launches. Data absent = a picket fence of
tiny `cuMemcpy` calls around every launch. Unmistakable once seen.

</details>

### 2.5 Async and overlap

**Q. `!$ACC UPDATE HOST(rho) ASYNC(1)` then `SUM(rho)` on the host. Symptom?**

<details><summary>Answer</summary>

Host races the transfer. Reads a mix of previous-step and current-step data.

Symptom: **plausible** numbers — right magnitude, smoothly varying, lagging ~a
step and jittering. Crucially they **differ between runs**. That run-to-run
difference is the tell, and it's why a reproducibility test earns its keep.

Fix: `!$ACC WAIT(1)` before the host read — and only transfer at the cadence of
the consumer (every 100 steps, not every step).

</details>

**Q. Better than transferring the field to sum it on the host?**

<details><summary>Answer</summary>

Do the reduction **on the device** with `REDUCTION(+:local_mass)`, then
`MPI_Allreduce` the scalar.

**8 bytes instead of 3.6 GB** — a factor of ~4.5×10⁸. Also faster in absolute
terms: the GPU has far more bandwidth to its own memory.

Reproducibility caveat: the device reduction order is non-deterministic. Fine if
it only feeds a printed diagnostic; needs a fixed-order reduction if it feeds a
*decision* (mass fixer, abort threshold).

> Reproducibility requirements follow from the **consumer**, not the quantity.

</details>

**Q. GPU-aware MPI — the directive, and the prerequisite?**

<details><summary>Answer</summary>

```fortran
!$ACC HOST_DATA USE_DEVICE(sendbuf, recvbuf)
CALL start_halo_exchange(sendbuf, recvbuf, req)
!$ACC END HOST_DATA
```

Inside that block the array names resolve to **device** addresses.

Prerequisite: MPI built with CUDA/ROCm support. **Verify, don't assume:**
`ompi_info --parsable --all | grep cuda_support`, often plus
`export OMPI_MCA_opal_cuda_support=1`.

If it isn't GPU-aware it typically **segfaults** rather than falling back.

</details>

**Q. `UPDATE HOST(sendbuf)` then `MPI_Waitall` on the previous step's requests.
Bug?**

<details><summary>Answer</summary>

**`sendbuf` is overwritten while the previous `MPI_Isend` from it may still be in
flight.** Classic buffer-reuse-before-completion. Corrupts the neighbour's data,
timing-dependent, appears as intermittent wrong halos.

Also: on the first iteration `req` is uninitialised → undefined behaviour.

Fix: `Waitall` **before** touching the buffer; initialise `req =
MPI_REQUEST_NULL`.

> Worth *leading with* — a correctness bug hiding inside a performance question.
> Noticing it says you read code rather than pattern-match directives.

</details>

**Q. How do you prove you actually overlapped anything?**

<details><summary>Answer</summary>

Measure `t_before`, `t_after`, and `t_comm_alone` (the exchange with no compute).

```
η = (t_before − t_after) / t_comm_alone
```

η ≈ 1 → fully hidden. **η ≈ 0 → you moved the wait but overlapped nothing** —
the common outcome.

Usual cause: **MPI makes no asynchronous progress**. A non-blocking call may not
transfer until you next call into MPI.

> **Non-blocking is not the same as overlapped.**

</details>

### 2.6 Coupling (what the posting actually describes)

**Q. Conservation vs consistency — one line each.**

<details><summary>Answer</summary>

- **Conservation:** `Σ(F_tgt·area_tgt) = Σ(f_src·area_src)`. Needed for
  **fluxes** — heat, freshwater, momentum.
- **Consistency:** a constant field remaps to the same constant (weights per
  target cell sum to 1). Needed for **state** — temperature, salinity.

First-order conservative remapping gives both. Most interpolation gives only the
second.

Why it matters: non-conservative coupling makes the system gain or lose energy
every coupling step. 48 steps/day over 100 years is a **fake climate trend**, not
a rounding error.

</details>

**Q. Ice sheet needs the surface mass balance every 100 atmosphere steps. Why
can't you just sample it?**

<details><summary>Answer</summary>

You'd miss every melt event between samples, and the error changes sign with
phase — worse than a consistent bias, because it looks like variability.

**Accumulate `flux·dt` over the window and hand over the time-average.** `dt`-weighted, not a plain mean — real models vary the timestep, and they shorten
it precisely when the flux is largest.

</details>

**Q. The classic restart bug in that scheme?**

<details><summary>Answer</summary>

**The accumulator is model state and gets left out of the restart file.**

It isn't "physics state", so it's forgotten. A restarted run then differs from a
continuous one by up to one coupling window of flux. Small, plausible, and only
appears at restart points that aren't multiples of the window — which is why it
survives for years.

</details>

**Q. Ice sheet costs 8× a step but runs 100× less often. Naive coupling?**

<details><summary>Answer</summary>

Every 100 steps **everything stalls** while it runs. Average looks like 8%, but
the whole machine is idle for it, and throughput drops more once you add the
synchronisation.

Fix: **concurrent + lagged.** Give it its own ranks; it receives the accumulated
forcing non-blockingly, computes while atmosphere and ocean carry on, and
delivers at the *next* window.

Justification for the lag: ice responds over centuries; one coupling window is
far below its timescale. **State the physical argument** — that's what makes it
legitimate engineering rather than a corner cut.

Failure mode to detect: if it never finishes inside a window, the lag grows
without bound. Report it loudly, don't let it drift.

</details>

---

## Part 3 — numbers to quote

- **3.6 GB** — one ICON field, 5M cells × 90 levels × real64
- **18 GB/step → 0.72 s at 25 GB/s** — a data region inside the timeloop
- **160,000 B** — thread stride with the loop order inverted at `nproma=20000`
- **4× / 16×** — coalescing amplification at sector / line granularity
- **0.0005%** — GPU utilisation when you gang over `jb` with `nblks=1`
- **25%** — cache-line efficiency of AoS with a 4-double struct
- **8 bytes vs 3.6 GB** — device-side reduction vs transferring the field
- **~5–10 µs** — cost of one `!$ACC WAIT` (device drain)
- **Amdahl:** 2× on 5% → 2.5% total. 10% on 60% → **5.5%**. Take the second.

---

## Part 4 — one-liners

Memorise. These are what get remembered after you leave the room.

1. "`kernels` is a question, `parallel loop` is a promise. The compiler only
   checks the question."
2. "Prefer gather to scatter — race-free, deterministic, coalesced writes."
3. "Make the slow path a loud failure, not a silent fallback."
4. "Time in a collective is time spent **waiting**. The rank with the *least* of
   it is the slow one."
5. "Non-blocking is not the same as overlapped."
6. "The layout is chosen so the natural physics loop order is also the
   coalescing-optimal one."
7. "Reproducibility requirements follow from the consumer, not the quantity."

---

## Part 5 — "what would you optimise?"

No snippet. Have a **procedure**, not an answer.

**First three questions back:**
1. What's the run configuration? Resolution, nodes, ranks × threads, CPU or GPU.
2. What does the profile say? Top five by inclusive time, and the MPI fraction.
3. What's the goal — throughput, capability, or cost? They lead to different work.

**"35% of time is in MPI." Four causes, four measurements:**

- **Load imbalance** — per-rank time in synchronising calls; *least* collective
  time = slowest rank
- **Bandwidth-limited** — histogram dominated by large messages; bytes/s near
  link limit
- **Message-rate-limited** — histogram dominated by tiny messages; msgs/s near
  NIC limit (~10⁶–10⁷/node). Fix: aggregate fields into one exchange
- **Bad decomposition** — halo bytes vs interior; compare partitioners

> "35% in MPI" is **not a diagnosis**. Most of it is usually imbalance showing up
> at the synchronisation point, not the network. Saying that is a strong answer
> on its own.

**"It got slower on the new machine." First hour:**
1. Is the comparison even the same run? Get both logs.
2. Did the optimisation flags survive the rebuild?
3. Pinning — ranks/node, `OMP_PROC_BIND`, NUMA, GPU affinity. Routinely 2×.
4. Single-node vs multi-node: separates compute from network.
5. Profile both, compare the *shape*.
6. Check I/O — different filesystem/striping. Commonly missed.

> Spend the first hour **narrowing, not fixing**.

**"3× faster but the last three bits change."**
- Don't ship it silently — bit-identity is the model's regression contract.
- Characterise: rounding-order (usually fine) vs algorithmic approximation
  (needs review).
- Evidence, not reassurance: ensemble comparison, statistical test.
- **The scientific owner decides, not you.** You bring speedup + mechanism +
  evidence.
- Practical: namelist flag, default off, reference config stays reproducible.

---

## Part 6 — ask them these

Have four. Pick from:

- Which model(s) would I be working on, and what's the current biggest
  performance pain point?
- Is the ice-sheet coupling work YAC-based — new development or extending
  existing infrastructure?
- How is the role split between GPU porting, CPU optimisation, and coupling
  infrastructure?
- Which European systems is the code being adapted for right now? Is the
  OpenACC/OpenMP-target question live, given LUMI is AMD?
- How do you handle the bit-reproducibility question when an optimisation
  changes results?

That fourth one is strong — it shows you've thought about the portability
problem their own job ad implies.

---

## Part 7 — be honest about

Prepare one sentence for each. "I hadn't worked with X, so I built Y" beats a
bluff every time.

- No production ICON experience — but you understand the grid, the `nproma`
  layout, and the coupling model
- CUDA before OpenACC — you know the hardware; the directive vocabulary is
  recent
- Coupling is study, not shipped work
- No German (if true) — the ad asks for English, so this is fine; say you're
  willing to learn

---

## Deep-dive index

Full material is in this repo:

- `PLAN.md` — the 6-day plan and gap analysis
- `25_icon_code_review/QUESTIONS.md` — 13 full code snippets
- `25_icon_code_review/ANSWERS.md` — worked answers with the arithmetic
- `17_openmp_target/` `18_openacc_port/` `19_gpu_aware_mpi_halo/` — GPU
- `08_nproma_blocking/` — the layout, measured
- `10_conservative_remap/` `12_coupling_timestep_lag/` `23_add_new_component/` — coupling
