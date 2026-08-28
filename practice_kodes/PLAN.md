# 6-Day Prep Plan — DKRZ Research Software Engineer, HPC Code Performance Optimisation

**Target role:** [Research Software Engineer (all genders) HPC Code Performance Optimisation](https://dkrz.softgarden.io/job/66326173?l=en) — Deutsches Klimarechenzentrum GmbH (DKRZ), Hamburg.
**Window:** Tue 25 Aug → Sun 30 Aug 2026 (6 days).

---

## ⚠️ Read this before Day 1

The posting says: *"Open until filled. Priority will be given to applications received by **23 August 2026**."*

That date was **two days ago**. Do not finish this plan and then apply — **send the application today (Tue 25 Aug)**, then keep working through these six days so you are sharp for the technical interview. The role is funded to end of 2028 with a September 2026 target start, which means they are screening *now*.

What "today" means concretely:
- CV + cover letter out the door before you write a single line of Exercise 01.
- Cover letter should name the three things they actually asked for: **MPI/OpenMP at scale**, **C/C++ and/or Fortran on Linux**, **scripting** — and lead with the GPU experience, which is on their "advantageous" list and is the thing most applicants for an Earth-system RSE job will not have.
- Link the repos. `cuda_stuff/jacobian_solver`, `low_level_optimisations/`, and `open_mpi_openmp_stuff/cannon_mat_mult` are already stronger evidence than anything you'll build this week.

This week's exercises are for the **interview**, not the application.

---

## ⚠️ What the interview actually looks like

Second-hand from someone who sat a similar DKRZ interview about a year ago: **they were shown a block of ICON code and asked what needed optimising.** The themes that came up were

- memory coalescing and exploiting data parallelism
- race conditions
- where to place data regions and host↔device transfers
- overlapping asynchronous operations

That is a **code-reading interview about GPU/OpenACC**, not a whiteboard-algorithms interview. It reorders this week:

| Priority | Material | Why |
|---|---|---|
| **1** | **[Exercise 25](25_icon_code_review/QUESTIONS.md)** — 13 ICON-style snippets with planted defects | This *is* the interview format |
| **2** | Day 5 (17, 18, 19) — directive-based GPU | The vocabulary those questions are asked in |
| **3** | Day 2 (08) — `nproma` layout, indirect addressing | Why coalescing works the way it does in ICON |
| **4** | Day 3 (09, 10, 12) + Exercise 23 | Still the job description; likely the second half of the interview |

Exercises 01–04, 13–16 and 21–22 remain worth doing, but if the week compresses, they are what gets cut. Note this is one data point about one interview a year ago — treat it as a strong hint, not a syllabus. The coupling material stays on the list because it is what the *posting* describes.

## What the job actually is

Stripping the HR language, the posting describes four concrete duties:

| Posting says | In practice this means |
|---|---|
| "collaborate with Earth System scientists to overcome technical challenges in HPC and software development" | You are the person a climate scientist brings a 500k-line Fortran model to when it runs 3× slower than it should. |
| "further developing, testing and optimising existing coupled Earth System models" | ICON / ICON-ESM. Fortran 2008, MPI + OpenMP, OpenACC on GPU partitions, unstructured icosahedral grid. |
| "extending models with new components (e.g. ice sheet models) without compromising the efficient execution time" | Bolting a new component into a **coupler** (DKRZ's is [YAC](https://dkrz-sw.gitlab-pages.dkrz.de/yac/)) and keeping the coupled timeline from collapsing — the hardest and most interesting part of the job. |
| "adapt complex models for European HPC systems" | Porting/tuning across Levante, LUMI, Leonardo, MareNostrum, JUPITER — different vendors, different compilers, different interconnects, different GPUs. |

**Required:** HPC + MPI/OpenMP, C/C++ and/or Fortran under UNIX/Linux, scripting, strong English.
**Advantageous:** GPU for scientific applications; algorithms, data structures, scientific software design.

### Your gap analysis

Your existing repo already covers a lot of this. Being honest about where you stand:

**Already strong — do not spend the week here:**
- MPI point-to-point, collectives, one-sided, communicators (`cannon_mat_mult`, `mpi_timer`, `job_scheduler`)
- CUDA and GPU kernel optimisation (`cuda_stuff/`, `low_level_optimisations/fft`)
- Roofline, blocking, cache optimisation (`low_level_optimisations/blocked_matrix_multiplication`)
- Parallel I/O awareness (`file_io_stuff/`)
- PETSc, containers, SLURM

**The real gaps for *this specific job*:**
1. **Fortran as a parallel language.** Your Fortran folder is serial algorithms and data structures. ICON is Fortran 2008 + MPI + OpenMP + OpenACC. You have never written `use mpi_f08` in anger. → **Day 1**
2. **Unstructured grids.** Every parallel thing you've built is on a structured Cartesian mesh with regular neighbours. ICON's icosahedral grid has indirect addressing, irregular halos, and `nproma` blocking. This is the single biggest technical gap. → **Day 2**
3. **Coupling.** You have never split a communicator into two independent models that exchange interpolated fields on mismatched grids at mismatched timesteps. This is literally the job description. → **Day 3**
4. **Directive-based GPU (OpenACC / OpenMP target).** You know CUDA. ESMs do not use CUDA — they use directives, because a climate scientist has to be able to read the loop. → **Day 5**
5. **Being the person who *diagnoses*, not just writes.** The job is optimising *someone else's* model. Profiling workflow, PMPI interception, load-imbalance forensics. → **Day 4**

---

## How to work through this

Each day is roughly **6–8 hours**, split into three blocks:

- **Block A** (~3 h) — exercises, the two harder ones while you're fresh
- **Block B** (~3 h) — exercises, then benchmark and write down numbers
- **Block C** (~1–1.5 h) — reading + a short written summary in `notes/dayN.md`

**Non-negotiable rule: every exercise ends with a number.** Not "it works" — a runtime, a speedup, a parallel efficiency, a bandwidth. DKRZ is hiring a *performance* engineer. In the interview you want to say "I got the halo exchange from 34% of runtime down to 11% by overlapping it with the interior update" — not "I implemented a halo exchange."

Keep a running `notes/` directory. On Day 6 it becomes your interview cheat-sheet.

### Getting started

```bash
cd /Users/prabhsharan/Desktop/mhpc/prabhkodes/open_mpi_openmp_stuff/practice_kodes && cat README.md
```

Every exercise directory has a source file whose header comment contains the goal, the DKRZ relevance, the tasks, and the acceptance criteria. Build and run with:

```bash
make run NP=4 OMP=2
```

---

## Day 1 — Tue 25 Aug · Fortran as a parallel language

*Rationale: ICON is Fortran. If you can't write `mpi_f08` fluently, nothing else this week matters.*

| # | Exercise | Time | Why it's here |
|---|---|---|---|
| 00 | `ring_allreduce.c` *(already present — warm-up)* | 30 m | Finish the TODOs you left. Then answer: why is ring allreduce bandwidth-optimal but latency-bad? |
| 01 | `01_fortran_mpi_halo_1d` | 90 m | `use mpi_f08`, derived types, `MPI_Isend/Irecv`, `MPI_Waitall`. The bread-and-butter pattern of every ESM. |
| 02 | `02_fortran_omp_reduction` | 60 m | OpenMP in Fortran: `collapse`, custom reductions, false sharing, `schedule` choice. |
| 03 | `03_hybrid_mpi_omp` | 90 m | `MPI_Init_thread`, FUNNELED vs MULTIPLE, overlapping halo comms with interior compute. |
| 04 | `04_c_fortran_interop` | 75 m | `iso_c_binding`, `bind(C)`, array descriptors, column- vs row-major. **YAC's API is C, ICON is Fortran — this boundary is a real part of the job.** |

**Block C reading:** Fortran 2008 `mpi_f08` vs the old `mpi` module (why `type(MPI_Comm)` instead of `integer` matters for type safety). Skim the [YAC documentation](https://dkrz-sw.gitlab-pages.dkrz.de/yac/) front page.

**End-of-day deliverable:** `notes/day1.md` — a table of halo-exchange time vs message size from Ex. 01, and the blocking-vs-nonblocking crossover point.

---

## Day 2 — Wed 26 Aug · Unstructured grids & the ICON data layout

*Rationale: this is your biggest gap. Everything you've parallelised has had regular neighbours. ICON does not.*

| # | Exercise | Time | Why it's here |
|---|---|---|---|
| 05 | `05_icosahedral_grid` | 90 m | Build cell/edge/vertex connectivity from an icosahedron, refine it, verify Euler's formula. You must be able to talk about this grid. |
| 06 | `06_unstructured_halo` | 120 m | Halo exchange with **no** regular neighbour pattern: build send/recv index lists, pack/unpack buffers, `MPI_Neighbor_alltoallv`. |
| 07 | `07_partitioning` | 90 m | Partition the mesh (space-filling curve vs naive). Measure **edge cut** and **load imbalance**. This is where coupled-model performance is won or lost. |
| 08 | `08_nproma_blocking` | 90 m | ICON's actual data layout: `field(nproma, nlev, nblks)`. Sweep `nproma`, find the cache/vector sweet spot, explain the curve. |

**Block C reading:** Why ICON uses `nproma` blocking at all (vectorisation + cache blocking + a single layout that works on CPU *and* GPU with different `nproma`). Look at the ICON grid description.

**End-of-day deliverable:** `notes/day2.md` — your `nproma` sweep plot, and the edge-cut comparison between the two partitioners.

---

## Day 3 — Thu 27 Aug · Coupling — the actual job

*Rationale: "extending models with new components without compromising execution time" is the headline duty. Today you build a miniature YAC.*

| # | Exercise | Time | Why it's here |
|---|---|---|---|
| 09 | `09_mpmd_components` | 90 m | Split `MPI_COMM_WORLD` by colour into "atmosphere" and "ocean" running concurrently on disjoint ranks; build an inter-communicator; exchange a field. |
| 10 | `10_conservative_remap` | 120 m | First-order **conservative** remapping between mismatched grids. Verify that global integral is preserved to machine precision. Non-conservative coupling drifts the climate — this is why couplers exist. |
| 11 | `11_sphere_interpolation` | 90 m | Source-cell search on the sphere (bucket/kd-tree), bilinear and nearest-neighbour interpolation, great-circle distance. Measure search cost vs field-application cost. |
| 12 | `12_coupling_timestep_lag` | 90 m | Atmosphere at dt=1, ocean at dt=4, ice sheet at dt=100. Accumulate/average fields across the mismatch, keep it conservative, keep it restartable. **This is the ice-sheet problem from the posting.** |

**Block C reading:** YAC's interpolation stack and the concept of a coupling configuration (what gets exchanged, on what grid, how often, with what interpolation). Start at the [YAC documentation](https://dkrz-sw.gitlab-pages.dkrz.de/yac/), then find the YAC paper in *Geoscientific Model Development* (search "YAC coupling software Earth system modelling Hanke Redler") and read the abstract and the interpolation section.

**End-of-day deliverable:** `notes/day3.md` — conservation error table (should be ~1e-16 relative), and a timeline diagram of who waits for whom in Ex. 12.

---

## Day 4 — Fri 28 Aug · Performance forensics

*Rationale: they are hiring you to make **other people's** code faster. That's a diagnostic skill, not a coding skill.*

| # | Exercise | Time | Why it's here |
|---|---|---|---|
| 13 | `13_scaling_harness` | 90 m | Reusable strong/weak scaling driver + SLURM script + plot. Compute parallel efficiency and Karp–Flatt. You'll reuse this all week. |
| 14 | `14_stencil_roofline` | 90 m | A memory-bound stencil. Measure arithmetic intensity, place it on a roofline, then optimise it and show it move. |
| 15 | `15_pmpi_profiler` | 120 m | Write a **PMPI interposition library** that counts and times every MPI call with zero source changes to the application. This is the single most impressive thing on this list. |
| 16 | `16_load_imbalance` | 90 m | A deliberately imbalanced workload. Diagnose it from timers alone, then fix it. Report the "before" and "after" imbalance ratio. |

**Block C reading:** How to read a Score-P / Vampir trace; what "late sender" and "wait at barrier" look like. Amdahl vs Gustafson framed for a climate model (why weak scaling is the metric that matters for ESMs).

**End-of-day deliverable:** `notes/day4.md` — a profile from your own PMPI tool applied to Exercise 06 or 09.

---

## Day 5 — Sat 29 Aug · Directive-based GPU + parallel I/O

*Rationale: you know CUDA; ESMs use OpenACC/OpenMP-target. Same hardware, completely different idiom. Also: a coupled run writes petabytes, so I/O is a first-class performance problem.*

> **Note:** Exercises 17–19 need a GPU. Each ships a CPU-host fallback so you can develop the code locally on macOS (`-fopenmp` host, `-acc=multicore`), but run the real thing on your cluster. If you have Leonardo/Marconi access, use it — `low_level_optimisations/leonardo_booster` has your existing setup.

| # | Exercise | Time | Why it's here |
|---|---|---|---|
| 17 | `17_openmp_target` | 90 m | Port the Day-4 stencil to `!$omp target teams distribute parallel do`. Get the **data regions** right — unnecessary host↔device transfer is the #1 ESM porting bug. |
| 18 | `18_openacc_port` | 90 m | Same kernel in OpenACC (`!$acc parallel loop`, `!$acc data`, `!$acc update`). **ICON's GPU port is OpenACC.** Compare against your CUDA instinct and against Ex. 17. |
| 19 | `19_gpu_aware_mpi_halo` | 90 m | Halo exchange straight from device pointers; overlap with interior compute; measure vs the staging-through-host version. |
| 20 | `20_parallel_io` | 90 m | Collective MPI-IO write of a decomposed field with a subarray datatype, then parallel HDF5. Compare against file-per-rank. Sweep stripe/collective-buffering. |

**Block C reading:** OpenACC vs OpenMP-target maturity for Fortran; why ESMs picked directives over CUDA (portability + scientists must be able to read the code). Unified memory vs explicit data movement.

**End-of-day deliverable:** `notes/day5.md` — kernel bandwidth CPU vs GPU, and the GPU-aware-MPI speedup number.

---

## Day 6 — Sun 30 Aug · Software engineering + capstones

*Rationale: "scientific software design" is on their wanted list. And an RSE who can't build, test, and integrate is just a programmer.*

| # | Exercise | Time | Why it's here |
|---|---|---|---|
| 21 | `21_cmake_mixed_build` | 75 m | CMake for mixed Fortran+C+MPI+OpenMP with Fortran module dependency ordering. Every ESM lives or dies by its build system. |
| 22 | `22_regression_harness` | 75 m | Bit-identical restart test + tolerance-based field comparison + a CI script. **"Testing" is in the job description.** |
| 23 | `23_add_new_component` | 150 m | **Capstone A.** Add a slow "ice sheet" component to the Ex. 09 coupled model *without slowing down the coupled timeline*. Asynchronous coupling, non-blocking exchange, sub-communicators. This is the posting's headline task, made concrete. |
| 24 | `24_mini_esm_capstone` | 150 m | **Capstone B.** Everything at once: unstructured atmosphere + ocean, partitioned, hybrid MPI+OpenMP, conservatively coupled, one offloaded kernel, scaling report, regression test. |

**Block C:** Write `notes/INTERVIEW.md`. Consolidate every number you produced this week into talking points.

---

## Exercise index

| # | Directory | Language | Core skill |
|---|---|---|---|
| 00 | `ring_allreduce.c` | C | MPI ring algorithms *(pre-existing)* |
| 01 | `01_fortran_mpi_halo_1d` | Fortran | `mpi_f08`, non-blocking halo |
| 02 | `02_fortran_omp_reduction` | Fortran | OpenMP reductions, false sharing |
| 03 | `03_hybrid_mpi_omp` | Fortran | Thread levels, comm/compute overlap |
| 04 | `04_c_fortran_interop` | Fortran + C | `iso_c_binding`, array layout |
| 05 | `05_icosahedral_grid` | Fortran | Unstructured connectivity |
| 06 | `06_unstructured_halo` | Fortran | Irregular halo, `Neighbor_alltoallv` |
| 07 | `07_partitioning` | Fortran | SFC partitioning, edge cut |
| 08 | `08_nproma_blocking` | Fortran | ICON memory layout, vectorisation |
| 09 | `09_mpmd_components` | Fortran | Communicator splitting, inter-comms |
| 10 | `10_conservative_remap` | Fortran | Conservative remapping |
| 11 | `11_sphere_interpolation` | C | Spatial search, interpolation |
| 12 | `12_coupling_timestep_lag` | Fortran | Multi-rate coupling, restart |
| 13 | `13_scaling_harness` | Fortran + Python | Scaling methodology |
| 14 | `14_stencil_roofline` | C | Roofline, memory-bound optimisation |
| 15 | `15_pmpi_profiler` | C | PMPI interposition |
| 16 | `16_load_imbalance` | Fortran | Imbalance diagnosis + repair |
| 17 | `17_openmp_target` | Fortran | OpenMP GPU offload |
| 18 | `18_openacc_port` | Fortran | OpenACC |
| 19 | `19_gpu_aware_mpi_halo` | Fortran | Device-pointer MPI |
| 20 | `20_parallel_io` | Fortran | MPI-IO, parallel HDF5 |
| 21 | `21_cmake_mixed_build` | CMake | Mixed-language builds |
| 22 | `22_regression_harness` | Fortran + Bash | Restart & tolerance testing |
| 23 | `23_add_new_component` | Fortran | **Capstone A** — async component coupling |
| 24 | `24_mini_esm_capstone` | Fortran | **Capstone B** — full mini ESM |
| **25** | **`25_icon_code_review`** | **Fortran/OpenACC** | **ICON-style code review — the interview format itself** |

---

## If you fall behind

Cut in this order — these are ranked by how likely they are to come up in a DKRZ interview:

**Never cut:** **25**, 17, 08, 09, 10, 12, 23 — 25 is the interview format, the rest are the job.
**Cut first:** 02 (you know OpenMP), 14 (you've done roofline already), 21, 11.
**Cut second:** 05 (read about the grid instead of building it), 18 *or* 17 (do one directive model, not both), 24 (23 is the better story), 01 (you know MPI halos).

A minimum viable week is **25, 17, 19, 08, 06, 09, 10, 12, 23** — about 16 hours, and it covers every theme the interview reportedly touched.

If you only have **one day**: Exercise 25 end to end, then Exercise 17, then skim Exercise 23's header. Exercise 25 alone is worth more than any three others for this specific interview.

---

## Theory questions

Coming in a separate pass once the exercises are underway. They'll cover: MPI semantics and progress, OpenMP memory model, cache/NUMA, unstructured-grid parallelism, coupling and conservation, GPU offload, Fortran-specific gotchas, and the DKRZ/ICON/YAC domain context.
