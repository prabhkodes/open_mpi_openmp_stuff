# Exercise 24 — CAPSTONE B: a miniature coupled Earth System Model

**This is a specification, not a skeleton.** By Day 6 you have built every
piece; this is where you assemble them and produce the artefact you talk
about in the interview.

Budget: ~2.5 hours. If you are short on time, **Exercise 23 is the better
interview story** — do that one properly and treat this as a stretch goal.

---

## What you are building

A coupled model with two components on an unstructured mesh, conservatively
coupled, hybrid-parallel, with one offloaded kernel, a scaling study, and a
regression test. Small enough to finish, complete enough to be real.

```
                  MPI_COMM_WORLD
        ┌──────────────────┴──────────────────┐
   ATMOSPHERE  (ranks 0..na-1)          OCEAN  (ranks na..)
   ─────────────────────────           ──────────────────────
   unstructured mesh, N_a cells        unstructured mesh, N_o cells
   partitioned by SFC        (Ex 07)   partitioned by SFC        (Ex 07)
   nproma-blocked fields     (Ex 08)   nproma-blocked fields     (Ex 08)
   halo exchange             (Ex 06)   halo exchange             (Ex 06)
   MPI + OpenMP hybrid       (Ex 03)   MPI + OpenMP hybrid       (Ex 03)
   one GPU-offloaded kernel  (Ex 17)
        │                                       │
        └──────── conservative coupling ────────┘
                  every ocean timestep
                  remap weights            (Ex 10)
                  time-accumulated flux    (Ex 12)
                  component split          (Ex 09)
```

## Required components

### 1. Mesh and decomposition
- Two unstructured meshes with **different cell counts** (e.g. 40000 and 25000).
  They must not be aligned — that is what makes the coupling real.
- Partition each across its component's ranks with a space-filling curve.
- Build the halo schedule once at setup, replay it every step.
- **Report:** edge cut and load imbalance for each component.

### 2. Fields and kernels
- Store prognostic fields as `f(nproma, nlev, nblks)`.
- One horizontal kernel using indirect neighbour addressing.
- One vertical/column kernel with a level dependency.
- Thread over blocks with OpenMP.
- **Report:** the `nproma` you chose and why.

### 3. Coupling
- Build the conservative remap weights **once** at setup.
- Atmosphere → ocean: a surface flux, time-accumulated over the ocean's
  timestep and handed over as a window average.
- Ocean → atmosphere: a surface temperature, remapped conservatively.
- **Report:** the conservation error. It must be ~1e-15 relative. If it is
  not, nothing else in this exercise counts.

### 4. GPU offload
- Port the column kernel to OpenMP target or OpenACC.
- Data-resident across the whole timeloop — no per-step mapping.
- CPU fallback so it still runs on macOS.
- **Report:** kernel time CPU vs GPU, and bytes crossing the bus per step.

### 5. Scaling study
- Strong and weak scaling with the Exercise 13 harness.
- **Report:** parallel efficiency, and the rank count where it breaks down.

### 6. Tests
- Restart bit-identity (Exercise 22).
- Conservation across the coupling interface, checked every step.
- **Report:** `./run_tests.sh` output showing everything green.

---

## Suggested file layout

```
24_mini_esm_capstone/
├── SPEC.md                  this file
├── Makefile
├── src/
│   ├── kinds.f90            precision
│   ├── mesh_mod.f90         mesh + SFC partition + halo schedule   (Ex 05-07)
│   ├── field_mod.f90        nproma-blocked field type              (Ex 08)
│   ├── kernels_mod.f90      horizontal + column kernels, offloaded (Ex 08,17)
│   ├── remap_mod.f90        conservative weights + apply           (Ex 10)
│   ├── couple_mod.f90       component split, accumulate, exchange  (Ex 09,12)
│   └── mini_esm.f90         the driver
├── run_tests.sh                                                    (Ex 22)
└── notes.md                 your report
```

---

## Order of work

Build it in this order. Each step is independently testable, and testing
each before moving on is the difference between finishing and debugging a
tangle at hour four.

1. **Mesh + partition** (30 min) — verify: edge cut is sane, every cell owned once.
2. **Fields + kernels, serial** (20 min) — verify: against a reference implementation.
3. **Halo exchange** (30 min) — verify: global-index poisoning check from Ex 06.
4. **Component split** (15 min) — verify: both components report the right rank sets.
5. **Remap weights** (30 min) — verify: constant field survives; global integral preserved.
6. **Coupled timeloop** (20 min) — verify: conservation holds every step.
7. **GPU offload** (20 min) — verify: checksum matches CPU.
8. **Scaling + tests** (25 min) — produce the numbers.

**Do not proceed past a step whose verification fails.** A conservation bug
introduced at step 5 and discovered at step 8 will cost you the rest of the
session.

---

## What "done" looks like

A `notes.md` containing:

| Metric | Your number |
|---|---|
| Cells: atmosphere / ocean | |
| Edge cut, SFC vs linear | |
| Load imbalance | |
| `nproma` chosen, and why | |
| Halo exchange, µs per step | |
| Conservation error (relative) | |
| Coupling overhead, % of step | |
| Column kernel, CPU vs GPU | |
| Strong scaling efficiency at max ranks | |
| Weak scaling efficiency at max ranks | |
| Restart bit-identical? | |

Plus one paragraph: **what you would optimise next, and why.** That question
is the one an interviewer actually cares about, because it shows whether you
can prioritise rather than just measure.

---

## The interview version

Be able to give this in 90 seconds:

> "I built a miniature coupled ESM — two unstructured meshes, different
> resolutions, partitioned with a space-filling curve, hybrid MPI+OpenMP,
> conservatively coupled with remap weights built at setup. Conservation
> holds to 1e-16 per step. The column kernel is offloaded with OpenMP target
> and data-resident across the timeloop. It weak-scales at X% to Y ranks,
> and the thing that limits it is Z."

Then have the numbers ready when they ask.

---

## Stretch goals

Only after everything above is green:

- Add the ice sheet from Exercise 23 as a third, lagged component.
- Replace the hand-rolled halo with `MPI_Neighbor_alltoallv` and compare.
- Add an I/O server: dedicate 2 ranks to writing, so compute never blocks (Ex 20).
- Run it on a real cluster and redo the scaling study at 100+ ranks.
