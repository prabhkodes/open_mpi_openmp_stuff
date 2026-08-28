# Exercise 25 — ICON-style code review: "what would you optimise here?"

**Do not read [ANSWERS.md](ANSWERS.md) until you have written your own answer for a snippet.**

---

## What this is

Intel from someone who sat a similar DKRZ interview: they were shown a block of
ICON code and asked what needed optimising. The answers that landed were about
**memory coalescing, exploiting data parallelism, race conditions, placement of
data regions and host↔device transfers, and overlapping asynchronous
operations**.

That is a *code reading* interview, not a trivia interview. These twelve
snippets are written in ICON's idiom with defects planted in them. Most have
more than one — as real code does.

> **On provenance:** these are ICON-*flavoured*, not copies of ICON source. The
> naming, data layout and directive style follow ICON's conventions so the
> reading transfers. Verify anything you plan to assert as fact about ICON
> against the actual repository at `gitlab.dkrz.de`.

## The ICON idiom you need to read fluently


```fortran
! Fields are ALWAYS blocked:   f(nproma, nlev, nblks)
!   jc  index within a block   (1 .. nproma)   <- fastest-varying
!   jk  vertical level         (1 .. nlev)
!   jb  block                  (1 .. nblks)
!
! On CPU: nproma ~ a few hundred, nblks large, OpenMP threads over jb.
! On GPU: nproma = all cells in the patch, so nblks is 1 or 2.
!         ^^^ this single fact is behind half the defects below.

DO jb = i_startblk, i_endblk
   CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                      i_startidx, i_endidx, rl_start, rl_end)
   ...
END DO

! Neighbour access is INDIRECT:
!   iidx => p_patch%cells%neighbor_idx    ! (nproma, nblks, 3)
!   iblk => p_patch%cells%neighbor_blk
!   f(iidx(jc,jb,1), jk, iblk(jc,jb,1))

! Directives:
!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$ACC END PARALLEL
!$ACC WAIT(1)
```

## How to answer out loud

The content matters, but so does the shape. Work through it in this order —
it signals that you optimise by measurement rather than by pattern-matching:

1. **What does this code do?** One sentence. Establish you read it.
2. **Where would the time go?** Bandwidth-bound or latency-bound? How much data
   moves? Say this *before* naming a fix.
3. **Correctness first.** A race or a stale-data bug outranks any speedup. If
   you see one, lead with it.
4. **Name the mechanism, not the label.** Not "it's not coalesced" but
   "consecutive threads access addresses `nproma*8` bytes apart, so each warp
   touches 32 cache lines instead of 4."
5. **Quantify.** Even roughly. "This is a 8–32× traffic amplification."
6. **Propose the fix, and state its cost.** Every fix has one.
7. **Say how you would verify.** `-Minfo=accel`, `nsys`, a checksum.

If you do not know something, say what you would measure to find out. That is a
much better answer than a confident guess, and interviewers can tell the
difference.

---

# Part A — memory coalescing and data layout

## Q1. The loop nest

```fortran
!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
DO jc = i_startidx, i_endidx
  DO jk = 1, nlev
    z_theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) &
      &                 * ( 1._wp + vwp1*p_prog%tracer(jc,jk,jb,iqv) &
      &                           - p_prog%tracer(jc,jk,jb,iqc) )
  END DO
END DO
!$ACC END PARALLEL
```

1. This is correct code. What is wrong with it on a GPU?
2. Two consecutive threads in a warp — what addresses do they touch, and how
   far apart? Give the answer in bytes for `nproma = 20000`, `wp = real64`.
3. How many cache lines does one 32-thread warp touch here, versus the ideal?
4. Write the fixed version.
5. On a **CPU** with `nproma = 256`, is the original loop order wrong too?
   Explain the difference — and what it tells you about writing one source for
   both machines.

---

## Q2. The neighbour gather

```fortran
iidx => p_patch%cells%neighbor_idx
iblk => p_patch%cells%neighbor_blk

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
DO jk = 1, nlev
  DO jc = i_startidx, i_endidx
    p_diag%div(jc,jk,jb) =                                            &
      &   p_diag%flux(iidx(jc,jb,1), jk, iblk(jc,jb,1)) * div_coeff(jc,jb,1) &
      & + p_diag%flux(iidx(jc,jb,2), jk, iblk(jc,jb,2)) * div_coeff(jc,jb,2) &
      & + p_diag%flux(iidx(jc,jb,3), jk, iblk(jc,jb,3)) * div_coeff(jc,jb,3)
  END DO
END DO
!$ACC END PARALLEL
```

1. The loop order here is right. But the access pattern still may not coalesce.
   Why not, and what does it depend on?
2. `iidx(jc,jb,1)` is itself an array read. Count *all* the loads needed to
   produce one output element. Which of them coalesce and which do not?
3. This is the fundamental tax of an unstructured grid. Name the setup-time
   decision that determines how bad it is here, and say which exercise in this
   repo you would point at to fix it.
4. A colleague proposes caching the three neighbour values in a local array to
   "avoid re-reading". Is that a good idea on a GPU? What about on a CPU?
5. How would you *measure* whether this kernel is coalescing-limited rather
   than just bandwidth-limited?

---

## Q3. The derived type

```fortran
TYPE t_cell_state
  REAL(wp) :: rho, theta_v, exner, w
END TYPE t_cell_state

TYPE(t_cell_state), ALLOCATABLE :: cells(:,:,:)   ! (nproma, nlev, nblks)

!$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1)
DO jk = 1, nlev
  DO jc = i_startidx, i_endidx
    cells(jc,jk,jb)%rho = cells(jc,jk,jb)%rho * (1._wp - dtime*divergence(jc,jk,jb))
  END DO
END DO
```

1. Loop order is right. Coalescing is still broken. Why?
2. What fraction of every fetched cache line does this kernel actually use?
3. Name the two layouts and say which one ICON uses and why.
4. Give one case where the layout in this snippet is the *better* choice.

---

# Part B — exploiting data parallelism

## Q4. The block loop

```fortran
!$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) ASYNC(1)
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)
  DO jk = 1, nlev
    DO jc = i_startidx, i_endidx
      p_diag%temp(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) * p_diag%exner(jc,jk,jb)
    END DO
  END DO
END DO
!$ACC END PARALLEL
```

1. On the GPU build `nproma` is set to the whole patch, so `nblks = 1`. What
   does this kernel actually launch? Quantify the utilisation on a device with
   108 SMs.
2. There is a second, independent problem with putting `!$ACC PARALLEL LOOP` on
   this particular loop. What is it? (Look at what is inside the loop body.)
3. Rewrite it correctly.
4. The CPU build uses `!$OMP PARALLEL DO` over exactly this `jb` loop, and that
   is right for CPU. Explain why the same loop is the correct parallel
   dimension on one machine and the wrong one on the other.

---

## Q5. The vertical solve

```fortran
!$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1)
DO jk = 2, nlev
  DO jc = i_startidx, i_endidx
    w(jc,jk,jb) = w(jc,jk-1,jb) + dz(jc,jk,jb) * src(jc,jk,jb)
  END DO
END DO
!$ACC END PARALLEL
```

1. What is wrong with this, and what will it produce at runtime?
2. Is the bug deterministic? What would you expect to see if you ran it twice?
3. Write the correct directive. Which loop carries the parallelism?
4. After your fix, how much parallelism is available? Is it enough to fill a
   GPU? What does that depend on?
5. `nlev` is 90 and the dependency is serial. Name one algorithmic change that
   would expose more parallelism, and say what it costs.
6. Would `-Minfo=accel` have caught this? Explain what the compiler can and
   cannot prove here, and why `COLLAPSE` is dangerous in a way that a bare
   `LOOP` is not.

---

# Part C — race conditions

## Q6. Edge to cell

```fortran
ieidx => p_patch%edges%cell_idx
ieblk => p_patch%edges%cell_blk

!$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1)
DO jk = 1, nlev
  DO je = i_startidx, i_endidx
    p_diag%div(ieidx(je,jb,1), jk, ieblk(je,jb,1)) =        &
      p_diag%div(ieidx(je,jb,1), jk, ieblk(je,jb,1)) + flux(je,jk,jb)
    p_diag%div(ieidx(je,jb,2), jk, ieblk(je,jb,2)) =        &
      p_diag%div(ieidx(je,jb,2), jk, ieblk(je,jb,2)) - flux(je,jk,jb)
  END DO
END DO
!$ACC END PARALLEL
```

1. Name the defect precisely.
2. Every interior edge is shared by exactly two cells. Use that to explain why
   the race is *guaranteed*, not merely possible.
3. Give three fixes. For each, state the cost and whether the result stays
   bit-reproducible.
4. ICON computes divergence by looping over **cells** and gathering from their
   edges (as in Q2), not by looping over edges and scattering. Given that the
   edge loop does half the flux arithmetic, why is the cell loop still the
   right choice?
5. Does the same bug exist in the OpenMP CPU version of this loop? Under what
   condition would it be hidden, and why is "it works on CPU" not evidence of
   correctness?

---

## Q7. The OpenMP block loop

```fortran
!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jk = 1, nlev
    DO jc = i_startidx, i_endidx
      z_tmp(jc,jk) = p_prog%rho(jc,jk,jb) * p_diag%w(jc,jk,jb)
    END DO
  END DO

  DO jk = 2, nlev
    DO jc = i_startidx, i_endidx
      p_diag%flux(jc,jk,jb) = 0.5_wp * (z_tmp(jc,jk) + z_tmp(jc,jk-1))
    END DO
  END DO

END DO
!$OMP END PARALLEL DO
```

1. There are **two** shared-variable defects here. Find both.
2. For each, describe the symptom: crash, wrong answer, or intermittent wrong
   answer? Which is worst and why?
3. Would either defect show up reliably in a 2-thread test? A 128-thread test?
4. Write the corrected directive.
5. `ICON_OMP_DEFAULT_SCHEDULE` is a preprocessor macro that expands to a
   `SCHEDULE(...)` clause. Given what you measured in Exercise 02 and 16, what
   would you set it to for this loop, and what would change your mind?

---

## Q8. The reduction

```fortran
max_vcfl = 0._wp

!$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1)
DO jk = 1, nlev
  DO jc = i_startidx, i_endidx
    vcfl = ABS(w(jc,jk,jb)) * dtime / dz(jc,jk,jb)
    IF (vcfl > max_vcfl) max_vcfl = vcfl
  END DO
END DO
!$ACC END PARALLEL

IF (max_vcfl > 1._wp) CALL finish('velocity_tendencies', 'CFL violated')
```

1. Find all three defects. (One is a race, one is a scoping bug, one is about
   *when* the value is read.)
2. Write the corrected version.
3. This is a CFL check, so it feeds a `finish()` that aborts the run. Explain
   why a race here is more dangerous than a race in a diagnostic that only gets
   printed.
4. Across MPI ranks, this needs to become a global maximum. Where would you put
   the `MPI_Allreduce`, and what does that do to the async pipeline you built
   in Part E?
5. Reductions on GPU are combined in a non-deterministic order. For `MAX`, does
   that threaten bit-reproducibility? What if it were `SUM`? Connect this to
   Exercise 02 and Exercise 22.

---

# Part D — data regions and host↔device transfers

## Q9. The timeloop

```fortran
DO jstep = 1, nsteps

  !$ACC DATA COPYIN(p_prog%rho, p_prog%theta_v, p_metrics%ddqz_z_full) &
  !$ACC      COPYOUT(p_diag%div, p_diag%temp)

  CALL compute_divergence(p_patch, p_prog, p_diag)
  CALL compute_temperature(p_patch, p_prog, p_diag)
  CALL vertical_diffusion(p_patch, p_prog, p_diag, p_metrics)

  !$ACC END DATA

  CALL nh_solve(p_patch, p_prog, p_diag)

END DO
```

1. What is the single biggest problem here?
2. `p_prog%rho` is `(nproma, nlev, nblks)` with 5 million cells and 90 levels.
   Compute the bytes moved per timestep, and the time at 25 GB/s.
3. Rewrite the data management. Where does the data region belong, and which
   directive would you use instead of `DATA`?
4. `nh_solve` is outside the data region and runs on the host. What does that
   do to your redesign, and what would you check before assuming you can move
   it inside?
5. ICON uses `!$ACC ENTER DATA` in each module's init rather than one giant
   region around the timeloop. Give two reasons, given that the model has
   hundreds of modules.

---

## Q10. The missing clause

```fortran
SUBROUTINE compute_exner(p_patch, p_prog, p_diag, i_startblk, i_endblk)
  ...
  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
        p_diag%exner(jc,jk,jb) = (p_prog%rho(jc,jk,jb) * rd &
          &                     * p_prog%theta_v(jc,jk,jb) / p0ref)**rd_o_cpd
      END DO
    END DO
    !$ACC END PARALLEL
  END DO
END SUBROUTINE
```

The arrays *are* already on the device via an earlier `ENTER DATA`.

1. One clause is missing. Which, and what does OpenACC do without it?
2. If the data is genuinely present, does the missing clause cost anything at
   runtime? So why does ICON insist on it in every kernel?
3. Describe the failure mode this clause is designed to turn into a loud error.
   Why is that failure mode so hard to find otherwise?
4. What does `nsys` output look like for this kernel with and without the data
   present?
5. This is a general principle in performance engineering. State it in one
   sentence.

---

## Q11. The diagnostic

```fortran
!$ACC ENTER DATA COPYIN(p_prog%rho, p_prog%w)

DO jstep = 1, nsteps

  CALL dynamics_step(p_patch, p_prog, p_diag)      ! all on device, ASYNC(1)

  !$ACC UPDATE HOST(p_prog%rho, p_prog%w) ASYNC(1)

  total_mass = SUM(p_prog%rho) * cell_volume
  max_w      = MAXVAL(ABS(p_prog%w))

  IF (MOD(jstep, 100) == 0) THEN
    WRITE(message_text,'(a,e18.10)') 'total mass: ', total_mass
    CALL message('dynamics', message_text)
  END IF

END DO
```

1. There is a correctness bug and a performance bug. Find both.
2. The correctness bug produces plausible-looking numbers rather than a crash.
   Explain the mechanism, and say what you would see if you printed
   `total_mass` every step.
3. Fix both. What is the ideal frequency of the transfer, and why?
4. `SUM(p_prog%rho)` runs on the host over a 3.6 GB array. Propose a better
   design that does not transfer the field at all.
5. Now the reproducibility question: your device-side sum is combined in a
   non-deterministic order. The model must reproduce bit-identically on
   restart. Reconcile those two facts.

---

# Part E — asynchronous execution and overlap

## Q12. The halo exchange

```fortran
DO jstep = 1, nsteps

  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT)
  DO jk = 1, nlev
    DO jc = 1, nproma
      p_diag%vt(jc,jk,jb) = compute_tangential(jc,jk,jb)
    END DO
  END DO
  !$ACC WAIT

  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
  DO i = 1, n_send
    sendbuf(i) = p_diag%vt(send_idx(i), send_lev(i), send_blk(i))
  END DO
  !$ACC WAIT

  !$ACC UPDATE HOST(sendbuf)
  CALL MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
  CALL start_halo_exchange(sendbuf, recvbuf, req)
  !$ACC UPDATE DEVICE(recvbuf)

  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
  DO i = 1, n_recv
    p_diag%vt(recv_idx(i), recv_lev(i), recv_blk(i)) = recvbuf(i)
  END DO
  !$ACC WAIT

  CALL update_interior(p_patch, p_diag)

END DO
```

1. Count the synchronisation points per timestep. Which are necessary?
2. Every `!$ACC WAIT` costs a full device drain. What is the effect on kernel
   launch latency, and roughly what does that cost per timestep at 90 levels
   and 30 fields?
3. There is a **correctness** bug in the MPI sequence, independent of
   performance. Find it. (Look carefully at the order of the `Waitall` and the
   `start_halo_exchange`.)
4. The two `UPDATE` calls stage the buffers through the host. Rewrite that part
   using GPU-aware MPI. What is the exact directive, and what must be true of
   the MPI build?
5. `update_interior` does not depend on the halo. Restructure the whole loop to
   overlap it with the exchange. Write the version with async queues.
6. Draw the timeline before and after, marking device idle time.
7. What single measurement would prove your restructuring actually overlapped
   anything, rather than just moving the wait?

---

# Part F — the open question

## Q13. "What would you optimise?"

This is how the interview actually opens. There is no snippet — you will be
handed one and asked this. Prepare a **procedure**, not an answer:

1. What are the first three questions you ask *before* proposing anything?
2. What do you need to know about the machine? About the run configuration?
3. Given the choice between a 2× speedup on a kernel that is 5% of runtime and
   a 10% speedup on one that is 60%, which do you take? Show the arithmetic.
4. You profile and find the model spends 35% of its time in MPI. Name four
   distinct causes, and the measurement that distinguishes each.
5. A scientist tells you their run "got slower after we moved to the new
   machine." What is your first hour of work?
6. You find a 3× win but it changes results in the last three bits. What do you
   do? Who decides?

Write these out. Then say them out loud, timed. Aim for 60–90 seconds each.

---

## Working through this set

Suggested order, if you are short on time:

**Do first:** Q1, Q5, Q6, Q9, Q12 — one from each part, covering every theme
your friend reported.
**Then:** Q4, Q8, Q11, Q13.
**If time:** Q2, Q3, Q7, Q10.

Write your answers in `notes/day5.md` *before* opening [ANSWERS.md](ANSWERS.md).
The gap between what you wrote and what is there is your actual study list.
