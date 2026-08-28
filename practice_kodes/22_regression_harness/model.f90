! ===========================================================================
! EXERCISE 22 — A regression test harness for a coupled model
! ===========================================================================
!
! GOAL
!   Build the three tests every ESM needs in CI, and make this deliberately
!   buggy toy model pass them:
!     1. RESTART      — stop at step N, restart, continue: bit-identical
!     2. REPRODUCE    — same run twice: bit-identical
!     3. INVARIANCE   — same run on a different rank count: identical to
!                       within a stated tolerance (NOT bit-identical -- and
!                       knowing why is half the exercise)
!
! WHY (DKRZ)
!   "further developing, TESTING and optimising existing coupled Earth
!   System Models" -- testing is named explicitly in the posting.
!
!   These three tests are the backbone of every climate model's CI, for a
!   reason that is specific to the domain: you cannot validate a climate
!   model against a known right answer, because there isn't one. What you
!   CAN do is pin down that the model has not changed unexpectedly. So the
!   entire testing strategy rests on reproducibility, and any change that
!   breaks bit-identity has to be justified.
!
!   This model contains a REAL restart bug of the kind that ships in
!   production code: a piece of state that is not written to the restart
!   file. Your job is to build the harness that catches it, then fix it.
!   That order matters — write the failing test first.
!
! THE BUG
!   Not going to tell you which variable. That is the exercise. The restart
!   test will fail; the diagnostic you build in TODO 2 should point at it.
!
! TASKS
!   TODO 1  write_restart / read_restart — find and fix the missing state
!   TODO 2  a diagnostic that localises WHICH variable diverged
!   TODO 3  make the rank-invariance test pass (or explain why it cannot)
!   TODO 4  run_tests.sh — wire it all into a CI script
!   TODO 5  answer the questions at the bottom
!
! ACCEPTANCE
!   - ./run_tests.sh reports 3/3 passing
!   - before your fix it reports the restart test FAILING (verify this
!     first -- a test that never fails is not a test)
!   - you can state which variable was missing and how your diagnostic
!     located it
!   - you can explain why rank-invariance is a tolerance test and not a
!     bit-identity test
!
! HINTS
!   - Unformatted stream I/O, and write the members explicitly. A formatted
!     write rounds, and then bit-identity is impossible by construction.
!   - For TODO 2: dump every state variable to a labelled text file at the
!     end of the run and diff those, rather than comparing one final number.
!     "The answer differs" is not actionable; "prev_tendency differs at step
!     501" is.
!   - The rank-invariance failure is not a bug. Global sums are computed by
!     MPI_Allreduce, whose combination order depends on rank count, and
!     floating-point addition is not associative. You met this in Exercise
!     02 -- the fix is the same, and it costs performance.
! ===========================================================================

program model
  use mpi_f08
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: NCELL = 2000

  type(MPI_Comm) :: comm
  integer :: rank, nprocs, ierr
  integer :: nsteps = 1000, restart_at = -1, i0
  character(len=32) :: mode = 'run'
  character(len=64) :: outfile = 'state.txt'

  ! ---- model state --------------------------------------------------------
  ! Everything here that evolves must be in the restart file. One of these
  ! is not. That is the bug.
  real(dp), allocatable :: u(:)          ! prognostic field
  real(dp) :: prev_tendency = 0.0_dp     ! last step's tendency (Adams-Bashforth)
  real(dp) :: accumulated_flux = 0.0_dp  ! running diagnostic
  integer  :: step_count = 0

  integer :: lo, hi, nloc

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)
  call read_cli()

  ! Decompose
  nloc = NCELL / nprocs
  lo = rank * nloc + 1
  if (rank == nprocs - 1) nloc = NCELL - lo + 1
  hi = lo + nloc - 1
  allocate(u(nloc))

  call init_state()
  i0 = 1

  if (trim(mode) == 'restart') then
     ! Run to restart_at, checkpoint, wipe, reload, continue.
     call step_range(1, restart_at)
     call write_restart('restart.bin')
     call init_state()                  ! wipe, as a fresh process would
     call read_restart('restart.bin')
     i0 = restart_at + 1
  end if

  call step_range(i0, nsteps)
  call dump_state(outfile)

  if (rank == 0) then
     print '(a,a,a,i0,a,i0)', 'mode=', trim(mode), '  steps=', nsteps, &
          '  ranks=', nprocs
     print '(a,a)', 'wrote ', trim(outfile)
  end if

  deallocate(u)
  call MPI_Finalize(ierr)

contains

  subroutine init_state()
    integer :: i
    do i = 1, nloc
       u(i) = sin(real(lo + i - 1, dp) * 0.01_dp)
    end do
    prev_tendency    = 0.0_dp
    accumulated_flux = 0.0_dp
    step_count       = 0
  end subroutine init_state

  !> One timestep. Deliberately uses a two-level scheme so that
  !> prev_tendency is genuine state and not merely a scratch variable.
  subroutine step_range(ifrom, ito)
    integer, intent(in) :: ifrom, ito
    integer  :: s, i, e
    real(dp) :: local_sum, global_sum, tendency

    do s = ifrom, ito
       local_sum = 0.0_dp
       do i = 1, nloc
          local_sum = local_sum + u(i)
       end do
       ! NOTE for TODO 3: the combination order inside this Allreduce depends
       ! on the rank count, and floating-point addition is not associative.
       call MPI_Allreduce(local_sum, global_sum, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, comm, e)

       tendency = 0.001_dp * global_sum / real(NCELL, dp)

       ! Adams-Bashforth style: needs LAST step's tendency as well as this
       ! one. That makes prev_tendency part of the model state.
       do i = 1, nloc
          u(i) = u(i) + 1.5_dp * tendency - 0.5_dp * prev_tendency
       end do

       prev_tendency    = tendency
       accumulated_flux = accumulated_flux + tendency
       step_count       = step_count + 1
    end do
  end subroutine step_range

  ! ------------------------------------------------------------------------
  ! TODO 1: the restart file.
  !
  ! Write EVERY evolving variable. Read them back in EXACTLY the same order.
  ! One is currently missing from the write; find it by making the restart
  ! test fail and then using your TODO 2 diagnostic.
  !
  ! Note the field u is DISTRIBUTED: each rank writes its own slice. For
  ! this exercise one file per rank is fine (name it with the rank), but say
  ! in your notes why a real model does not do that, and what it does
  ! instead. You built the machinery for the alternative in Exercise 20.
  ! ------------------------------------------------------------------------
  subroutine write_restart(fname)
    character(len=*), intent(in) :: fname
    character(len=80) :: rf
    integer :: unit
    write(rf, '(a,a,i4.4)') trim(fname), '.', rank
    open(newunit=unit, file=rf, form='unformatted', access='stream', &
         status='replace', action='write')
    write(unit) nloc
    write(unit) u
    write(unit) accumulated_flux
    write(unit) step_count
    ! TODO 1: something evolving is missing here. Which?
    close(unit)
  end subroutine write_restart

  subroutine read_restart(fname)
    character(len=*), intent(in) :: fname
    character(len=80) :: rf
    integer :: unit, n
    write(rf, '(a,a,i4.4)') trim(fname), '.', rank
    open(newunit=unit, file=rf, form='unformatted', access='stream', &
         status='old', action='read')
    read(unit) n
    if (n /= nloc) then
       print '(a,i0,a,i0)', 'restart size mismatch: file ', n, ' expected ', nloc
       call MPI_Abort(comm, 1, ierr)
    end if
    read(unit) u
    read(unit) accumulated_flux
    read(unit) step_count
    ! TODO 1: read back whatever you added above, in the same order.
    close(unit)
  end subroutine read_restart

  ! ------------------------------------------------------------------------
  ! TODO 2: the diagnostic dump.
  !
  ! Write every state variable with a LABEL, so a diff points at a name
  ! rather than at a byte offset. Use a format with full precision
  ! (es24.16) -- anything less and two genuinely different states can look
  ! identical, which is the worst possible failure mode for a test.
  !
  ! Rank 0 writes the global reductions; that keeps the file independent of
  ! rank count, which is what makes the invariance test in TODO 3 possible
  ! at all.
  ! ------------------------------------------------------------------------
  subroutine dump_state(fname)
    character(len=*), intent(in) :: fname
    integer  :: unit, e
    real(dp) :: gmin, gmax, gsum, lsum

    lsum = sum(u)
    call MPI_Reduce(lsum,      gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, e)
    call MPI_Reduce(minval(u), gmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, e)
    call MPI_Reduce(maxval(u), gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, e)

    if (rank /= 0) return
    open(newunit=unit, file=fname, status='replace', action='write')
    write(unit,'(a,es24.16)') 'u_sum            ', gsum
    write(unit,'(a,es24.16)') 'u_min            ', gmin
    write(unit,'(a,es24.16)') 'u_max            ', gmax
    write(unit,'(a,es24.16)') 'accumulated_flux ', accumulated_flux
    write(unit,'(a,es24.16)') 'prev_tendency    ', prev_tendency
    write(unit,'(a,i0)')      'step_count       ', step_count
    ! TODO 2: add anything else that would help localise a divergence.
    close(unit)
  end subroutine dump_state

  subroutine read_cli()
    character(len=64) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); mode = trim(arg)
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) nsteps
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg); read(arg,*) restart_at
    end if
    if (command_argument_count() >= 4) then
       call get_command_argument(4, arg); outfile = trim(arg)
    end if
    if (trim(mode) == 'restart' .and. restart_at < 1) restart_at = nsteps / 2
  end subroutine read_cli

end program model

! ===========================================================================
! TODO 5 — write your answers in notes/day6.md
!
! (a) Which variable was missing from the restart file? How did your
!     diagnostic localise it, and how long would it have taken you without
!     one? (Be honest -- this is the argument for building the diagnostic.)
!
! (b) The rank-invariance test cannot pass bit-identically. Explain why in
!     terms of MPI_Allreduce and floating-point associativity. Then: what
!     tolerance IS appropriate, and how would you choose it defensibly
!     rather than by picking a number that makes the test pass?
!
! (c) Implement a reproducible global sum (Exercise 02, TODO 4) and re-run
!     the invariance test. Does it pass bit-identically now? What did it
!     cost in runtime? Would you turn it on in production, or only in CI?
!
! (d) Restart at a step that is NOT a multiple of any internal period, and
!     at one that is. Does the bug appear in both cases? Why are bugs that
!     appear only at some restart points especially dangerous in a model
!     that restarts every 30 simulated days?
!
! (e) Design the CI suite you would want for a coupled ESM. List the tests,
!     what each catches, and roughly what each costs in core-hours. Which
!     run on every commit, which nightly, which only before a release?
!
! (f) DKRZ context: a scientist changes a physics parameterisation, and the
!     restart test now fails bit-identity. That is EXPECTED -- the physics
!     changed. So how do you distinguish "expected change" from "new bug"?
!     Describe the workflow, including what gets reviewed and by whom.
! ===========================================================================
