! ===========================================================================
! EXERCISE 13 — A reusable scaling harness
! ===========================================================================
!
! GOAL
!   Build the measurement infrastructure you will use for the rest of the
!   week: a timer with proper statistics, a strong/weak scaling driver, and
!   the derived metrics that turn raw times into an argument — speedup,
!   parallel efficiency, and the Karp-Flatt serial fraction.
!
! WHY (DKRZ)
!   The posting's job title is "HPC Code Performance Optimisation". Almost
!   everything in that role starts with someone saying "the model is slow"
!   and you having to turn that into a defensible number.
!
!   Two things separate a useful measurement from a useless one:
!
!     1. STATISTICS. One run is not a measurement. Report min/median/max
!        across repetitions, and report the MAX across ranks, not the mean —
!        a parallel program runs at the speed of its slowest rank, so the
!        mean systematically flatters you.
!
!     2. THE RIGHT METRIC. For a climate model, WEAK scaling is what
!        matters: nobody runs the same grid on more nodes, they run a finer
!        grid. Quoting strong-scaling efficiency for an ESM is a category
!        error, and an interviewer will notice.
!
!   Karp-Flatt is the one people do not know. It backs the serial fraction
!   out of measured speedup, and unlike Amdahl it exposes whether your
!   scaling loss is a fixed serial section (f stays constant) or growing
!   overhead (f rises with P). That distinction tells you where to look.
!
! TASKS
!   TODO 1  timer_t          — accumulate, min/median/max, MPI reduction
!   TODO 2  parallel_metrics — speedup, efficiency, Karp-Flatt
!   TODO 3  strong_scaling   — fixed global problem
!   TODO 4  weak_scaling     — fixed problem PER RANK
!   TODO 5  emit CSV, and plot it with plot_scaling.py
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - `make strong` and `make weak` produce clean CSV
!   - `make plot` produces scaling.png with both curves and the ideal lines
!   - the Karp-Flatt fraction is computed and you can interpret its trend
!   - you reuse this harness in at least two other exercises this week
!
! HINTS
!   - Always MPI_Barrier before starting a timer, or you measure the
!     previous phase's load imbalance instead of this one's.
!   - Use MPI_Wtime, not cpu_time or system_clock. cpu_time sums over
!     threads; system_clock can have terrible resolution.
!   - Discard the first repetition. First touch, page faults, and cache
!     warm-up are real and they are not what you are trying to measure.
!   - Karp-Flatt:  e = (1/S - 1/P) / (1 - 1/P)   for measured speedup S on
!     P processors. If e is flat, you have a fixed serial section. If e
!     grows with P, you have overhead that scales with P — communication.
! ===========================================================================

module scaling_mod
  use mpi_f08
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: MAX_REPS = 64

  !> A timer that keeps every sample, so you can report a distribution
  !> instead of a single number.
  type :: timer_t
     character(len=32) :: name = ''
     integer  :: n = 0
     real(dp) :: t0 = 0.0_dp
     real(dp) :: samples(MAX_REPS) = 0.0_dp
  end type timer_t

contains

  subroutine timer_start(t)
    type(timer_t), intent(inout) :: t
    t%t0 = MPI_Wtime()
  end subroutine timer_start

  ! ------------------------------------------------------------------------
  ! TODO 1: record one sample.
  !   t%n = t%n + 1;  t%samples(t%n) = MPI_Wtime() - t%t0
  ! Guard against overrunning MAX_REPS.
  ! ------------------------------------------------------------------------
  subroutine timer_stop(t)
    type(timer_t), intent(inout) :: t
    ! TODO 1
    if (t%n < 0) continue
  end subroutine timer_stop

  !> Median of the recorded samples. Reported instead of the mean because a
  !> single OS hiccup skews a mean badly and a median not at all.
  real(dp) function timer_median(t) result(m)
    type(timer_t), intent(in) :: t
    real(dp) :: s(MAX_REPS)
    integer  :: i, j, n
    real(dp) :: tmp
    n = t%n
    if (n == 0) then
       m = 0.0_dp; return
    end if
    s(1:n) = t%samples(1:n)
    do i = 1, n-1                       ! insertion sort; n is tiny
       do j = i+1, n
          if (s(j) < s(i)) then
             tmp = s(i); s(i) = s(j); s(j) = tmp
          end if
       end do
    end do
    if (mod(n,2) == 1) then
       m = s((n+1)/2)
    else
       m = 0.5_dp * (s(n/2) + s(n/2+1))
    end if
  end function timer_median

  !> The number that matters in parallel: the SLOWEST rank's time. A
  !> collective at the end of the timestep makes everyone wait for it, so
  !> the mean across ranks is not the cost you pay.
  real(dp) function timer_max_over_ranks(t, comm) result(tmax)
    type(timer_t),  intent(in) :: t
    type(MPI_Comm), intent(in) :: comm
    real(dp) :: local
    integer  :: e
    local = timer_median(t)
    call MPI_Allreduce(local, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, e)
  end function timer_max_over_ranks

  ! ------------------------------------------------------------------------
  ! TODO 2: the derived metrics.
  !
  !   speedup     S = t_ref / t_p
  !   efficiency  E = S / (P / P_ref)
  !   Karp-Flatt  e = (1/S - 1/P) / (1 - 1/P)      (P > 1)
  !
  ! For weak scaling the definitions change: "efficiency" is t_ref / t_p
  ! directly, because the work per rank is constant and ideal is a FLAT
  ! line, not a rising one. Handle both and label them clearly — mixing them
  ! up is the most common way to produce a misleading scaling plot.
  ! ------------------------------------------------------------------------
  subroutine parallel_metrics(t_ref, t_p, p_ref, p, weak, speedup, eff, karp)
    real(dp), intent(in)  :: t_ref, t_p
    integer,  intent(in)  :: p_ref, p
    logical,  intent(in)  :: weak
    real(dp), intent(out) :: speedup, eff, karp
    speedup = 0.0_dp; eff = 0.0_dp; karp = 0.0_dp
    ! TODO 2
    if (t_ref < 0.0_dp .or. t_p < 0.0_dp .or. p_ref < 0 .or. p < 0) continue
    if (weak) continue
  end subroutine parallel_metrics

  !> The workload being scaled: a distributed 5-point stencil with a halo
  !> exchange. Deliberately simple -- the point of this exercise is the
  !> measurement, not the kernel.
  subroutine workload(n_local, n_sweeps, comm)
    integer,        intent(in) :: n_local, n_sweeps
    type(MPI_Comm), intent(in) :: comm
    real(dp), allocatable :: u(:), unew(:)
    integer :: rank, nprocs, left, right, s, i, e
    type(MPI_Status) :: st

    call MPI_Comm_rank(comm, rank, e)
    call MPI_Comm_size(comm, nprocs, e)
    left  = modulo(rank-1, nprocs)
    right = modulo(rank+1, nprocs)

    allocate(u(0:n_local+1), unew(0:n_local+1))
    do i = 0, n_local+1
       u(i) = sin(real(i, dp) * 1.0e-4_dp)
    end do
    unew = u

    do s = 1, n_sweeps
       call MPI_Sendrecv(u(n_local), 1, MPI_DOUBLE_PRECISION, right, 1, &
                         u(0),       1, MPI_DOUBLE_PRECISION, left,  1, &
                         comm, st, e)
       call MPI_Sendrecv(u(1),         1, MPI_DOUBLE_PRECISION, left,  2, &
                         u(n_local+1), 1, MPI_DOUBLE_PRECISION, right, 2, &
                         comm, st, e)
       !$omp parallel do
       do i = 1, n_local
          unew(i) = 0.25_dp*u(i-1) + 0.5_dp*u(i) + 0.25_dp*u(i+1)
       end do
       !$omp end parallel do
       u = unew
    end do

    deallocate(u, unew)
  end subroutine workload

end module scaling_mod


program scaling
  use scaling_mod
  use mpi_f08
  implicit none

  type(MPI_Comm) :: comm
  type(timer_t)  :: tm
  integer :: rank, nprocs, ierr
  integer :: n_global = 4000000     ! strong scaling: fixed globally
  integer :: n_per_rank = 500000    ! weak scaling:   fixed per rank
  integer :: n_sweeps = 200, nreps = 7
  character(len=8) :: mode = 'strong'
  integer :: n_local, r
  real(dp) :: t_max

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)
  call read_cli()

  ! ---- TODO 3 / TODO 4: choose the problem size for the mode -------------
  if (trim(mode) == 'weak') then
     n_local = n_per_rank                      ! work per rank is CONSTANT
  else
     n_local = n_global / nprocs               ! work per rank SHRINKS
     if (rank < mod(n_global, nprocs)) n_local = n_local + 1
  end if

  tm%name = 'workload'

  ! Warm-up, discarded: first touch and page faults are not the measurement.
  call workload(n_local, 10, comm)

  do r = 1, nreps
     call MPI_Barrier(comm, ierr)         ! or you time the PREVIOUS imbalance
     call timer_start(tm)
     call workload(n_local, n_sweeps, comm)
     call timer_stop(tm)
  end do

  t_max = timer_max_over_ranks(tm, comm)

  ! ---- TODO 5: emit one CSV row per run ----------------------------------
  ! The Makefile concatenates these into scaling_strong.csv / scaling_weak.csv
  ! and plot_scaling.py turns them into a figure.
  if (rank == 0) then
     if (t_max <= 0.0_dp) then
        print '(a)', '# timer_stop is still a stub (TODO 1) -- times are zero'
     end if
     print '(a,a,a,i0,a,i0,a,es14.6,a,i0)', &
          trim(mode), ',', 'ranks,', nprocs, ',nlocal,', n_local, &
          ',time,', t_max, ',sweeps,', n_sweeps
  end if

  call MPI_Finalize(ierr)

contains

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); mode = trim(arg)
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) n_sweeps
    end if
  end subroutine read_cli

end program scaling

! ===========================================================================
! TODO 6 — write your answers in notes/day4.md
!
! (a) Plot both curves. At what rank count does strong scaling break down?
!     Compute n_local at that point and compare it to your L2 cache size and
!     to the halo:interior ratio. Which one explains the breakdown?
!
! (b) Karp-Flatt: is the serial fraction constant or rising with P? Say what
!     each case would mean, and which one you measured. If it is rising,
!     what specifically is growing?
!
! (c) Weak scaling should be a flat line and is not. Where does the loss
!     come from, given that each rank does identical work? (Two candidates:
!     the collective inside the halo exchange, and the fact that MPI_Sendrecv
!     latency does not shrink. Design a measurement that separates them.)
!
! (d) Amdahl vs Gustafson for a climate model. A scientist says "we get 60%
!     parallel efficiency at 1000 nodes, that's bad." What do you need to
!     know before agreeing? Frame your answer around what they actually do
!     with more nodes.
!
! (e) You reported the MAX over ranks. Re-plot using the MEAN and put both
!     on the same axes. How different are they, and what does the gap
!     measure? (That gap is the subject of Exercise 16.)
!
! (f) DKRZ context: Levante has 128 cores per node. Design the scaling study
!     you would run for a new ESM component: which node counts, which
!     ranks-per-node splits, how many repetitions, and what you would hold
!     fixed. Justify each choice in one line — this is a real deliverable
!     shape for the job.
! ===========================================================================
