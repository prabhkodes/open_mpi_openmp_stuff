! ===========================================================================
! EXERCISE 16 — Diagnosing and fixing load imbalance
! ===========================================================================
!
! GOAL
!   A workload where every rank owns the same NUMBER of cells but not the
!   same amount of WORK. Diagnose it from timers alone, quantify how much of
!   the machine it wastes, then fix it three ways and measure each:
!     A  equal cell count          (the broken baseline)
!     B  cost-weighted partition   (static, needs a cost model)
!     C  measured-cost repartition (adaptive, uses last step's timings)
!     D  over-decomposition        (many small chunks, dynamically assigned)
!
! WHY (DKRZ)
!   This is the most common real performance bug in an Earth System Model,
!   and the reason is physical, not sloppy programming:
!
!     - radiation is only computed on the sunlit hemisphere
!     - deep convection only fires where the column is unstable
!     - sea ice thermodynamics only run where there is ice
!     - land surface does nothing over ocean
!
!   Every one of those makes cost depend on WHERE a cell is and WHAT the
!   weather is doing there — so a partition that is perfectly balanced by
!   cell count is badly imbalanced by time, and it changes as the model runs.
!
!   The killer detail: imbalance is invisible in a normal profile. The slow
!   rank looks busy and healthy; the fast ranks look like they are "spending
!   time in MPI_Allreduce", so a naive reading blames the network. You have
!   to invert it: time in a collective is time spent WAITING, and the rank
!   with the least of it is the culprit. Exercise 15 built the tool; this
!   exercise is the diagnosis.
!
! TASKS
!   TODO 1  imbalance_metrics — max/mean, wasted core-seconds, % lost
!   TODO 2  partition_weighted — split by cumulative COST, not count
!   TODO 3  measured_repartition — use last step's timings as the weights
!   TODO 4  over-decomposition with a dynamic supervisor
!   TODO 5  answer the questions at the bottom
!
! ACCEPTANCE
!   - baseline shows an imbalance ratio well above 1.0 (with the default
!     cost model, expect roughly 1.5-2x)
!   - weighted partitioning brings the ratio close to 1.0 and you can quote
!     the wall-time improvement
!   - you can state which rank was slow and how you knew, WITHOUT looking at
!     the cost model
!   - you can say when the adaptive version beats the static one, and what
!     it costs
!
! HINTS
!   - The right metric is not the ratio, it is the WASTE:
!         wasted = sum over ranks of (t_max - t_rank)
!     because that is core-seconds you paid for and did not use.
!   - Weighted partitioning is a prefix-sum problem: compute the cumulative
!     cost, then cut it into P equal-cost pieces. O(n), no iteration needed.
!   - The adaptive version has a chicken-and-egg problem: you need timings
!     to repartition, and repartitioning invalidates the timings. Repartition
!     every N steps, not every step, and say how you chose N.
!   - Over-decomposition is the robust answer when cost is unpredictable:
!     make many more chunks than ranks and hand them out on demand. You have
!     already built this pattern in open_mpi_openmp_stuff/job_scheduler --
!     reuse the idea.
! ===========================================================================

module imbalance_mod
  use mpi_f08
  implicit none
  integer, parameter :: dp = kind(1.0d0)

contains

  !> Cost of one cell, in arbitrary work units.
  !>
  !> Modelled on radiation: expensive on the sunlit hemisphere, cheap on the
  !> night side, with a sharp terminator. The important property is that cost
  !> is a function of POSITION, so a contiguous block partition gives some
  !> ranks all the expensive cells.
  pure function cell_cost(gidx, nglobal) result(w)
    integer, intent(in) :: gidx, nglobal
    real(dp) :: w, x
    x = real(gidx - 1, dp) / real(nglobal, dp)     ! position in [0,1)
    if (x < 0.5_dp) then
       w = 10.0_dp        ! day side: radiation runs
    else
       w = 1.0_dp         ! night side: it does not
    end if
    ! A little extra structure so the answer is not trivially "cut at 0.5".
    w = w + 4.0_dp * exp(-((x - 0.75_dp)/0.05_dp)**2)   ! a convecting band
  end function cell_cost

  !> Burn the given number of work units. Stands in for the physics.
  subroutine do_cell_work(units)
    real(dp), intent(in) :: units
    real(dp) :: acc
    integer  :: i, n
    n = int(units * 300.0_dp)
    acc = 0.0_dp
    do i = 1, n
       acc = acc + sqrt(real(i, dp))
    end do
    if (acc < 0.0_dp) print *, acc
  end subroutine do_cell_work

  !> Baseline: equal CELL COUNT per rank. Balanced by the wrong quantity.
  subroutine partition_by_count(nglobal, nprocs, rank, lo, hi)
    integer, intent(in)  :: nglobal, nprocs, rank
    integer, intent(out) :: lo, hi
    integer :: base, rem
    base = nglobal / nprocs
    rem  = mod(nglobal, nprocs)
    lo = rank * base + min(rank, rem) + 1
    hi = lo + base - 1
    if (rank < rem) hi = hi + 1
  end subroutine partition_by_count

  ! ------------------------------------------------------------------------
  ! TODO 2: partition by equal COST.
  !
  !   1. build the cumulative cost  C(i) = sum of cell_cost(1..i)
  !   2. target for rank r is  C_total * r / nprocs
  !   3. find the cut points by scanning (or binary searching) C
  !
  ! O(n) with a prefix sum. Note that ranks now own DIFFERENT numbers of
  ! cells -- the day-side ranks get far fewer. That is the whole point, and
  ! it is also why this breaks any code that assumed n_local was uniform.
  ! ------------------------------------------------------------------------
  subroutine partition_by_cost(nglobal, nprocs, rank, weights, lo, hi)
    integer,  intent(in)  :: nglobal, nprocs, rank
    real(dp), intent(in)  :: weights(:)
    integer,  intent(out) :: lo, hi

    ! TODO 2: replace this fallback with the real cost-based split.
    call partition_by_count(nglobal, nprocs, rank, lo, hi)
    if (size(weights) < 0) continue
  end subroutine partition_by_cost

  ! ------------------------------------------------------------------------
  ! TODO 1: turn per-rank timings into the numbers you would show a
  ! scientist.
  !
  !   t_max      slowest rank -- this IS the step time
  !   t_mean     average
  !   ratio      t_max / t_mean          (1.0 = perfect)
  !   wasted     sum over ranks of (t_max - t_rank)   [core-seconds]
  !   pct_lost   100 * (t_max - t_mean) / t_max
  !
  ! Report `wasted` in core-seconds, then convert to core-HOURS for a
  ! realistic run length. That conversion is what makes the argument land:
  ! "this costs 40000 core-hours per simulated year" is a budget line, and
  ! "the imbalance ratio is 1.6" is a statistic nobody acts on.
  ! ------------------------------------------------------------------------
  subroutine imbalance_metrics(t_local, comm, t_max, t_mean, ratio, wasted, &
                               pct_lost, slow_rank)
    real(dp),       intent(in)  :: t_local
    type(MPI_Comm), intent(in)  :: comm
    real(dp),       intent(out) :: t_max, t_mean, ratio, wasted, pct_lost
    integer,        intent(out) :: slow_rank
    integer  :: nprocs, e
    real(dp) :: t_sum

    call MPI_Comm_size(comm, nprocs, e)
    call MPI_Allreduce(t_local, t_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, e)
    call MPI_Allreduce(t_local, t_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, e)
    t_mean = t_sum / real(nprocs, dp)

    ! TODO 1: ratio, wasted, pct_lost
    ratio = 0.0_dp; wasted = 0.0_dp; pct_lost = 0.0_dp

    ! TODO 1b: find WHICH rank was slowest, with MPI_MAXLOC.
    slow_rank = -1
  end subroutine imbalance_metrics

end module imbalance_mod


program imbalance
  use imbalance_mod
  use mpi_f08
  implicit none

  type(MPI_Comm) :: comm
  integer :: rank, nprocs, ierr
  integer :: nglobal = 4000, nsteps = 10
  integer :: lo, hi, i, step
  real(dp), allocatable :: weights(:)
  real(dp) :: t0, t_local, t_bar

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)
  call read_cli()

  allocate(weights(nglobal))
  do i = 1, nglobal
     weights(i) = cell_cost(i, nglobal)
  end do

  if (rank == 0) then
     print '(a)', '=== Exercise 16: load imbalance ==='
     print '(a,i0,a,i0)', 'cells ', nglobal, '   ranks ', nprocs
     print '(a,f10.1)',   'total cost units : ', sum(weights)
     print '(a,f10.1,a,f6.1)', 'cost per cell    : min ', minval(weights), &
          '   max ', maxval(weights)
     print '(a)', '(cost model: radiation on the day side + a convecting band)'
     print '(a)', ''
  end if

  ! ---- variant A: equal cell count (the broken baseline) -----------------
  call partition_by_count(nglobal, nprocs, rank, lo, hi)
  call run_and_report('A: equal cell count   ', lo, hi)

  ! ---- variant B: equal cost ---------------------------------------------
  call partition_by_cost(nglobal, nprocs, rank, weights, lo, hi)
  call run_and_report('B: cost-weighted      ', lo, hi)

  ! ---- variant C: repartition from MEASURED time -------------------------
  call measured_repartition(lo, hi)
  call run_and_report('C: measured-cost      ', lo, hi)

  if (rank == 0) then
     print '(a)', ''
     print '(a)', '  A should be badly imbalanced; B close to 1.0.'
     print '(a)', '  If B is no better than A, partition_by_cost is still'
     print '(a)', '  the count-based fallback (TODO 2).'
     print '(a)', ''
     print '(a)', '  Then: run this under your Exercise 15 profiler and check'
     print '(a)', '  that its imbalance analysis fingers the same rank.'
  end if

  deallocate(weights)
  call MPI_Finalize(ierr)

contains

  subroutine run_and_report(label, lo_, hi_)
    character(len=*), intent(in) :: label
    integer,          intent(in) :: lo_, hi_
    real(dp) :: t_max, t_mean, ratio, wasted, pct_lost
    integer  :: slow_rank, e

    call MPI_Barrier(comm, e)
    t0 = MPI_Wtime()
    do step = 1, nsteps
       do i = lo_, hi_
          call do_cell_work(weights(i))
       end do
    end do
    t_local = MPI_Wtime() - t0

    ! Time spent in the barrier IS the wasted time -- the fast ranks sit
    ! here. This is the signal Exercise 15's profiler picks up.
    t0 = MPI_Wtime()
    call MPI_Barrier(comm, e)
    t_bar = MPI_Wtime() - t0

    call imbalance_metrics(t_local, comm, t_max, t_mean, ratio, wasted, &
                           pct_lost, slow_rank)

    if (rank == 0) then
       print '(a,a,a,f8.4,a,f8.4,a,f7.3,a)', '  ', label, &
            '  t_max ', t_max, '  t_mean ', t_mean, '  ratio ', ratio, ''
       print '(a,f9.4,a,f6.1,a,i0)', &
            '                          wasted ', wasted, &
            ' core-s   lost ', pct_lost, ' %   slowest rank ', slow_rank
    end if
  end subroutine run_and_report

  ! ------------------------------------------------------------------------
  ! TODO 3: repartition using MEASURED per-cell time rather than a model.
  !
  ! The cost model in cell_cost() is a fiction — in a real model you do not
  ! have one, because cost depends on the weather. What you DO have is last
  ! step's timings.
  !
  ! Approach:
  !   - time each cell (or each small group of cells) during a normal step
  !   - MPI_Allgatherv the per-cell times so every rank has the global cost
  !     vector
  !   - run the same prefix-sum split as TODO 2 on the MEASURED costs
  !
  ! Then answer: how often should you do this? Every step is too expensive
  ! (the Allgatherv and the data migration are not free) and never is too
  ! rarely (the day side moves). Pick an interval and justify it.
  ! ------------------------------------------------------------------------
  subroutine measured_repartition(lo_, hi_)
    integer, intent(inout) :: lo_, hi_
    ! TODO 3: measure, Allgatherv, prefix-sum split.
    ! Falls through unchanged for now, so variant C == variant B.
    continue
  end subroutine measured_repartition

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) nglobal
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) nsteps
    end if
  end subroutine read_cli

end program imbalance

! ===========================================================================
! TODO 5 — write your answers in notes/day4.md
!
! (a) For variant A, predict the imbalance ratio from the cost model BEFORE
!     running: rank 0 owns the first nglobal/P cells, all at cost 10; the
!     last rank owns cells at cost 1. Work out t_max/t_mean analytically and
!     compare with the measurement.
!
! (b) Convert the wasted core-seconds into core-hours for a realistic run:
!     1000 ranks, 365 simulated days at 100 steps/day. Would that number get
!     a scientist's attention? This conversion is the actual deliverable of
!     an imbalance investigation.
!
! (c) You are handed a model with no cost model and no source access to the
!     physics. Describe how you would establish that imbalance is the
!     problem, using ONLY MPI-level information. (You built the tool in
!     Exercise 15 -- what exactly would you look at, and what would rule
!     imbalance OUT?)
!
! (d) Cost-weighted partitioning gives ranks different cell counts. Name two
!     things elsewhere in a model that break when n_local stops being
!     uniform. (Think about anything sized at compile time, and about the
!     nproma blocking from Exercise 08.)
!
! (e) The day side MOVES. A static cost-weighted partition computed at model
!     start is wrong 12 hours later. Options: (i) repartition periodically,
!     (ii) over-decompose and assign dynamically, (iii) split radiation onto
!     its own decomposition. Give one advantage and one cost for each. ICON
!     actually does (iii) -- explain why that is attractive despite needing
!     a whole extra communication step.
!
! (f) Now the coupled version, which is Exercise 23's problem: imbalance
!     BETWEEN components, not within one. If the atmosphere finishes its
!     step 20% before the ocean, you cannot fix that by repartitioning
!     inside either component. What CAN you do? Name two approaches.
! ===========================================================================
