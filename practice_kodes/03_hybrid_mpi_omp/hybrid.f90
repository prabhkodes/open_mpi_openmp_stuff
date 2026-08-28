! ===========================================================================
! EXERCISE 03 — Hybrid MPI+OpenMP: thread levels and communication overlap
! ===========================================================================
!
! GOAL
!   Take the halo exchange from Exercise 01 and hide it behind computation.
!   Three variants, timed against each other:
!     A  exchange, then compute everything        (no overlap — the baseline)
!     B  post Irecv/Isend, compute the INTERIOR with OpenMP, Waitall,
!        then compute the two boundary points     (overlap)
!     C  MPI_THREAD_MULTIPLE: several fields exchanged concurrently by
!        different threads
!
! WHY (DKRZ)
!   "Parallel programming with MPI/OpenMP" in the posting means hybrid, not
!   one or the other. On Levante a node has 128 cores; running 128 pure-MPI
!   ranks per node multiplies halo surface area and MPI internal buffers,
!   so ESMs run a few ranks per node with OpenMP inside. Variant B is the
!   canonical latency-hiding trick and the first thing you would reach for
!   when a scientist says "my model stops scaling past 512 nodes".
!
! TASKS
!   TODO 1  MPI_Init_thread and report what level you actually got
!   TODO 2  variant A — blocking exchange then full-domain compute
!   TODO 3  variant B — split interior/boundary compute around the Waitall
!   TODO 4  variant C — MPI_THREAD_MULTIPLE, one thread per field
!   TODO 5  compute and print the overlap efficiency
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - all three variants produce the same field checksum (printed; must match
!     to 1e-12 relative)
!   - variant B is measurably faster than A; you can quote the overlap
!     efficiency  eta = (tA - tB) / t_comm_alone
!   - you know which thread level your OpenMPI build actually provides
!
! HINTS
!   - MPI_THREAD_FUNNELED means only the thread that called MPI_Init_thread
!     may call MPI. That covers variant B: the MPI calls sit OUTSIDE the
!     parallel region, or inside a `!$omp master`.
!   - Variant B only pays off if the interior work is big enough to cover
!     the message latency. If B ~= A, increase nsweeps or n_global.
!   - MPI_THREAD_MULTIPLE is often slower per-call than FUNNELED because the
!     library must lock internally. Measure it, don't assume it.
!   - Watch out: with `!$omp parallel do` on the interior, the compiler will
!     happily also parallelise a loop containing an MPI call. Don't.
! ===========================================================================

program hybrid
  use mpi_f08
  use omp_lib
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: NFIELD = 4        ! fields exchanged in variant C

  integer :: n_global = 8000000
  integer :: nsweeps  = 50                ! stencil sweeps between exchanges
  integer :: n_iter   = 20                ! timestep count

  type(MPI_Comm) :: comm
  integer :: rank, nprocs, ierr, provided, required
  integer :: n_local, i_start, left, right
  real(dp), allocatable :: u(:), unew(:)
  real(dp), allocatable :: multi(:,:)     ! multi(1-1:n_local+1, NFIELD)
  real(dp) :: tA, tB, tC, t_comm, chkA, chkB, chkC

  ! ------------------------------------------------------------------------
  ! TODO 1: replace MPI_Init with MPI_Init_thread.
  !
  !   required = MPI_THREAD_FUNNELED
  !   call MPI_Init_thread(required, provided, ierr)
  !
  ! then check `provided >= required` and abort with a clear message if not.
  ! Later, for variant C, you need MPI_THREAD_MULTIPLE — decide whether to
  ! always request MULTIPLE or to skip variant C when it is unavailable.
  ! ------------------------------------------------------------------------
  required = MPI_THREAD_FUNNELED
  provided = -1
  call MPI_Init(ierr)                     ! TODO 1: -> MPI_Init_thread

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)
  call read_cli()

  n_local = n_global / nprocs
  if (rank < mod(n_global, nprocs)) n_local = n_local + 1
  i_start = rank * (n_global / nprocs) + min(rank, mod(n_global, nprocs)) + 1

  left  = modulo(rank - 1, nprocs)
  right = modulo(rank + 1, nprocs)

  allocate(u(0:n_local+1), unew(0:n_local+1))
  allocate(multi(0:n_local+1, NFIELD))

  if (rank == 0) then
     print '(a)', '=== Exercise 03: hybrid MPI+OpenMP overlap ==='
     print '(a,i0,a,i0,a,i0)', 'ranks=', nprocs, '  threads/rank=', &
          omp_get_max_threads(), '  n_local=', n_local
     print '(a,i0,a,i0)', 'sweeps per step=', nsweeps, '  steps=', n_iter
     print '(a,a)', 'thread level provided: ', thread_level_name(provided)
     print '(a)', ''
  end if

  ! ---- measure the bare communication cost, for the efficiency metric ----
  t_comm = time_comm_alone()

  ! ---- variant A ---------------------------------------------------------
  call init_field()
  call MPI_Barrier(comm, ierr)
  tA = MPI_Wtime()
  call variant_a()
  tA = MPI_Wtime() - tA
  chkA = checksum()
  call report('A: exchange then compute ', tA, chkA)

  ! ---- variant B ---------------------------------------------------------
  call init_field()
  call MPI_Barrier(comm, ierr)
  tB = MPI_Wtime()
  call variant_b()
  tB = MPI_Wtime() - tB
  chkB = checksum()
  call report('B: overlapped           ', tB, chkB)

  ! ---- variant C ---------------------------------------------------------
  call init_field()
  call MPI_Barrier(comm, ierr)
  tC = MPI_Wtime()
  call variant_c()
  tC = MPI_Wtime() - tC
  chkC = checksum()
  call report('C: THREAD_MULTIPLE      ', tC, chkC)

  ! ---- TODO 5: the number that goes in your notes ------------------------
  if (rank == 0) then
     print '(a)', ''
     print '(a,f10.5,a)', '  bare communication time  : ', t_comm, ' s'
     print '(a,f10.5,a)', '  A - B (time saved)       : ', tA - tB, ' s'
     ! TODO 5: overlap efficiency eta = (tA - tB) / t_comm, as a percentage.
     !         eta = 1.0 means communication is completely hidden.
     !         eta > 1.0 means something else changed too -- find out what.
     print '(a)',         '  overlap efficiency       : TODO 5'
     if (abs(chkA - chkB) > 1.0e-12_dp * abs(chkA) .or. &
         abs(chkA - chkC) > 1.0e-12_dp * abs(chkA)) then
        print '(a)', '  CHECKSUM MISMATCH -- the variants are not equivalent!'
     else
        print '(a)', '  checksums agree: PASS'
     end if
  end if

  deallocate(u, unew, multi)
  call MPI_Finalize(ierr)

contains

  ! ------------------------------------------------------------------------
  ! TODO 2: variant A — the baseline.
  !
  !   do step = 1, n_iter
  !      blocking halo exchange on u          (MPI_Sendrecv, as Exercise 01)
  !      do sweep = 1, nsweeps
  !         OpenMP-parallel stencil over the WHOLE local domain 1..n_local
  !      end do
  !   end do
  !
  ! Note the exchange is outside the sweep loop on purpose: nsweeps controls
  ! the compute-to-communication ratio, which is the knob that decides
  ! whether overlap can help at all.
  ! ------------------------------------------------------------------------
  subroutine variant_a()
    integer :: step, sweep

    do step = 1, n_iter
       ! TODO 2a: blocking halo exchange (two MPI_Sendrecv calls)
       do sweep = 1, nsweeps
          ! TODO 2b: !$omp parallel do
          call stencil(1, n_local)
          call swap()
       end do
    end do
    sweep = 0
  end subroutine variant_a

  ! ------------------------------------------------------------------------
  ! TODO 3: variant B — overlap.
  !
  !   do step = 1, n_iter
  !      post Irecv x2, Isend x2                     (no waiting)
  !      compute the INTERIOR: stencil(2, n_local-1)  <- needs no ghost cells
  !      MPI_Waitall
  !      compute the BOUNDARY: point 1 and point n_local  <- needs ghosts
  !      remaining sweeps over the whole domain
  !   end do
  !
  ! The key insight: stencil point i reads u(i-1..i+1), so only i=1 and
  ! i=n_local touch a ghost cell. Everything else can proceed immediately.
  ! ------------------------------------------------------------------------
  subroutine variant_b()
    type(MPI_Request) :: req(4)
    integer :: step, sweep, e

    do step = 1, n_iter
       ! TODO 3a: post 2 x Irecv (into u(0) and u(n_local+1)) then 2 x Isend

       ! TODO 3b: interior compute -- this runs WHILE the messages fly
       call stencil(2, n_local - 1)

       ! TODO 3c: MPI_Waitall(4, req, MPI_STATUSES_IGNORE, e)

       ! TODO 3d: boundary compute -- now the ghosts are valid
       call stencil(1, 1)
       call stencil(n_local, n_local)
       call swap()

       do sweep = 2, nsweeps
          call stencil(1, n_local)
          call swap()
       end do
    end do
    req = MPI_REQUEST_NULL; e = 0
  end subroutine variant_b

  ! ------------------------------------------------------------------------
  ! TODO 4: variant C — MPI_THREAD_MULTIPLE.
  !
  ! An ESM exchanges many fields per timestep. With THREAD_MULTIPLE each
  ! OpenMP thread can drive its own exchange concurrently:
  !
  !   !$omp parallel do
  !   do f = 1, NFIELD
  !      MPI_Sendrecv on multi(:, f)     <- legal ONLY at THREAD_MULTIPLE
  !   end do
  !
  ! Requirements: distinct tags per field (or distinct communicators — the
  ! safer choice, since tag matching across threads is a classic race).
  !
  ! Then do the same compute as variant A so the timings are comparable.
  ! If your MPI only provides FUNNELED, print a clear "skipped" message
  ! rather than producing a wrong number.
  ! ------------------------------------------------------------------------
  subroutine variant_c()
    integer :: step, sweep, f

    do step = 1, n_iter
       ! TODO 4a: parallel-over-fields exchange on multi(:, 1..NFIELD)
       do f = 1, NFIELD
          ! placeholder: serial, no MPI
       end do
       do sweep = 1, nsweeps
          call stencil(1, n_local)
          call swap()
       end do
    end do
    f = 0; sweep = 0
  end subroutine variant_c

  ! ------------------------------------------------------------------------
  ! Harness below — you should not need to change it.
  ! ------------------------------------------------------------------------

  !> Three-point smoother. Reads u, writes unew, over [lo, hi] only.
  !> Add the OpenMP directive in TODO 2b; it is intentionally absent here.
  subroutine stencil(lo, hi)
    integer, intent(in) :: lo, hi
    integer :: i
    do i = lo, hi
       unew(i) = 0.25_dp * u(i-1) + 0.5_dp * u(i) + 0.25_dp * u(i+1)
    end do
  end subroutine stencil

  subroutine swap()
    real(dp), allocatable :: tmp(:)
    call move_alloc(u, tmp)
    call move_alloc(unew, u)
    call move_alloc(tmp, unew)
  end subroutine swap

  subroutine init_field()
    integer :: i, f
    do i = 0, n_local + 1
       u(i) = sin(real(i_start + i - 1, dp) * 1.0e-5_dp)
    end do
    unew = u
    do f = 1, NFIELD
       multi(:, f) = u * real(f, dp)
    end do
  end subroutine init_field

  !> Global sum over interior points only — must be identical for A, B, C.
  real(dp) function checksum() result(s)
    real(dp) :: local
    integer  :: e
    local = sum(u(1:n_local))
    call MPI_Allreduce(local, s, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, e)
  end function checksum

  !> Time n_iter halo exchanges with no computation at all, for the
  !> denominator of the overlap-efficiency metric.
  real(dp) function time_comm_alone() result(t)
    integer :: step, e
    type(MPI_Status) :: st
    call init_field()
    call MPI_Barrier(comm, e)
    t = MPI_Wtime()
    do step = 1, n_iter
       call MPI_Sendrecv(u(n_local), 1, MPI_DOUBLE_PRECISION, right, 10, &
                         u(0),       1, MPI_DOUBLE_PRECISION, left,  10, &
                         comm, st, e)
       call MPI_Sendrecv(u(1),          1, MPI_DOUBLE_PRECISION, left,  11, &
                         u(n_local+1),  1, MPI_DOUBLE_PRECISION, right, 11, &
                         comm, st, e)
    end do
    t = MPI_Wtime() - t
    call MPI_Allreduce(MPI_IN_PLACE, t, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       comm, e)
  end function time_comm_alone

  subroutine report(label, t, chk)
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: t, chk
    real(dp) :: tmax
    integer  :: e
    call MPI_Reduce(t, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, e)
    if (rank == 0) print '(a,a,a,f10.5,a,es20.12)', '  ', label, &
         '  time ', tmax, ' s   checksum ', chk
  end subroutine report

  function thread_level_name(lvl) result(nm)
    integer, intent(in) :: lvl
    character(len=20) :: nm
    select case (lvl)
    case (MPI_THREAD_SINGLE);     nm = 'SINGLE'
    case (MPI_THREAD_FUNNELED);   nm = 'FUNNELED'
    case (MPI_THREAD_SERIALIZED); nm = 'SERIALIZED'
    case (MPI_THREAD_MULTIPLE);   nm = 'MULTIPLE'
    case default;                 nm = 'not set (TODO 1)'
    end select
  end function thread_level_name

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) n_global
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) nsweeps
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg); read(arg,*) n_iter
    end if
  end subroutine read_cli

end program hybrid

! ===========================================================================
! TODO 6 — write your answers in notes/day1.md
!
! (a) Sweep nsweeps = 1, 2, 5, 10, 50, 200 and plot (tA - tB) / tA. At which
!     compute-to-communication ratio does overlap stop paying? Explain the
!     shape of the curve.
!
! (b) Does OpenMPI make asynchronous progress? Post an Isend and then sleep
!     without calling any MPI function — does the message move? Try
!     `export OMPI_MCA_mpi_yield_when_idle=1` and an async progress thread.
!     This is why "non-blocking" does not automatically mean "overlapped".
!
! (c) Fix total cores at 8 and sweep the split: 8x1, 4x2, 2x4, 1x8
!     (ranks x threads). Which wins, and why? Relate your answer to halo
!     surface area: how does the total bytes exchanged change with rank
!     count for a fixed global problem?
!
! (d) Variant C: is MPI_THREAD_MULTIPLE faster or slower per call than
!     FUNNELED on your build? Why might an MPI library be slower at
!     MULTIPLE even when only one thread is calling?
!
! (e) DKRZ context: a coupled ICON run puts atmosphere and ocean on
!     different rank sets. Sketch how variant B's idea (compute what does
!     not depend on remote data first) applies at the COUPLING level, not
!     just the halo level. You will build this in Exercise 23.
! ===========================================================================
