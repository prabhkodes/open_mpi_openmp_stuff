! ===========================================================================
! EXERCISE 02 — OpenMP in Fortran: reductions, false sharing, scheduling
! ===========================================================================
!
! GOAL
!   Four things every ESM developer trips over, measured rather than assumed:
!     (a) why a hand-rolled per-thread accumulator array is slower than
!         `reduction(+:)`  — false sharing, and how padding fixes it
!     (b) how loop schedule choice changes wall time on an imbalanced loop
!     (c) `collapse` when the outer loop is shorter than the thread count
!     (d) why a parallel floating-point sum is NOT bit-reproducible, and what
!         to do about it
!
! WHY (DKRZ)
!   (d) is the one that gets you hired. Climate models must produce
!   bit-identical results when restarted, and often when the thread count
!   changes, otherwise you cannot tell a physics bug from a rounding
!   artefact. Every ESM group has been burnt by a non-reproducible global
!   sum in a diagnostic. Being the person who can explain and fix that is
!   exactly the "further developing, testing and optimising" in the posting.
!
! TASKS
!   TODO 1  sum_reduction     — the idiomatic OpenMP reduction
!   TODO 2  sum_false_sharing — deliberately bad per-thread accumulators
!   TODO 3  sum_padded        — same, but padded to cache lines
!   TODO 4  sum_ordered       — deterministic result independent of nthreads
!   TODO 5  imbalanced_loop   — compare static / dynamic / guided
!   TODO 6  collapse_demo     — collapse(2) on a short outer loop
!   TODO 7  answer the questions at the bottom
!
! ACCEPTANCE
!   - all four sums agree with the serial reference to within 1e-9 relative
!   - sum_ordered gives a BIT-IDENTICAL result for OMP=1,2,4,8
!   - you can quote the false-sharing slowdown factor (expect 2x-20x)
!   - you can say which schedule wins on the triangular loop and why
!
! HINTS
!   - A cache line is 64 B on x86 and 128 B on Apple silicon. Padding to
!     8 doubles covers x86; check whether you need 16 on your machine.
!   - `omp_get_wtime()` not `cpu_time()` — cpu_time sums over all threads.
!   - For TODO 4: sum into per-thread partials, then add the partials in a
!     fixed order on one thread. Same order every run => same bits.
! ===========================================================================

program omp_reduce
  use omp_lib
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: PAD = 8        ! doubles per cache line (64 B / 8 B)

  integer  :: n = 50000000
  real(dp), allocatable :: a(:)
  real(dp) :: ref, s
  real(dp) :: t0, t_ref, t_red, t_false, t_pad, t_ord
  integer  :: nthreads, i

  call read_cli()
  nthreads = omp_get_max_threads()

  allocate(a(n))
  ! Values chosen so the sum is large but individual terms are small: this
  ! maximises the rounding difference between summation orders.
  !$omp parallel do
  do i = 1, n
     a(i) = 1.0_dp / real(i, dp)
  end do
  !$omp end parallel do

  print '(a)',       '=== Exercise 02: OpenMP reductions and scheduling ==='
  print '(a,i0,a,i0)', 'n = ', n, '   threads = ', nthreads
  print '(a)',       ''

  ! ---- serial reference --------------------------------------------------
  t0 = omp_get_wtime()
  ref = 0.0_dp
  do i = 1, n
     ref = ref + a(i)
  end do
  t_ref = omp_get_wtime() - t0
  print '(a,f14.10,a,f8.4,a)', '  serial          sum = ', ref, &
       '   time ', t_ref, ' s'

  ! ---- part (a): reduction vs false sharing vs padding -------------------
  t0 = omp_get_wtime(); s = sum_reduction(a);     t_red   = omp_get_wtime() - t0
  call show('reduction(+:)  ', s, t_red, ref, t_ref)

  t0 = omp_get_wtime(); s = sum_false_sharing(a); t_false = omp_get_wtime() - t0
  call show('false sharing  ', s, t_false, ref, t_ref)

  t0 = omp_get_wtime(); s = sum_padded(a);        t_pad   = omp_get_wtime() - t0
  call show('padded partials', s, t_pad, ref, t_ref)

  t0 = omp_get_wtime(); s = sum_ordered(a);       t_ord   = omp_get_wtime() - t0
  call show('ordered (repro)', s, t_ord, ref, t_ref)

  print '(a)', ''
  if (t_pad > 0.0_dp) then
     print '(a,f6.2,a)', '  false-sharing penalty: ', t_false / t_pad, &
          'x slower than padded'
  end if
  print '(a,z16)', '  ordered sum bit pattern: ', sum_ordered(a)
  print '(a)',     '  ^ re-run with OMP=1,2,4,8 -- these hex digits must not change.'

  ! ---- part (b) and (c) --------------------------------------------------
  print '(a)', ''
  call imbalanced_loop()
  print '(a)', ''
  call collapse_demo()

  deallocate(a)

contains

  ! ------------------------------------------------------------------------
  ! TODO 1: one line of OpenMP. This is the version you should almost always
  ! write. The runtime allocates a private accumulator per thread in a
  ! register and combines them at the end — no shared cache line is touched
  ! in the hot loop.
  ! ------------------------------------------------------------------------
  real(dp) function sum_reduction(x) result(total)
    real(dp), intent(in) :: x(:)
    integer :: j
    total = 0.0_dp
    ! TODO 1: add  !$omp parallel do reduction(+:total)
    do j = 1, size(x)
       total = total + x(j)
    end do
  end function sum_reduction

  ! ------------------------------------------------------------------------
  ! TODO 2: the anti-pattern. Every thread writes to partial(tid), and
  ! partial(1..nthreads) all live in ONE cache line. Each write invalidates
  ! that line in every other core's cache, so the loop degenerates into a
  ! cache-line ping-pong across the memory system.
  !
  ! Write it the naive way on purpose — you need to see the number.
  ! ------------------------------------------------------------------------
  real(dp) function sum_false_sharing(x) result(total)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: partial(:)
    integer :: j, tid, nt

    nt = omp_get_max_threads()
    allocate(partial(nt))
    partial = 0.0_dp

    ! TODO 2a: parallel region, private(tid, j), shared(partial)
    ! TODO 2b: tid = omp_get_thread_num() + 1
    ! TODO 2c: !$omp do  over j, accumulating into partial(tid)
    do j = 1, size(x)
       partial(1) = partial(1) + x(j)      ! placeholder: serial, replace me
    end do

    total = 0.0_dp
    do j = 1, nt
       total = total + partial(j)
    end do
    deallocate(partial)
    tid = 0
  end function sum_false_sharing

  ! ------------------------------------------------------------------------
  ! TODO 3: same algorithm, but give each thread its own cache line by
  ! indexing partial(1, tid) in an array dimensioned partial(PAD, nt).
  ! Compare against TODO 2. The arithmetic is identical; only the memory
  ! layout changed.
  ! ------------------------------------------------------------------------
  real(dp) function sum_padded(x) result(total)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: partial(:,:)
    integer :: j, tid, nt

    nt = omp_get_max_threads()
    allocate(partial(PAD, nt))
    partial = 0.0_dp

    ! TODO 3: as TODO 2, but accumulate into partial(1, tid)
    do j = 1, size(x)
       partial(1,1) = partial(1,1) + x(j)  ! placeholder
    end do

    total = 0.0_dp
    do j = 1, nt
       total = total + partial(1, j)
    end do
    deallocate(partial)
    tid = 0
  end function sum_padded

  ! ------------------------------------------------------------------------
  ! TODO 4: REPRODUCIBLE parallel sum.
  !
  ! `reduction(+:)` combines the per-thread partials in an unspecified order,
  ! and the number of partials depends on the thread count. Floating-point
  ! addition is not associative, so the last bits move when you change
  ! OMP_NUM_THREADS. For a climate model that means a "different" answer
  ! from a run that is physically identical.
  !
  ! Strategy: split the array into a FIXED number of chunks (independent of
  ! thread count), sum each chunk, then combine the chunk sums in index
  ! order on a single thread. Same partition + same combine order => same
  ! bits, for any number of threads.
  !
  ! Use NCHUNK below. Note this is exactly how reproducible global sums are
  ! done in production ESMs.
  ! ------------------------------------------------------------------------
  real(dp) function sum_ordered(x) result(total)
    real(dp), intent(in) :: x(:)
    integer, parameter :: NCHUNK = 256
    real(dp) :: chunk(NCHUNK)
    integer  :: c, lo, hi, m, j

    m = size(x)
    chunk = 0.0_dp

    ! TODO 4a: !$omp parallel do private(lo, hi, j) schedule(static)
    do c = 1, NCHUNK
       lo = (c - 1) * (m / NCHUNK) + 1
       hi = c * (m / NCHUNK)
       if (c == NCHUNK) hi = m           ! last chunk mops up the remainder
       ! TODO 4b: sum x(lo:hi) into chunk(c)
       do j = lo, hi
          chunk(c) = chunk(c) + x(j)
       end do
    end do

    ! Fixed serial combine order -> deterministic.
    total = 0.0_dp
    do c = 1, NCHUNK
       total = total + chunk(c)
    end do
  end function sum_ordered

  ! ------------------------------------------------------------------------
  ! TODO 5: scheduling on an imbalanced loop.
  !
  ! The inner work grows linearly with i (a triangular loop — think of a
  ! column-wise physics routine where deep columns cost more than shallow
  ! ones). With schedule(static) thread 0 gets the cheap rows and idles
  ! while the last thread grinds.
  !
  ! Fill in the three timed variants and print the ratio.
  ! ------------------------------------------------------------------------
  subroutine imbalanced_loop()
    integer, parameter :: m = 20000
    real(dp) :: out(m), t_static, t_dynamic, t_guided, t1
    integer  :: i2

    print '(a)', '  --- schedule comparison on a triangular loop ---'

    t1 = omp_get_wtime()
    ! TODO 5a: !$omp parallel do schedule(static)
    do i2 = 1, m
       out(i2) = busy_work(i2)
    end do
    t_static = omp_get_wtime() - t1

    t1 = omp_get_wtime()
    ! TODO 5b: !$omp parallel do schedule(dynamic, 16)
    do i2 = 1, m
       out(i2) = busy_work(i2)
    end do
    t_dynamic = omp_get_wtime() - t1

    t1 = omp_get_wtime()
    ! TODO 5c: !$omp parallel do schedule(guided)
    do i2 = 1, m
       out(i2) = busy_work(i2)
    end do
    t_guided = omp_get_wtime() - t1

    print '(a,f8.4,a)', '    static        ', t_static,  ' s'
    print '(a,f8.4,a)', '    dynamic,16    ', t_dynamic, ' s'
    print '(a,f8.4,a)', '    guided        ', t_guided,  ' s'
    print '(a,es14.6,a)', '    (checksum ', sum(out), ')'
  end subroutine imbalanced_loop

  !> Work proportional to i — this is what makes the loop imbalanced.
  real(dp) function busy_work(i) result(r)
    integer, intent(in) :: i
    integer :: k
    r = 0.0_dp
    do k = 1, i
       r = r + sqrt(real(k, dp))
    end do
  end function busy_work

  ! ------------------------------------------------------------------------
  ! TODO 6: collapse.
  !
  ! The outer loop has only 3 iterations (think: 3 tracer species) but the
  ! inner has 200000. With plain `parallel do` on the outer loop you can
  ! only ever use 3 threads. `collapse(2)` fuses the iteration spaces so all
  ! 3 * 200000 iterations are distributed.
  !
  ! ICON hits this constantly — short outer dimensions (levels, tracers)
  ! wrapped around long inner ones (cells).
  ! ------------------------------------------------------------------------
  subroutine collapse_demo()
    integer, parameter :: nouter = 3, ninner = 200000
    real(dp) :: f(ninner, nouter), t1, t_plain, t_coll
    integer  :: io, ii

    print '(a)', '  --- collapse(2) with a short outer loop ---'

    t1 = omp_get_wtime()
    ! TODO 6a: !$omp parallel do  (outer only -- limited to nouter threads)
    do io = 1, nouter
       do ii = 1, ninner
          f(ii, io) = sqrt(real(ii * io, dp))
       end do
    end do
    t_plain = omp_get_wtime() - t1

    t1 = omp_get_wtime()
    ! TODO 6b: !$omp parallel do collapse(2)
    do io = 1, nouter
       do ii = 1, ninner
          f(ii, io) = sqrt(real(ii * io, dp))
       end do
    end do
    t_coll = omp_get_wtime() - t1

    print '(a,f8.4,a)', '    outer only    ', t_plain, ' s'
    print '(a,f8.4,a)', '    collapse(2)   ', t_coll,  ' s'
    print '(a,es14.6,a)', '    (checksum ', sum(f), ')'
  end subroutine collapse_demo

  ! ------------------------------------------------------------------------
  subroutine show(label, val, t, reference, tser)
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: val, t, reference, tser
    real(dp) :: relerr
    relerr = abs(val - reference) / abs(reference)
    print '(a,a,a,f14.10,a,f8.4,a,f6.2,a,es9.2,a)', '  ', label, ' sum = ', &
         val, '   time ', t, ' s  speedup ', tser / max(t, 1.0e-12_dp), &
         '   relerr ', relerr, merge('  OK  ', ' FAIL ', relerr < 1.0e-9_dp)
  end subroutine show

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) n
    end if
  end subroutine read_cli

end program omp_reduce

! ===========================================================================
! TODO 7 — write your answers in notes/day1.md
!
! (a) What slowdown factor did false sharing cost you? Compute the number of
!     cache lines touched per second and compare it to your memory
!     bandwidth. Does the penalty match a coherence-traffic explanation?
!
! (b) Change PAD from 8 to 1,2,4,8,16 and plot the time. Where does the
!     curve flatten? What does that tell you about the cache line size of
!     the machine you are on?
!
! (c) Run `make repro` — sum_ordered must print identical hex for every
!     thread count. Does reduction(+:) manage that? Why not?
!
! (d) DKRZ context: an ESM must reproduce bit-identically after a restart.
!     Name three other places besides a global sum where thread or rank
!     count can silently change the last bits of a result.
!
! (e) Which schedule won the triangular loop, and what is the theoretical
!     best-case imbalance for schedule(static) on a triangular loop with
!     P threads? (Work out the ratio of the busiest to the average thread.)
! ===========================================================================
