! ===========================================================================
! EXERCISE 17 — GPU offload with OpenMP target
! ===========================================================================
!
! GOAL
!   Offload a 3D stencil to the GPU with OpenMP target directives, and get
!   the DATA MOVEMENT right. Three variants, timed:
!     A  naive     — map arrays in and out on every kernel launch
!     B  resident  — one target data region around the whole timeloop
!     C  selective — resident, plus `target update` only when the host
!                    genuinely needs to see the data
!
! WHY (DKRZ)
!   You already know CUDA. This exercise is not about learning GPUs — it is
!   about learning the idiom Earth System Models actually use, and about the
!   one mistake that dominates every real porting effort.
!
!   That mistake is variant A. It is not a strawman: it is what you get by
!   default when you decorate a loop nest with `!$omp target` and move on.
!   The kernel gets faster, the application gets slower, because the arrays
!   now cross PCIe (or even NVLink) twice per timestep. Ports have been
!   abandoned over this. Being able to say "your kernel is fine, your data
!   residency is not" is a large part of what a performance RSE does on a
!   GPU port.
!
!   Why directives instead of CUDA: a climate scientist has to be able to
!   read and modify the physics. Rewriting 500k lines of Fortran into CUDA C
!   would fork the codebase permanently and lock out the domain experts. So
!   ESMs use directives, keep one source, and accept some performance on the
!   table. Understand that tradeoff — it comes up in interviews.
!
! RUNNING WITHOUT A GPU
!   gfortran compiles target regions and runs them on the HOST when no
!   device is present, so you can write and debug this on macOS. The
!   program prints omp_get_num_devices(); if it says 0 you are measuring
!   host execution, and the A-vs-B difference will be small because there
!   is no real transfer. Develop here, then run on a cluster with
!   nvfortran -mp=gpu (or a GCC with -foffload=nvptx-none) for real numbers.
!
! TASKS
!   TODO 1  stencil_host      — the CPU reference
!   TODO 2  stencil_naive_gpu — target with per-call mapping (variant A)
!   TODO 3  resident data region (variant B)
!   TODO 4  target update, and a correctness trap (variant C)
!   TODO 5  tune teams/threads and the loop clauses
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - all three variants match the host reference to 1e-12 relative
!   - on a real GPU, B is several times faster than A; you can quote it
!   - you can state how many bytes cross the bus per timestep in each
!     variant, computed by hand and matching a profiler
!
! HINTS
!   - `map(to:)` copies host->device, `map(from:)` device->host,
!     `map(tofrom:)` both, `map(alloc:)` neither (device scratch).
!     The default for an array in a target region is tofrom — which is
!     exactly how variant A happens by accident.
!   - Inside a `target data` region, a `target` construct reuses the already
!     -present mapping instead of copying again. That is the entire fix.
!   - `!$omp target update to(a)` / `from(a)` moves data explicitly while
!     staying inside the data region. Use it for I/O and halo exchange.
!   - collapse(2) or collapse(3) on the loop nest: a GPU needs far more
!     parallelism than a CPU, and the outer loop alone will not fill it.
! ===========================================================================

program omp_target_stencil
  use omp_lib
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer :: n = 192, nsteps = 50
  real(dp), allocatable :: a(:,:,:), b(:,:,:)
  real(dp) :: chk_ref, chk_a, chk_b, chk_c
  real(dp) :: t_host, t_naive, t_res, t_sel
  integer  :: ndev

  call read_cli()
  ndev = omp_get_num_devices()

  print '(a)', '=== Exercise 17: OpenMP target offload ==='
  print '(a,i0,a,f8.2,a)', 'grid ', n, '^3   ', &
       real(n,dp)**3 * 8.0_dp / 1048576.0_dp, ' MiB per array'
  print '(a,i0,a,i0)', 'steps ', nsteps, '   omp_get_num_devices() = ', ndev
  if (ndev == 0) then
     print '(a)', 'NO GPU DETECTED -- target regions will run on the host.'
     print '(a)', 'Develop here; run on a cluster for meaningful numbers.'
  end if
  print '(a)', ''

  allocate(a(n,n,n), b(n,n,n))

  ! ---- host reference ----------------------------------------------------
  call init(a); b = 0.0_dp
  t_host = omp_get_wtime()
  call run_host()
  t_host = omp_get_wtime() - t_host
  chk_ref = checksum(a)
  call report('host (OpenMP CPU)     ', t_host, chk_ref, chk_ref)

  ! ---- variant A: naive mapping ------------------------------------------
  call init(a); b = 0.0_dp
  t_naive = omp_get_wtime()
  call run_naive_gpu()
  t_naive = omp_get_wtime() - t_naive
  chk_a = checksum(a)
  call report('A: map every call     ', t_naive, chk_a, chk_ref)

  ! ---- variant B: resident data ------------------------------------------
  call init(a); b = 0.0_dp
  t_res = omp_get_wtime()
  call run_resident()
  t_res = omp_get_wtime() - t_res
  chk_b = checksum(a)
  call report('B: resident data      ', t_res, chk_b, chk_ref)

  ! ---- variant C: resident + selective update ----------------------------
  call init(a); b = 0.0_dp
  t_sel = omp_get_wtime()
  call run_selective()
  t_sel = omp_get_wtime() - t_sel
  chk_c = checksum(a)
  call report('C: + target update    ', t_sel, chk_c, chk_ref)

  print '(a)', ''
  print '(a,f10.2,a)', '  bytes crossing the bus per step, variant A: ', &
       2.0_dp * real(n,dp)**3 * 8.0_dp * 2.0_dp / 1048576.0_dp, ' MiB'
  print '(a)', '  variant B: TODO 6a -- work it out and check with a profiler'
  if (t_res > 0.0_dp) &
       print '(a,f8.2,a)', '  speedup B over A: ', t_naive / t_res, 'x'

  deallocate(a, b)

contains

  ! ------------------------------------------------------------------------
  ! TODO 1: the host reference. Plain OpenMP, no offload.
  ! Swap a and b each step so the result depends on every step.
  ! ------------------------------------------------------------------------
  subroutine run_host()
    integer :: s, i, j, k
    do s = 1, nsteps
       !$omp parallel do collapse(2) private(i)
       do k = 2, n-1
          do j = 2, n-1
             do i = 2, n-1
                ! TODO 1: b(i,j,k) = 0.5*a(i,j,k) + (1/12)*(6 neighbours)
                b(i,j,k) = a(i,j,k)
             end do
          end do
       end do
       !$omp end parallel do
       call swap()
    end do
  end subroutine run_host

  ! ------------------------------------------------------------------------
  ! TODO 2: variant A — the mistake, written deliberately.
  !
  ! Put `!$omp target teams distribute parallel do collapse(2)` directly on
  ! the loop nest, with NO enclosing data region. Every launch will map a
  ! and b tofrom by default: 4 array-copies per step.
  !
  ! Write it, measure it, and keep the number. It is the "before" in every
  ! GPU-porting story you will ever tell.
  ! ------------------------------------------------------------------------
  subroutine run_naive_gpu()
    integer :: s, i, j, k
    do s = 1, nsteps
       ! TODO 2: !$omp target teams distribute parallel do collapse(2) private(i)
       do k = 2, n-1
          do j = 2, n-1
             do i = 2, n-1
                b(i,j,k) = a(i,j,k)
             end do
          end do
       end do
       call swap()
    end do
  end subroutine run_naive_gpu

  ! ------------------------------------------------------------------------
  ! TODO 3: variant B — the fix.
  !
  !   !$omp target data map(tofrom: a) map(alloc: b)
  !   do s = 1, nsteps
  !      !$omp target teams distribute parallel do collapse(2) private(i)
  !      ... same kernel ...
  !      call swap()
  !   end do
  !   !$omp end target data
  !
  ! One transfer in, one out, for the WHOLE timeloop.
  !
  ! WARNING, and this is the subtle part: `swap` uses move_alloc, which
  ! changes the host descriptors while the device mapping is keyed on the
  ! ORIGINAL addresses. Inside a target data region that is a bug waiting
  ! to happen. Either swap by copying (wasteful), or restructure to a
  ! ping-pong with an explicit parity flag and two kernels. Work out which
  ! and say why in your notes — this exact issue bites every real Fortran
  ! GPU port.
  ! ------------------------------------------------------------------------
  subroutine run_resident()
    integer :: s, i, j, k
    ! TODO 3: wrap the loop in !$omp target data ... !$omp end target data
    do s = 1, nsteps
       do k = 2, n-1
          do j = 2, n-1
             do i = 2, n-1
                b(i,j,k) = a(i,j,k)
             end do
          end do
       end do
       call swap()
    end do
  end subroutine run_resident

  ! ------------------------------------------------------------------------
  ! TODO 4: variant C — resident, plus explicit updates.
  !
  ! Real models must occasionally look at the data: to write output, to
  ! exchange halos, to compute a diagnostic. Add a `!$omp target update
  ! from(a)` every 10 steps and a host-side diagnostic (a global sum), then
  ! measure what those updates cost.
  !
  ! The trap: if you compute the diagnostic on the host WITHOUT the update,
  ! you read stale data and get a plausible wrong answer with no error
  ! message. Try it deliberately once so you recognise the symptom.
  ! ------------------------------------------------------------------------
  subroutine run_selective()
    integer :: s, i, j, k
    real(dp) :: diag
    do s = 1, nsteps
       do k = 2, n-1
          do j = 2, n-1
             do i = 2, n-1
                b(i,j,k) = a(i,j,k)
             end do
          end do
       end do
       call swap()
       if (mod(s, 10) == 0) then
          ! TODO 4: !$omp target update from(a)
          diag = sum(a)
          if (diag /= diag) print *, 'NaN'   ! keeps the sum alive
       end if
    end do
  end subroutine run_selective

  subroutine swap()
    real(dp), allocatable :: tmp(:,:,:)
    call move_alloc(a, tmp)
    call move_alloc(b, a)
    call move_alloc(tmp, b)
  end subroutine swap

  subroutine init(x)
    real(dp), intent(out) :: x(:,:,:)
    integer :: i, j, k
    !$omp parallel do collapse(2) private(i)
    do k = 1, n
       do j = 1, n
          do i = 1, n
             x(i,j,k) = real(modulo(i*7 + j*13 + k*23, 97), dp) * 0.01_dp
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine init

  real(dp) function checksum(x) result(s)
    real(dp), intent(in) :: x(:,:,:)
    s = sum(x)
  end function checksum

  subroutine report(label, t, chk, ref)
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: t, chk, ref
    real(dp) :: rel
    rel = abs(chk - ref) / max(abs(ref), 1.0e-30_dp)
    print '(a,a,a,f9.4,a,es12.4,a)', '  ', label, '  time ', t, &
         ' s   relerr ', rel, merge('   OK  ', '  FAIL ', rel < 1.0e-12_dp)
  end subroutine report

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) n
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) nsteps
    end if
  end subroutine read_cli

end program omp_target_stencil

! ===========================================================================
! TODO 6 — write your answers in notes/day5.md
!
! (a) Count the bytes crossing the bus per timestep for each variant, by
!     hand. Then confirm with nsys (nsys profile --stats=true ./omp_target)
!     on a real GPU. Do your numbers match the profiler's?
!
! (b) At n=192, one array is 54 MiB. Over PCIe gen4 (~25 GB/s effective),
!     how long does variant A's transfer take per step, and how does that
!     compare to the kernel time? Compute the arithmetic intensity of the
!     stencil and say whether the GPU can ever win with variant A.
!
! (c) The move_alloc/swap problem in TODO 3: describe exactly what goes
!     wrong, and give two fixes. Which one would you put in a code that
!     domain scientists have to maintain, and why?
!
! (d) Remove the `map(alloc: b)` clause so b defaults to tofrom. Measure.
!     Then explain why `alloc` is correct here: what does b contain at the
!     start of the region, and what does the host need from it at the end?
!
! (e) Tune it: try `num_teams`, `thread_limit`, and collapse(2) vs
!     collapse(3). What is the occupancy story? Which combination wins, and
!     is the winner the same on CPU-fallback and on a real device?
!
! (f) DKRZ context: ICON's GPU port uses OpenACC, not OpenMP target. Both
!     express the same thing. Give two practical reasons a large Fortran ESM
!     might have picked OpenACC around 2018, and say what has changed since.
!     (You write the OpenACC version next, in Exercise 18 -- answer this
!     after you have both in front of you.)
! ===========================================================================
