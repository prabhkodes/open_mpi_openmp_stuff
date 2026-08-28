! ===========================================================================
! EXERCISE 18 — The same kernel in OpenACC (what ICON actually uses)
! ===========================================================================
!
! GOAL
!   Port the Exercise 17 stencil to OpenACC and put the two side by side.
!   Same hardware, same maths, different directive vocabulary. Learn the
!   translation table well enough to move between them without looking it
!   up, and find the places where they genuinely differ.
!
! WHY (DKRZ)
!   ICON's GPU port is OpenACC. If you end up working on ICON's GPU
!   performance, this is the language you will be reading and writing, every
!   day. Coming from CUDA, the translation is mechanical but the mental
!   model is different: you describe the parallelism available, and the
!   compiler decides the mapping. That is a feature (portability, and
!   physicists can read it) and a curse (you must read -Minfo output to find
!   out what it actually did).
!
!   The `kernels` vs `parallel` distinction below is the single most common
!   source of "why is my OpenACC slow" — worth having a firm opinion on.
!
! TRANSLATION TABLE (learn this)
!
!   OpenMP target                        OpenACC
!   ------------------------------------ ----------------------------------
!   !$omp target data map(tofrom: a)     !$acc data copy(a)
!     map(to: a)                           copyin(a)
!     map(from: a)                         copyout(a)
!     map(alloc: a)                        create(a)
!   !$omp target teams distribute        !$acc parallel loop gang
!        parallel do                       vector
!   !$omp target update from(a)          !$acc update host(a)
!   !$omp target update to(a)            !$acc update device(a)
!   !$omp declare target                 !$acc routine
!   collapse(2)                          collapse(2)
!   is_device_ptr / use_device_ptr       host_data use_device
!
!   The one with no clean equivalent: OpenACC's `kernels` construct, which
!   asks the compiler to find the parallelism itself. Powerful when it
!   works, silently serial when it does not.
!
! RUNNING WITHOUT A GPU
!   gfortran -fopenacc compiles and runs these on the host. Verified on
!   macOS with GCC 15. Develop here, run on a cluster with nvfortran -acc
!   -Minfo=accel for real numbers and, crucially, for the compiler feedback.
!
! TASKS
!   TODO 1  acc_parallel_loop — explicit gang/vector, the workhorse
!   TODO 2  acc_kernels       — let the compiler decide; compare
!   TODO 3  data region       — copyin/copyout/create around the timeloop
!   TODO 4  acc update        — host-side diagnostics without leaving GPU
!   TODO 5  async + wait      — overlap two independent kernels
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - all variants match the host reference to 1e-12 relative
!   - you can quote whether `kernels` matched `parallel loop` on your
!     compiler, and what -Minfo said about each
!   - your OpenACC data-region timing is within ~10% of your OpenMP target
!     timing from Exercise 17 (if not, one of them is doing extra transfers)
!
! HINTS
!   - ALWAYS compile with -Minfo=accel (nvfortran) or -fopt-info-optimized
!     (gfortran). OpenACC without compiler feedback is guesswork: the
!     compiler will silently serialise a loop it cannot prove independent,
!     and you will never know from the timing alone.
!   - `gang` maps to thread blocks, `vector` to threads within a block,
!     `worker` sits between. Start with `!$acc parallel loop gang vector
!     collapse(2)` and only get more specific if -Minfo shows a problem.
!   - If the compiler reports a loop-carried dependence you know is not
!     real, `!$acc loop independent` asserts it. Be sure before you do.
!   - `!$acc async(1)` queues work; `!$acc wait(1)` synchronises. This is
!     how you overlap, and it is what Exercise 19 needs.
! ===========================================================================

program acc_port
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer :: n = 192, nsteps = 50
  real(dp), allocatable :: a(:,:,:), b(:,:,:)
  real(dp) :: chk_ref, chk
  real(dp) :: t_host, t_par, t_kern, t_data
  real(dp) :: t0

  call read_cli()

  print '(a)', '=== Exercise 18: OpenACC port ==='
  print '(a,i0,a,f8.2,a)', 'grid ', n, '^3   ', &
       real(n,dp)**3 * 8.0_dp / 1048576.0_dp, ' MiB per array'
  print '(a,i0)', 'steps ', nsteps
  print '(a)', '(build with `make nvidia` and read -Minfo=accel -- the'
  print '(a)', ' compiler feedback is the whole point of OpenACC)'
  print '(a)', ''

  allocate(a(n,n,n), b(n,n,n))

  call init(a); b = 0.0_dp
  t0 = wall(); call run_host();          t_host = wall() - t0
  chk_ref = sum(a)
  call report('host reference        ', t_host, chk_ref, chk_ref)

  call init(a); b = 0.0_dp
  t0 = wall(); call run_acc_parallel();  t_par = wall() - t0
  chk = sum(a)
  call report('acc parallel loop     ', t_par, chk, chk_ref)

  call init(a); b = 0.0_dp
  t0 = wall(); call run_acc_kernels();   t_kern = wall() - t0
  chk = sum(a)
  call report('acc kernels           ', t_kern, chk, chk_ref)

  call init(a); b = 0.0_dp
  t0 = wall(); call run_acc_data();      t_data = wall() - t0
  chk = sum(a)
  call report('acc data region       ', t_data, chk, chk_ref)

  print '(a)', ''
  print '(a)', '  Compare `acc data region` against Exercise 17 variant B.'
  print '(a)', '  They express the same thing; if they differ by more than'
  print '(a)', '  ~10%, find out which one is moving extra data.'

  deallocate(a, b)

contains

  subroutine run_host()
    integer :: s, i, j, k
    do s = 1, nsteps
       !$omp parallel do collapse(2) private(i)
       do k = 2, n-1
          do j = 2, n-1
             do i = 2, n-1
                ! TODO 1 (reference): the 7-point stencil, as Exercise 17
                b(i,j,k) = a(i,j,k)
             end do
          end do
       end do
       !$omp end parallel do
       call swap()
    end do
  end subroutine run_host

  ! ------------------------------------------------------------------------
  ! TODO 1: explicit gang/vector parallelism.
  !
  !   !$acc parallel loop gang vector collapse(2) private(i) &
  !   !$acc          copyin(a) copyout(b)
  !
  ! Note the data clauses ON the parallel construct: without an enclosing
  ! data region, this copies both arrays every step — the OpenACC spelling
  ! of Exercise 17's variant A. Write it that way first, so you have the
  ! "before" number, then fix it in TODO 3.
  ! ------------------------------------------------------------------------
  subroutine run_acc_parallel()
    integer :: s, i, j, k
    do s = 1, nsteps
       ! TODO 1: !$acc parallel loop gang vector collapse(2) private(i)
       do k = 2, n-1
          do j = 2, n-1
             do i = 2, n-1
                b(i,j,k) = a(i,j,k)
             end do
          end do
       end do
       call swap()
    end do
  end subroutine run_acc_parallel

  ! ------------------------------------------------------------------------
  ! TODO 2: `kernels` — let the compiler find the parallelism.
  !
  !   !$acc kernels
  !   do k = ... ; do j = ... ; do i = ...
  !   !$acc end kernels
  !
  ! The compiler analyses the nest and decides what is parallel. When it
  ! works it is wonderfully concise. When it cannot PROVE independence it
  ! silently generates serial device code, which is catastrophically slow
  ! and produces no warning unless you asked for -Minfo.
  !
  ! Compile with -Minfo=accel and record what it says for BOTH TODO 1 and
  ! TODO 2. That comparison is the deliverable here, more than the timing.
  ! ------------------------------------------------------------------------
  subroutine run_acc_kernels()
    integer :: s, i, j, k
    do s = 1, nsteps
       ! TODO 2: !$acc kernels ... !$acc end kernels
       do k = 2, n-1
          do j = 2, n-1
             do i = 2, n-1
                b(i,j,k) = a(i,j,k)
             end do
          end do
       end do
       call swap()
    end do
  end subroutine run_acc_kernels

  ! ------------------------------------------------------------------------
  ! TODO 3 + TODO 4: the version you would actually ship.
  !
  !   !$acc data copy(a) create(b)
  !   do s = 1, nsteps
  !      !$acc parallel loop gang vector collapse(2) private(i)
  !      ... kernel ...
  !      call swap()
  !      if (mod(s,10) == 0) then
  !         !$acc update host(a)            ! TODO 4
  !         diag = sum(a)
  !      end if
  !   end do
  !   !$acc end data
  !
  ! Same move_alloc caveat as Exercise 17 TODO 3 — the device mapping is
  ! keyed on addresses that swap() changes underneath it. Solve it the same
  ! way you solved it there, and note whether OpenACC's `present` clause
  ! makes the failure easier or harder to detect than OpenMP's did.
  !
  ! TODO 5 (stretch): add `async(1)` to the kernel and `!$acc wait(1)`
  ! before the update. On a real device, does that overlap anything here?
  ! Why not — and what would you have to restructure for it to help?
  ! ------------------------------------------------------------------------
  subroutine run_acc_data()
    integer :: s, i, j, k
    real(dp) :: diag
    ! TODO 3: !$acc data copy(a) create(b)
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
          ! TODO 4: !$acc update host(a)
          diag = sum(a)
          if (diag /= diag) print *, 'NaN'
       end if
    end do
    ! TODO 3: !$acc end data
  end subroutine run_acc_data

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

  real(dp) function wall() result(t)
    integer(8) :: c, r
    call system_clock(c, r)
    t = real(c, dp) / real(r, dp)
  end function wall

  subroutine report(label, t, chk_, ref)
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: t, chk_, ref
    real(dp) :: rel
    rel = abs(chk_ - ref) / max(abs(ref), 1.0e-30_dp)
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

end program acc_port

! ===========================================================================
! TODO 6 — write your answers in notes/day5.md
!
! (a) Paste the -Minfo=accel output for `parallel loop` and for `kernels`
!     side by side. Did the compiler parallelise both? What gang/vector
!     geometry did it choose, and did it choose the same for each?
!
! (b) Write out the OpenMP-target <-> OpenACC translation table from memory,
!     then check it against the header. Which construct has NO clean
!     equivalent in the other model, in each direction?
!
! (c) Compare your Exercise 17 variant B timing with this exercise's data
!     region. If they differ by more than 10%, profile both and find the
!     extra transfer. (Common cause: a clause defaulting to copy where you
!     meant create.)
!
! (d) `kernels` is more concise and sometimes matches `parallel loop`. Give
!     a concrete loop where it would fail to parallelise, and say how you
!     would detect that in a 500k-line model you did not write.
!
! (e) ICON's GPU port: find one ICON source file with !$acc directives (the
!     repo is public at gitlab.dkrz.de). Describe the pattern they use for
!     data residency across a timestep. Do they use one big data region, or
!     per-module `!$acc enter data`? What does that choice buy them, given
!     the model has hundreds of modules?
!
! (f) The strategic question, and a likely interview one: OpenACC is
!     effectively NVIDIA-driven, while OpenMP target is a broader standard
!     and matters for AMD (LUMI) and Intel GPUs. The posting says "adapt
!     complex models for European HPC systems" — and Europe's flagships
!     include AMD-based LUMI and NVIDIA-based JUPITER. What would you
!     recommend for a NEW component today, and how would you handle a model
!     that is already all-OpenACC? Give a concrete migration path.
! ===========================================================================
