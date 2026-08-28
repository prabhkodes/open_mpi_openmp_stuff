! ===========================================================================
! EXERCISE 08 — nproma blocking and indirect addressing: ICON's data layout
! ===========================================================================
!
! GOAL
!   Measure why ICON stores every field as  f(nproma, nlev, nblks)  instead
!   of the obvious f(ncells, nlev). Sweep nproma from 1 to ncells, plot the
!   performance curve, and explain both ends of it. Then measure what
!   indirect addressing costs, and how much of that cost is recoverable by
!   renumbering the cells.
!
! WHY (DKRZ)
!   `nproma` is the first thing anyone asks you about ICON performance. The
!   idea: split ncells into nblks blocks of nproma cells, so the inner loop
!   over nproma is a clean vectorisable stride-1 loop, the block's working
!   set fits in cache across all nlev levels, and OpenMP threads over blocks.
!
!   One layout serves both machines: on CPU you want nproma ~ a few hundred
!   (cache-resident block); on GPU you want nproma = ncells (one huge
!   parallel loop). Same source code, one namelist parameter. That is why
!   ICON survived the GPU port at all — and tuning it per machine is
!   literally "adapt complex models for European HPC systems".
!
!   The second half is the cost of indirect addressing. On a structured grid
!   a neighbour is c+1. On the icosahedral grid it is neigh(c,1..3) — a
!   gather. That gather's cost depends entirely on whether neighbouring
!   cells have nearby indices, which is decided by the renumbering you chose
!   in Exercise 07. This exercise closes that loop.
!
! TASKS
!   TODO 1  column_flat      — reference kernel on f(ncells, nlev)
!   TODO 2  column_blocked   — same maths on f(nproma, nlev, nblks)
!   TODO 3  divergence_kernel— indirect 3-neighbour gather, blocked layout
!   TODO 4  nproma sweep + report GB/s and GFLOP/s
!   TODO 5  locality experiment — sorted vs randomly permuted numbering
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - blocked and flat kernels agree on the checksum to 1e-12 relative
!   - a printed nproma sweep with a clear interior maximum
!   - you can name the mechanism limiting performance at BOTH ends
!     (nproma = 1 and nproma = ncells) and they are different mechanisms
!   - you can quote the slowdown from randomly permuting the numbering
!
! HINTS
!   - Loop order in the blocked layout is  jb (blocks) -> jk (levels) ->
!     jc (cells within block). The innermost loop must be over nproma for
!     stride-1 access. Getting this backwards makes the whole thing pointless.
!   - The last block is partially full. Either carry an `npromz` for it
!     (ICON's approach) or pad and compute garbage you discard. Pick one and
!     say why in your notes — this is a real design decision with a real
!     cost, and interviewers ask about it.
!   - Thread with `!$omp parallel do` over jb. Each block is independent, so
!     there is no reduction and no false sharing — that is the third reason
!     for the layout.
!   - Compile with -O3 -march=native and check vectorisation with
!     `-fopt-info-vec` (make vecreport). If the inner loop is not vectorised,
!     find out what stopped it before you tune anything else.
! ===========================================================================

module nproma_mod
  use omp_lib
  implicit none
  integer, parameter :: dp = kind(1.0d0)

contains

  ! ------------------------------------------------------------------------
  ! TODO 1: reference kernel, flat layout f(ncells, nlev).
  !
  ! A stand-in for a physics parametrisation: a downward sweep with a
  ! vertical dependency, so level jk needs level jk-1. That dependency is
  ! what makes the COLUMN the natural unit of work in an ESM, and it is why
  ! the layout question is interesting at all.
  !
  !   f(jc, 1) stays as initialised
  !   f(jc, jk) = a*f(jc, jk-1) + b*g(jc, jk)      for jk = 2..nlev
  !
  ! Loop order: jk outer, jc inner (stride-1 in the first dimension).
  ! ------------------------------------------------------------------------
  subroutine column_flat(f, g, ncells, nlev)
    real(dp), intent(inout) :: f(:,:)      ! (ncells, nlev)
    real(dp), intent(in)    :: g(:,:)
    integer,  intent(in)    :: ncells, nlev
    real(dp), parameter :: a = 0.97_dp, b = 0.03_dp
    integer :: jc, jk

    ! TODO 1: the two nested loops. Which one is outer?
    do jk = 2, nlev
       do jc = 1, ncells
          f(jc, jk) = f(jc, jk)      ! TODO 1: a*f(jc,jk-1) + b*g(jc,jk)
       end do
    end do
    if (a < 0.0_dp .or. b < 0.0_dp) continue
  end subroutine column_flat

  ! ------------------------------------------------------------------------
  ! TODO 2: the same maths on the blocked layout f(nproma, nlev, nblks).
  !
  ! Loop nest:
  !     !$omp parallel do private(jk, jc, len)
  !     do jb = 1, nblks
  !        len = nproma;  if (jb == nblks) len = npromz     ! ragged last block
  !        do jk = 2, nlev
  !           do jc = 1, len
  !              f(jc, jk, jb) = a*f(jc, jk-1, jb) + b*g(jc, jk, jb)
  !           end do
  !        end do
  !     end do
  !
  ! Note what this buys: the inner loop is stride-1 over nproma, the whole
  ! block (nproma * nlev doubles) can sit in L2 for the entire jk sweep, and
  ! jb is an embarrassingly parallel OpenMP dimension.
  ! ------------------------------------------------------------------------
  subroutine column_blocked(f, g, nproma, nlev, nblks, npromz)
    real(dp), intent(inout) :: f(:,:,:)    ! (nproma, nlev, nblks)
    real(dp), intent(in)    :: g(:,:,:)
    integer,  intent(in)    :: nproma, nlev, nblks, npromz
    real(dp), parameter :: a = 0.97_dp, b = 0.03_dp
    integer :: jb, jk, jc, len

    ! TODO 2: add !$omp parallel do private(jk, jc, len)
    do jb = 1, nblks
       len = nproma
       if (jb == nblks) len = npromz
       do jk = 2, nlev
          do jc = 1, len
             f(jc, jk, jb) = f(jc, jk, jb)   ! TODO 2: the real update
          end do
       end do
    end do
    if (a < 0.0_dp .or. b < 0.0_dp) continue
  end subroutine column_blocked

  ! ------------------------------------------------------------------------
  ! TODO 3: the indirect-addressing kernel.
  !
  ! A divergence-like horizontal operator: each cell combines its own value
  ! with three neighbours reached through an index array. This is the shape
  ! of every dycore kernel on an unstructured grid.
  !
  !   div(jc,jk,jb) = 3*f(jc,jk,jb) - sum over n=1..3 of f(at neighbour n)
  !
  ! The neighbour arrives as a (block, index) pair, precomputed:
  !   nb_idx(jc, jb, n)  and  nb_blk(jc, jb, n)
  ! so the access is  f(nb_idx(jc,jb,n), jk, nb_blk(jc,jb,n)).
  !
  ! Look hard at that expression. Two loads just to compute the ADDRESS of
  ! the load you actually want, and the target address is unpredictable, so
  ! the hardware prefetcher is useless. That is the tax an unstructured grid
  ! charges, and it is why the renumbering in TODO 5 matters so much.
  ! ------------------------------------------------------------------------
  subroutine divergence_blocked(f, div, nb_idx, nb_blk, nproma, nlev, nblks, npromz)
    real(dp), intent(in)    :: f(:,:,:)
    real(dp), intent(out)   :: div(:,:,:)
    integer,  intent(in)    :: nb_idx(:,:,:), nb_blk(:,:,:)
    integer,  intent(in)    :: nproma, nlev, nblks, npromz
    integer :: jb, jk, jc, n, len

    ! TODO 3: add !$omp parallel do, then the real gather.
    do jb = 1, nblks
       len = nproma
       if (jb == nblks) len = npromz
       do jk = 1, nlev
          do jc = 1, len
             div(jc, jk, jb) = 3.0_dp * f(jc, jk, jb)
             ! TODO 3: subtract the three neighbour values
             do n = 1, 3
                ! div(jc,jk,jb) = div(jc,jk,jb) &
                !   - f(nb_idx(jc,jb,n), jk, nb_blk(jc,jb,n))
             end do
          end do
       end do
    end do
    if (n < 0) continue
  end subroutine divergence_blocked

  !> Flat index -> (index within block, block). ICON does this everywhere.
  pure subroutine idx_to_blk(c, nproma, jc, jb)
    integer, intent(in)  :: c, nproma
    integer, intent(out) :: jc, jb
    jb = (c - 1) / nproma + 1
    jc = c - (jb - 1) * nproma
  end subroutine idx_to_blk

end module nproma_mod


program nproma_sweep
  use nproma_mod
  use omp_lib
  implicit none

  integer :: ncells = 500000, nlev = 90
  integer :: nrep = 20
  logical :: permute = .false.

  integer, allocatable :: neigh(:,:)          ! (ncells, 3) flat neighbours
  integer, allocatable :: perm(:), iperm(:)
  real(dp) :: chk_flat
  integer  :: sweep(11) = [1, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 4096]
  integer  :: s

  call read_cli()

  print '(a)', '=== Exercise 08: nproma blocking and indirect addressing ==='
  print '(a,i0,a,i0,a,i0)', 'ncells = ', ncells, '   nlev = ', nlev, &
       '   threads = ', omp_get_max_threads()
  print '(a,f8.2,a)', 'one field = ', &
       real(ncells,dp)*real(nlev,dp)*8.0_dp/1024.0_dp**2, ' MiB'
  if (permute) then
     print '(a)', 'cell numbering: RANDOMLY PERMUTED (locality destroyed)'
  else
     print '(a)', 'cell numbering: sequential (good locality)'
  end if
  print '(a)', ''

  call build_neighbours()

  ! ---- reference: flat layout --------------------------------------------
  call run_flat(chk_flat)

  ! ---- the sweep ---------------------------------------------------------
  print '(a)', ''
  print '(a)', '  nproma   nblks    column GB/s   divergence GB/s   checksum ok'
  print '(a)', '  --------------------------------------------------------------'
  do s = 1, size(sweep)
     if (sweep(s) <= ncells) call run_blocked(sweep(s), chk_flat)
  end do
  call run_blocked(ncells, chk_flat)     ! nproma = ncells: the GPU setting

  print '(a)', ''
  print '(a)', '  Both ends of this curve are slow, for DIFFERENT reasons.'
  print '(a)', '  Name them before you read anything (TODO 6a).'
  print '(a)', ''
  print '(a)', '  Then:  make permute   -- same run, randomised cell numbering.'
  print '(a)', '  The column kernel should barely change; the divergence'
  print '(a)', '  kernel should fall off a cliff. Explain the difference.'

  deallocate(neigh)
  if (allocated(perm))  deallocate(perm)
  if (allocated(iperm)) deallocate(iperm)

contains

  !> Build a 3-neighbour topology. Sequential numbering gives neighbours at
  !> c-1, c+1 and c+stride — realistic for a well-renumbered unstructured
  !> mesh. With `permute`, the same topology is relabelled randomly, which
  !> is what a BAD renumbering looks like.
  subroutine build_neighbours()
    integer :: c, stride, i, j, tmp
    real :: r

    allocate(neigh(ncells, 3))
    stride = max(1, int(sqrt(real(ncells))))

    do c = 1, ncells
       neigh(c, 1) = modulo(c - 2,        ncells) + 1
       neigh(c, 2) = modulo(c,            ncells) + 1
       neigh(c, 3) = modulo(c + stride - 1, ncells) + 1
    end do

    if (permute) then
       allocate(perm(ncells), iperm(ncells))
       do i = 1, ncells
          perm(i) = i
       end do
       ! Fisher-Yates with a fixed seed so runs are comparable.
       call random_seed_fixed()
       do i = ncells, 2, -1
          call random_number(r)
          j = int(r * real(i)) + 1
          tmp = perm(i); perm(i) = perm(j); perm(j) = tmp
       end do
       do i = 1, ncells
          iperm(perm(i)) = i
       end do
       ! Relabel: cell perm(c) plays the role that cell c used to.
       block
         integer, allocatable :: nb2(:,:)
         allocate(nb2(ncells, 3))
         do c = 1, ncells
            do i = 1, 3
               nb2(iperm(c), i) = iperm(neigh(c, i))
            end do
         end do
         call move_alloc(nb2, neigh)
       end block
    end if
  end subroutine build_neighbours

  subroutine run_flat(chk)
    real(dp), intent(out) :: chk
    real(dp), allocatable :: f(:,:), g(:,:)
    real(dp) :: t0, t, bytes
    integer  :: rep

    allocate(f(ncells, nlev), g(ncells, nlev))
    call init2(f, g)

    t0 = omp_get_wtime()
    do rep = 1, nrep
       call column_flat(f, g, ncells, nlev)
    end do
    t = omp_get_wtime() - t0

    ! Same quantity the blocked version measures: top and bottom level of
    ! every column. Summing the whole array would not be comparable, because
    ! the blocked array has padding in its last block.
    chk = 0.0_dp
    do rep = 1, ncells
       chk = chk + f(rep, 1) + f(rep, nlev)
    end do

    bytes = 3.0_dp * real(ncells,dp) * real(nlev,dp) * 8.0_dp * real(nrep,dp)
    print '(a,f8.4,a,a,a)', '  flat f(ncells,nlev) column kernel: ', t, &
         ' s   ', rate(bytes, t), ' GB/s'
    deallocate(f, g)
  end subroutine run_flat

  !> Format a bandwidth, but refuse to print a number when the kernel took
  !> essentially no time -- that means the compiler deleted a stub kernel,
  !> and a headline "1.9e9 GB/s" helps nobody.
  function rate(bytes, t) result(s)
    real(dp), intent(in) :: bytes, t
    character(len=12) :: s
    if (t < 1.0e-5_dp) then
       s = '      stub?'
    else
       write(s, '(f12.2)') bytes / t / 1.0e9_dp
    end if
  end function rate

  subroutine run_blocked(nproma, chk_ref)
    integer,  intent(in) :: nproma
    real(dp), intent(in) :: chk_ref
    real(dp), allocatable :: f(:,:,:), g(:,:,:), div(:,:,:)
    integer,  allocatable :: nb_idx(:,:,:), nb_blk(:,:,:)
    integer  :: nblks, npromz, c, jc, jb, n, rep
    real(dp) :: t0, t_col, t_div, bytes_col, bytes_div, chk
    logical  :: ok

    nblks  = (ncells + nproma - 1) / nproma
    npromz = ncells - (nblks - 1) * nproma

    allocate(f(nproma, nlev, nblks), g(nproma, nlev, nblks))
    allocate(div(nproma, nlev, nblks))
    allocate(nb_idx(nproma, nblks, 3), nb_blk(nproma, nblks, 3))

    call init3(f, g, nproma, nlev, nblks)

    ! Translate the flat neighbour table into (index, block) pairs.
    nb_idx = 1; nb_blk = 1
    do c = 1, ncells
       call idx_to_blk(c, nproma, jc, jb)
       do n = 1, 3
          call idx_to_blk(neigh(c, n), nproma, nb_idx(jc, jb, n), nb_blk(jc, jb, n))
       end do
    end do

    t0 = omp_get_wtime()
    do rep = 1, nrep
       call column_blocked(f, g, nproma, nlev, nblks, npromz)
    end do
    t_col = omp_get_wtime() - t0

    t0 = omp_get_wtime()
    do rep = 1, nrep
       call divergence_blocked(f, div, nb_idx, nb_blk, nproma, nlev, nblks, npromz)
    end do
    t_div = omp_get_wtime() - t0

    ! Checksum over real cells only — padding in the last block is garbage.
    chk = 0.0_dp
    do c = 1, ncells
       call idx_to_blk(c, nproma, jc, jb)
       chk = chk + f(jc, 1, jb) + f(jc, nlev, jb)
    end do
    ok = abs(chk - chk_ref) <= 1.0e-12_dp * max(abs(chk_ref), 1.0_dp)

    bytes_col = 3.0_dp * real(ncells,dp) * real(nlev,dp) * 8.0_dp * real(nrep,dp)
    bytes_div = 5.0_dp * real(ncells,dp) * real(nlev,dp) * 8.0_dp * real(nrep,dp)

    print '(a,i8,i8,a,a,a,a)', '  ', nproma, nblks, &
         '   ', rate(bytes_col, t_col), &
         '      ', rate(bytes_div, t_div) // merge('        yes', '     NO (!)', ok)

    deallocate(f, g, div, nb_idx, nb_blk)
  end subroutine run_blocked

  subroutine init2(f, g)
    real(dp), intent(out) :: f(:,:), g(:,:)
    integer :: i, k
    !$omp parallel do private(i)
    do k = 1, size(f, 2)
       do i = 1, size(f, 1)
          f(i, k) = real(modulo(i * 7 + k, 100), dp) * 0.01_dp
          g(i, k) = real(modulo(i * 3 + k, 100), dp) * 0.01_dp
       end do
    end do
    !$omp end parallel do
  end subroutine init2

  subroutine init3(f, g, nproma, nlev, nblks)
    real(dp), intent(out) :: f(:,:,:), g(:,:,:)
    integer,  intent(in)  :: nproma, nlev, nblks
    integer :: c, jc, jb, k
    f = 0.0_dp; g = 0.0_dp
    do c = 1, ncells
       call idx_to_blk(c, nproma, jc, jb)
       do k = 1, nlev
          f(jc, k, jb) = real(modulo(c * 7 + k, 100), dp) * 0.01_dp
          g(jc, k, jb) = real(modulo(c * 3 + k, 100), dp) * 0.01_dp
       end do
    end do
    if (nblks < 0) continue
  end subroutine init3

  subroutine random_seed_fixed()
    integer, allocatable :: seed(:)
    integer :: sz
    call random_seed(size=sz)
    allocate(seed(sz)); seed = 20260825
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine random_seed_fixed

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) ncells
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) nlev
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg)
       permute = (trim(arg) == 'permute' .or. trim(arg) == '1')
    end if
  end subroutine read_cli

end program nproma_sweep

! ===========================================================================
! TODO 6 — write your answers in notes/day2.md
!
! (a) The sweep curve is slow at BOTH ends. Name the mechanism at each end:
!       nproma = 1        -> ?
!       nproma = ncells   -> ?
!     They are different. Where is your peak, and what is special about that
!     size on your machine? (Compute nproma * nlev * 8 bytes * how many
!     fields the kernel touches, and compare to your L2 size.)
!
! (b) Compile with -O3 -march=native -fopt-info-vec (make vecreport). Is the
!     inner loop vectorised at every nproma? At nproma = 4 on a machine with
!     8-wide AVX-512, what fraction of each vector operation is wasted?
!
! (c) Run `make permute`. Quantify: how much slower is the divergence kernel
!     with random numbering, and how much slower is the column kernel?
!     Explain why the two kernels respond so differently. Connect this back
!     to Exercise 07 — which partitioner would you want feeding this?
!
! (d) The ragged last block: this code carries npromz. The alternative is
!     padding to a full nproma and computing garbage. For ncells = 500000
!     and nproma = 4096, how many wasted cells? Now for nproma = ncells on
!     a GPU? Which approach would you choose for each machine and why?
!
! (e) Thread scaling: run OMP_NUM_THREADS = 1,2,4,8 at your best nproma.
!     Does the column kernel scale? Should it? (Careful — work out whether
!     you are compute-bound or bandwidth-bound BEFORE you look at the
!     numbers, then check whether you were right.)
!
! (f) DKRZ context: ICON runs the same source on Levante (CPU, nproma~
!     a few hundred) and on GPU partitions (nproma = number of cells in the
!     patch). Explain how one code can want two settings three orders of
!     magnitude apart, and what that implies about where the parallelism
!     comes from in each case.
! ===========================================================================
