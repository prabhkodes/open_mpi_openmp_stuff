! ===========================================================================
! EXERCISE 07 — Mesh partitioning: where coupled-model performance is won
! ===========================================================================
!
! GOAL
!   Partition the same mesh four ways and measure what actually matters:
!     1. linear      — split by global index (what Exercise 06 used)
!     2. blocked 2D  — recursive coordinate bisection into rectangles
!     3. Morton      — Z-order space-filling curve
!     4. Hilbert     — Hilbert space-filling curve
!   For each, report edge cut, load imbalance, neighbours per rank, and the
!   maximum halo size. Then run the Exercise 06 halo exchange on each and
!   confirm the predicted ordering with a real timing.
!
! WHY (DKRZ)
!   A coupled ESM lives or dies here. Two facts drive it:
!
!     - The halo you must exchange is proportional to the EDGE CUT. A bad
!       partition can multiply your communication volume by 10x while
!       computing exactly the same thing.
!     - Every timestep ends in a collective. A collective runs at the speed
!       of the SLOWEST rank, so a 5% load imbalance is a 5% loss on the
!       whole machine, forever.
!
!   ICON partitions the sphere with a space-filling curve for exactly these
!   reasons. When the posting says "adapt complex models for European HPC
!   systems", a large part of that is re-tuning the decomposition for a new
!   node count and topology. This exercise is that skill in miniature.
!
! TASKS
!   TODO 1  part_blocked   — recursive coordinate bisection
!   TODO 2  morton_index   — bit-interleave (Z-order) key
!   TODO 3  hilbert_index  — Hilbert curve key (harder, better locality)
!   TODO 4  edge_cut / imbalance / neighbour count metrics
!   TODO 5  answer the questions at the bottom
!
! ACCEPTANCE
!   - a table comparing all four partitioners on the same mesh
!   - Hilbert beats Morton beats blocked beats linear on edge cut
!     (if it does not, your curve is wrong — check with make picture)
!   - load imbalance under 1% for all four (they all split evenly by count;
!     the difference is purely in the CUT)
!   - you can quote the edge-cut ratio linear:Hilbert at 64 partitions
!
! HINTS
!   - Edge cut = number of (cell, neighbour) pairs whose owners differ.
!     Count each undirected edge once, or be consistent and say so.
!   - Morton: interleave the bits of (i, j). Cheap, and already much better
!     than linear, but it has long jumps at power-of-two boundaries.
!   - Hilbert: the standard iterative algorithm rotates and reflects the
!     quadrant as it descends. ~20 lines. Wikipedia's `xy2d` is the usual
!     reference implementation; port it, then verify with make picture that
!     consecutive indices are always adjacent cells.
!   - To partition by a curve: sort cells by their curve index, then hand
!     out equal-sized contiguous runs. That is the whole algorithm — the
!     quality is entirely in the curve.
! ===========================================================================

module partition_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)

contains

  !> Partition 1: contiguous blocks of global index. Fast, and bad.
  !> Each rank gets a horizontal strip nx wide, so the cut is 2*nx per rank
  !> regardless of how many cells it owns.
  subroutine part_linear(nx, ny, nparts, owner)
    integer, intent(in)  :: nx, ny, nparts
    integer, intent(out) :: owner(:)
    integer :: n, base, rem, r, lo, hi
    n = nx * ny
    base = n / nparts; rem = mod(n, nparts)
    lo = 1
    do r = 0, nparts - 1
       hi = lo + base - 1
       if (r < rem) hi = hi + 1
       owner(lo:hi) = r
       lo = hi + 1
    end do
  end subroutine part_linear

  ! ------------------------------------------------------------------------
  ! TODO 1: recursive coordinate bisection.
  !
  ! Split the cell set along its LONGEST axis into two halves of equal
  ! count, recurse until you have nparts pieces. For nparts a power of two
  ! this gives near-square blocks; for other counts, split proportionally
  ! (e.g. 6 parts -> 3 and 3, or 2 and 4 — pick and justify).
  !
  ! Near-square is the point: a square of area A has perimeter 4*sqrt(A),
  ! versus 2*nx + 2*A/nx for a strip. That perimeter IS your halo.
  ! ------------------------------------------------------------------------
  subroutine part_blocked(nx, ny, nparts, owner)
    integer, intent(in)  :: nx, ny, nparts
    integer, intent(out) :: owner(:)
    ! TODO 1: implement RCB. Suggested helper:
    !   recursive subroutine rcb(i0,i1, j0,j1, p0, np, owner)
    !     if (np == 1) { owner(cells in box) = p0; return }
    !     split the longer of (i1-i0), (j1-j0) so the two halves get
    !     np/2 and np-np/2 parts, proportional to their cell counts
    owner = 0
    if (nx < 0 .or. ny < 0 .or. nparts < 0) continue
  end subroutine part_blocked

  ! ------------------------------------------------------------------------
  ! TODO 2: Morton (Z-order) index.
  !
  ! Interleave the bits of i and j:  i = i2 i1 i0, j = j2 j1 j0
  !                                  -> j2 i2 j1 i1 j0 i0
  ! 16 bits each is plenty here; use integer(8) for the result.
  !
  ! The classic trick is the "magic number" bit-spreading:
  !   x = (x | (x << 8))  & 0x00FF00FF
  !   x = (x | (x << 4))  & 0x0F0F0F0F
  !   x = (x | (x << 2))  & 0x33333333
  !   x = (x | (x << 1))  & 0x55555555
  ! then key = spread(i) | (spread(j) << 1).
  ! In Fortran: ior, ishft, iand.
  ! ------------------------------------------------------------------------
  integer(8) function morton_index(i, j) result(key)
    integer, intent(in) :: i, j
    ! TODO 2: bit-interleave i and j.
    key = int(j, 8) * 100000_8 + int(i, 8)   ! placeholder == row-major
  end function morton_index

  ! ------------------------------------------------------------------------
  ! TODO 3: Hilbert curve index.
  !
  ! Better than Morton because it has NO long jumps: consecutive indices are
  ! always geometrically adjacent. That is exactly the property you want
  ! from a partition — it bounds the perimeter of every contiguous run.
  !
  ! Iterative form (order = number of bits, side = 2^order):
  !
  !   rx, ry, d = 0
  !   s = side/2
  !   do while (s > 0)
  !      rx = merge(1, 0, iand(x, s) > 0)
  !      ry = merge(1, 0, iand(y, s) > 0)
  !      d  = d + s * s * ieor(3 * rx, ry)
  !      call rot(side, x, y, rx, ry)       ! rotate/reflect the quadrant
  !      s = s / 2
  !   end do
  !
  ! with rot() swapping x and y and reflecting when ry == 0. Write rot too.
  ! ------------------------------------------------------------------------
  integer(8) function hilbert_index(i, j, order) result(key)
    integer, intent(in) :: i, j, order
    ! TODO 3: implement the Hilbert d2xy/xy2d transform.
    key = morton_index(i, j)   ! placeholder: falls back to Morton
    if (order < 0) continue
  end function hilbert_index

  !> Partition by sorting on a space-filling-curve key, then handing out
  !> equal contiguous runs. Shared by the Morton and Hilbert variants.
  subroutine part_by_key(nx, ny, nparts, key, owner)
    integer,    intent(in)  :: nx, ny, nparts
    integer(8), intent(in)  :: key(:)
    integer,    intent(out) :: owner(:)
    integer, allocatable :: perm(:)
    integer :: n, base, rem, r, lo, hi, k

    n = nx * ny
    allocate(perm(n))
    call argsort(key, perm)

    base = n / nparts; rem = mod(n, nparts)
    lo = 1
    do r = 0, nparts - 1
       hi = lo + base - 1
       if (r < rem) hi = hi + 1
       do k = lo, hi
          owner(perm(k)) = r
       end do
       lo = hi + 1
    end do
    deallocate(perm)
  end subroutine part_by_key

  !> Index sort (simple merge sort — O(n log n), stable, good enough here).
  subroutine argsort(a, perm)
    integer(8), intent(in)  :: a(:)
    integer,    intent(out) :: perm(:)
    integer, allocatable :: tmp(:)
    integer :: n, width, i, l, m, r
    n = size(a)
    do i = 1, n
       perm(i) = i
    end do
    allocate(tmp(n))
    width = 1
    do while (width < n)
       i = 1
       do while (i <= n)
          l = i; m = min(i + width - 1, n); r = min(i + 2*width - 1, n)
          if (m < r) call merge_run(a, perm, tmp, l, m, r)
          i = i + 2*width
       end do
       width = 2*width
    end do
    deallocate(tmp)
  end subroutine argsort

  subroutine merge_run(a, perm, tmp, l, m, r)
    integer(8), intent(in)    :: a(:)
    integer,    intent(inout) :: perm(:), tmp(:)
    integer,    intent(in)    :: l, m, r
    integer :: i, j, k
    i = l; j = m + 1; k = l
    do while (i <= m .and. j <= r)
       if (a(perm(i)) <= a(perm(j))) then
          tmp(k) = perm(i); i = i + 1
       else
          tmp(k) = perm(j); j = j + 1
       end if
       k = k + 1
    end do
    do while (i <= m)
       tmp(k) = perm(i); i = i + 1; k = k + 1
    end do
    do while (j <= r)
       tmp(k) = perm(j); j = j + 1; k = k + 1
    end do
    perm(l:r) = tmp(l:r)
  end subroutine merge_run

  ! ------------------------------------------------------------------------
  ! TODO 4: the metrics. These four numbers are the whole point.
  !
  !   edge_cut     — neighbour pairs whose owners differ (count once each)
  !   imbalance    — max_cells_per_part / mean_cells_per_part
  !   max_neigh    — the largest number of distinct peer partitions
  !   max_halo     — the largest number of remote cells any partition needs
  !
  ! max_halo is the one that predicts your message VOLUME; max_neigh
  ! predicts your message COUNT. On a modern fabric those two are limited
  ! by different hardware, so you need both.
  ! ------------------------------------------------------------------------
  subroutine metrics(nx, ny, nparts, owner, cut, imbal, max_neigh, max_halo)
    integer,  intent(in)  :: nx, ny, nparts, owner(:)
    integer,  intent(out) :: cut, max_neigh, max_halo
    real(dp), intent(out) :: imbal
    integer :: i, j, c, cn, k
    integer, allocatable :: cnt(:)
    integer :: ni(4), nj(4)

    allocate(cnt(0:nparts-1)); cnt = 0
    do c = 1, nx*ny
       cnt(owner(c)) = cnt(owner(c)) + 1
    end do
    imbal = real(maxval(cnt), dp) / (real(nx*ny, dp) / real(nparts, dp))

    ! TODO 4a: count the cut. For each cell, look at its east and north
    ! neighbour only (counting each undirected edge exactly once) and
    ! increment when the owners differ.
    cut = 0
    do j = 1, ny
       do i = 1, nx
          c = (j-1)*nx + i
          ni = [modulo(i, nx) + 1, i, 0, 0]
          nj = [j, modulo(j, ny) + 1, 0, 0]
          do k = 1, 2
             cn = (nj(k)-1)*nx + ni(k)
             ! TODO 4a: if (owner(c) /= owner(cn)) cut = cut + 1
          end do
       end do
    end do

    ! TODO 4b: max_neigh and max_halo. For each partition, collect the set
    ! of distinct remote owners it touches, and the set of distinct remote
    ! cells it needs. Report the maxima over partitions.
    max_neigh = 0
    max_halo  = 0

    deallocate(cnt)
  end subroutine metrics

end module partition_mod


program partition
  use partition_mod
  implicit none

  integer :: nx = 256, ny = 256, nparts = 16
  integer, allocatable :: owner(:)
  integer(8), allocatable :: key(:)
  integer :: i, j, c, order

  call read_cli()
  allocate(owner(nx*ny), key(nx*ny))

  order = 1
  do while (2**order < max(nx, ny))
     order = order + 1
  end do

  print '(a)', '=== Exercise 07: partition quality ==='
  print '(a,i0,a,i0,a,i0,a)', 'mesh ', nx, ' x ', ny, ' into ', nparts, ' parts'
  print '(a,i0,a)', 'ideal cells per part: ', nx*ny/nparts, &
       '   (periodic mesh, so no boundary effects)'
  print '(a)', ''
  print '(a)', '  partitioner    edge cut   ratio   imbal   max nbrs   max halo'
  print '(a)', '  ------------------------------------------------------------'

  call part_linear(nx, ny, nparts, owner)
  call show('linear     ', owner)

  call part_blocked(nx, ny, nparts, owner)
  call show('blocked RCB', owner)

  do j = 1, ny
     do i = 1, nx
        c = (j-1)*nx + i
        key(c) = morton_index(i-1, j-1)
     end do
  end do
  call part_by_key(nx, ny, nparts, key, owner)
  call show('Morton SFC ', owner)

  do j = 1, ny
     do i = 1, nx
        c = (j-1)*nx + i
        key(c) = hilbert_index(i-1, j-1, order)
     end do
  end do
  call part_by_key(nx, ny, nparts, key, owner)
  call show('Hilbert SFC', owner)

  print '(a)', ''
  print '(a)', '  Ratio is edge cut relative to the linear baseline.'
  print '(a)', '  Target: Hilbert should cut roughly sqrt(nparts) times less'
  print '(a)', '  than linear on a square mesh. Work out why before you run it.'
  print '(a)', ''
  print '(a)', '  make picture  writes partition.ppm so you can SEE the curve.'

  ! The Hilbert partition is still in `owner` here — dump it to an image.
  call write_ppm('partition.ppm', owner)
  print '(a)', '  wrote partition.ppm (Hilbert). Open it and check that every'
  print '(a)', '  colour region is COMPACT -- if a colour is scattered across'
  print '(a)', '  the image, your curve is wrong.'

  deallocate(owner, key)

contains

  !> Dump the partition as a colour image. Seeing the partition is by far
  !> the fastest way to debug a space-filling curve: a correct Hilbert
  !> partition looks like a set of blobby but connected regions, a broken
  !> one looks like confetti.
  subroutine write_ppm(fname, own)
    character(len=*), intent(in) :: fname
    integer,          intent(in) :: own(:)
    integer :: u, ii, jj, cc, p
    integer :: r, g, b
    open(newunit=u, file=fname, status='replace', action='write')
    write(u,'(a)') 'P3'
    write(u,'(i0,1x,i0)') nx, ny
    write(u,'(a)') '255'
    do jj = ny, 1, -1
       do ii = 1, nx
          cc = (jj-1)*nx + ii
          p  = own(cc)
          ! Cheap hash to spread partition ids across the colour wheel.
          r = modulo(p * 97,  256)
          g = modulo(p * 57 + 80, 256)
          b = modulo(p * 151 + 160, 256)
          write(u,'(i0,1x,i0,1x,i0)') r, g, b
       end do
    end do
    close(u)
  end subroutine write_ppm

  subroutine show(label, own)
    character(len=*), intent(in) :: label
    integer,          intent(in) :: own(:)
    integer  :: cut, mn, mh
    real(dp) :: imb
    integer, save :: baseline = 0

    call metrics(nx, ny, nparts, own, cut, imb, mn, mh)
    if (baseline == 0) baseline = max(cut, 1)

    print '(a,a,i10,f9.2,f8.3,i11,i11)', '  ', label, cut, &
         real(cut, dp) / real(baseline, dp), imb, mn, mh
  end subroutine show

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) nx
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) ny
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg); read(arg,*) nparts
    end if
  end subroutine read_cli

end program partition

! ===========================================================================
! TODO 5 — write your answers in notes/day2.md
!
! (a) Derive the edge cut for the linear partition on an nx-by-ny mesh split
!     into P horizontal strips, and for a blocked partition into sqrt(P) x
!     sqrt(P) squares. Confirm your formulas against the measured numbers.
!     What is the asymptotic ratio as P grows?
!
! (b) Run `make scaling` (nparts = 4..256). Plot edge cut vs nparts for all
!     four. Which curves have the same slope? What does the slope mean
!     physically?
!
! (c) Feed each partition into Exercise 06's halo exchange and time it. Does
!     the measured time ordering match the edge-cut ordering? If not, what
!     else is in play? (Think about message COUNT vs message VOLUME, and
!     about which ranks end up on the same node.)
!
! (d) ICON partitions a SPHERE, not a rectangle. A 2D space-filling curve
!     does not directly apply. Look up how ICON does it (hint: the
!     icosahedron's 20 faces are each a triangle that can carry its own
!     curve) and describe the approach in three sentences.
!
! (e) Everything here assumes uniform work per cell. In a real ESM it is
!     not: radiation only runs on the day side, convection only fires where
!     it is unstable, and sea ice only exists at the poles. Describe how you
!     would partition for that, and what it costs. (Search term: weighted
!     partitioning, and why ICON separates the "radiation" decomposition
!     from the dynamics one.)
!
! (f) Coupled-model version of the same question: atmosphere and ocean have
!     different cell counts and different cost per cell, and the ocean has
!     land points that do NO work at all. If you give the ocean 30% of the
!     ranks, how do you decide the split? What measurement would you take
!     first? (This is a real DKRZ interview question shape — you meet it
!     for real in Exercise 23.)
! ===========================================================================
