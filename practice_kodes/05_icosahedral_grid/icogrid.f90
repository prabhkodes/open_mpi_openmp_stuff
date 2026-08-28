! ===========================================================================
! EXERCISE 05 — Build an icosahedral grid and its connectivity
! ===========================================================================
!
! GOAL
!   Construct ICON's grid from scratch: start from a 20-face icosahedron,
!   bisect every edge k times (each triangle -> 4), project the new vertices
!   onto the sphere, and build the connectivity tables a dycore actually
!   needs — cell->vertex, cell->edge, edge->cell, cell->neighbour.
!   Then verify it with Euler's formula and a spherical-area sum.
!
! WHY (DKRZ)
!   The "ICON" in ICON is *ICOsahedral Nonhydrostatic*. Every performance
!   property of the model follows from this grid:
!     - cells have no (i,j) index, so every neighbour access is INDIRECT
!       (cell_neighbour(c, 1..3)) — that is what kills vectorisation and
!       forces the nproma blocking you meet in Exercise 08
!     - the grid is nearly uniform, unlike lat-lon, so there is no polar
!       timestep restriction — this is why ESMs moved to it
!     - exactly 12 vertices have degree 5 instead of 6; those pentagon
!       points are a perennial source of special-case bugs
!   You cannot discuss ICON performance in an interview without being able
!   to describe this grid. Build it once and it is yours.
!
!   ICON names resolutions RnBk: R2B4 is ~160 km, R2B9 (~5 km) is what
!   nextGEMS-class runs use, R2B11 (~1 km) is the current frontier. Here
!   `nrefine` = k with n = 2, so nrefine=4 gives 20*4^4 = 5120 cells.
!
! TASKS
!   TODO 1  build_edges     — deduplicate the 3 edges of every triangle
!   TODO 2  refine          — bisect edges, project midpoints to the sphere
!   TODO 3  build_neighbours— cell -> up to 3 edge-adjacent cells
!   TODO 4  cell_area       — spherical excess (Girard's theorem)
!   TODO 5  grid_quality    — min/max area ratio and vertex degree histogram
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - Euler check V - E + F = 2 passes at every refinement level
!   - counts match: F = 20*4^k, E = 30*4^k, V = 10*4^k + 2
!   - total area = 4*pi to within 1e-10 (unit sphere)
!   - exactly 12 vertices have degree 5, all others degree 6
!   - every cell has exactly 3 neighbours (closed surface, no boundary)
!
! HINTS
!   - Deduplicate an edge by its sorted vertex pair (min,max). A sort of
!     3*ncells keys, or a hash map, both work. Sorting is simpler and this
!     is a setup-time cost, not an inner loop.
!   - When you bisect an edge, the midpoint must be created ONCE and shared
!     by both adjacent triangles — otherwise you get cracks and the vertex
!     count is wrong. Key the midpoint on the edge, not on the triangle.
!   - Project to the sphere by normalising: v = v / |v|. Note this does NOT
!     give equal-area cells; that is why TODO 5 asks for the ratio.
!   - Girard: area of a spherical triangle = (A + B + C - pi) * R^2, where
!     A,B,C are the interior angles. Get them from the dot products of the
!     edge tangent vectors, or use the more robust l'Huilier formula.
! ===========================================================================

module icogrid_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: PI = 3.14159265358979323846_dp

  !> Unstructured triangular grid on the unit sphere.
  type :: grid_t
     integer :: nverts = 0, nedges = 0, ncells = 0
     real(dp), allocatable :: vlon(:), vlat(:)     ! vertex coords (radians)
     real(dp), allocatable :: vxyz(:,:)            ! vxyz(3, nverts) cartesian
     integer,  allocatable :: cell_vert(:,:)       ! (3, ncells)
     integer,  allocatable :: cell_edge(:,:)       ! (3, ncells)
     integer,  allocatable :: edge_vert(:,:)       ! (2, nedges)
     integer,  allocatable :: edge_cell(:,:)       ! (2, nedges)  0 = none
     integer,  allocatable :: cell_neigh(:,:)      ! (3, ncells)  0 = none
     real(dp), allocatable :: area(:)              ! (ncells)
  end type grid_t

contains

  !> The 12 vertices and 20 faces of a regular icosahedron. This is just
  !> geometry bookkeeping — provided so you can spend your time on the
  !> refinement and connectivity, which is where the ideas are.
  subroutine base_icosahedron(g)
    type(grid_t), intent(out) :: g
    real(dp) :: phi, s
    integer :: i
    integer, parameter :: faces(3,20) = reshape([ &
         1,12, 6,   1, 6, 2,   1, 2, 8,   1, 8,11,   1,11,12, &
         2, 6,10,   6,12, 5,  12,11, 3,  11, 8, 7,   8, 2, 9, &
         4,10, 5,   4, 5, 3,   4, 3, 7,   4, 7, 9,   4, 9,10, &
         5,10, 6,   3, 5,12,   7, 3,11,   9, 7, 8,  10, 9, 2  ], [3,20])

    phi = (1.0_dp + sqrt(5.0_dp)) / 2.0_dp
    g%nverts = 12
    g%ncells = 20
    allocate(g%vxyz(3,12), g%cell_vert(3,20))

    g%vxyz(:, 1) = [-1.0_dp,  phi, 0.0_dp]
    g%vxyz(:, 2) = [ 1.0_dp,  phi, 0.0_dp]
    g%vxyz(:, 3) = [-1.0_dp, -phi, 0.0_dp]
    g%vxyz(:, 4) = [ 1.0_dp, -phi, 0.0_dp]
    g%vxyz(:, 5) = [0.0_dp, -1.0_dp,  phi]
    g%vxyz(:, 6) = [0.0_dp,  1.0_dp,  phi]
    g%vxyz(:, 7) = [0.0_dp, -1.0_dp, -phi]
    g%vxyz(:, 8) = [0.0_dp,  1.0_dp, -phi]
    g%vxyz(:, 9) = [ phi, 0.0_dp, -1.0_dp]
    g%vxyz(:,10) = [ phi, 0.0_dp,  1.0_dp]
    g%vxyz(:,11) = [-phi, 0.0_dp, -1.0_dp]
    g%vxyz(:,12) = [-phi, 0.0_dp,  1.0_dp]

    do i = 1, 12
       s = norm2(g%vxyz(:,i))
       g%vxyz(:,i) = g%vxyz(:,i) / s
    end do

    g%cell_vert = faces
  end subroutine base_icosahedron

  ! ------------------------------------------------------------------------
  ! TODO 1: build the unique edge list.
  !
  ! Each triangle contributes 3 edges; each interior edge is shared by
  ! exactly 2 triangles, so a closed triangulation has E = 3F/2 edges.
  !
  ! Fill in:
  !   g%nedges, g%edge_vert(2, nedges)     -- sorted pair (min, max)
  !   g%cell_edge(3, ncells)               -- which edge is opposite each vertex
  !   g%edge_cell(2, nedges)               -- the (up to) 2 adjacent cells
  !
  ! Suggested approach: build a 3*ncells list of (vmin, vmax, cell, slot)
  ! keys, sort it lexicographically by (vmin, vmax), then walk the sorted
  ! list — identical adjacent entries are the same edge.
  ! ------------------------------------------------------------------------
  subroutine build_edges(g)
    type(grid_t), intent(inout) :: g

    ! TODO 1: replace this stub.
    g%nedges = 0
    if (allocated(g%edge_vert)) deallocate(g%edge_vert)
    if (allocated(g%edge_cell)) deallocate(g%edge_cell)
    if (allocated(g%cell_edge)) deallocate(g%cell_edge)
    allocate(g%edge_vert(2, max(g%nedges,1)))
    allocate(g%edge_cell(2, max(g%nedges,1)))
    allocate(g%cell_edge(3, g%ncells))
    g%edge_vert = 0; g%edge_cell = 0; g%cell_edge = 0
  end subroutine build_edges

  ! ------------------------------------------------------------------------
  ! TODO 2: one level of refinement.
  !
  ! For each triangle (a,b,c):
  !   - find/create the midpoint of each edge, PROJECTED onto the sphere
  !   - replace the triangle with 4:  (a,ab,ca) (b,bc,ab) (c,ca,bc) (ab,bc,ca)
  !
  ! The critical detail: midpoints are shared between adjacent triangles.
  ! Create each one exactly once, keyed on the edge. If you run build_edges
  ! first, the edge list gives you that key for free — new vertex for edge e
  ! gets index nverts_old + e.
  !
  ! That trick makes this routine short. Use it.
  ! ------------------------------------------------------------------------
  subroutine refine(g)
    type(grid_t), intent(inout) :: g

    ! TODO 2: replace this stub with real refinement.
    ! Structure:
    !   call build_edges(g)                     ! gives you the midpoint keys
    !   nv_new = g%nverts + g%nedges
    !   nc_new = g%ncells * 4
    !   allocate new arrays, copy old vertices, append edge midpoints
    !   (normalise each midpoint onto the unit sphere!)
    !   emit 4 child triangles per parent
    !   move_alloc the new arrays into g, then call build_edges(g) again
    continue
  end subroutine refine

  ! ------------------------------------------------------------------------
  ! TODO 3: cell -> neighbour connectivity.
  !
  ! Two cells are neighbours if they share an edge. Once edge_cell is built
  ! this is a transpose: for edge e with cells (c1, c2), c2 is a neighbour
  ! of c1 and vice versa. Place the neighbour in the slot matching the
  ! shared edge, so that cell_neigh(i, c) is across cell_edge(i, c).
  !
  ! That slot alignment is not pedantry — it is what lets a flux kernel use
  ! the same loop index for the edge and the neighbour, which is the whole
  ! point of the layout.
  ! ------------------------------------------------------------------------
  subroutine build_neighbours(g)
    type(grid_t), intent(inout) :: g
    if (allocated(g%cell_neigh)) deallocate(g%cell_neigh)
    allocate(g%cell_neigh(3, g%ncells))
    g%cell_neigh = 0
    ! TODO 3: fill from g%edge_cell and g%cell_edge
  end subroutine build_neighbours

  ! ------------------------------------------------------------------------
  ! TODO 4: spherical triangle areas via Girard's theorem.
  !
  !   area = (A + B + C - pi) * R^2          (R = 1 here)
  !
  ! where A, B, C are the interior angles at the three vertices. For vertex
  ! a with neighbours b and c, the interior angle is the angle between the
  ! tangent directions a->b and a->c, which you get by projecting b and c
  ! into the tangent plane at a.
  !
  ! Watch for catastrophic cancellation: at high refinement the triangles
  ! are tiny, A+B+C is barely more than pi, and a naive acos loses most of
  ! your significant digits. l'Huilier's formula is the numerically stable
  ! alternative — try both and compare the area sum against 4*pi.
  ! ------------------------------------------------------------------------
  subroutine compute_areas(g)
    type(grid_t), intent(inout) :: g
    integer :: c
    if (allocated(g%area)) deallocate(g%area)
    allocate(g%area(g%ncells))
    g%area = 0.0_dp
    do c = 1, g%ncells
       ! TODO 4: g%area(c) = spherical_excess(v1, v2, v3)
       g%area(c) = 4.0_dp * PI / real(g%ncells, dp)   ! placeholder: uniform
    end do
  end subroutine compute_areas

  !> Angle at vertex `a` in the spherical triangle (a, b, c), all unit vectors.
  !> Used by TODO 4 — you may need to write the tangent-plane projection.
  real(dp) function sphere_angle(a, b, c) result(ang)
    real(dp), intent(in) :: a(3), b(3), c(3)
    real(dp) :: tb(3), tc(3)
    ! Project b and c into the tangent plane at a, then take the angle.
    tb = b - a * dot_product(a, b)
    tc = c - a * dot_product(a, c)
    tb = tb / max(norm2(tb), tiny(1.0_dp))
    tc = tc / max(norm2(tc), tiny(1.0_dp))
    ang = acos(max(-1.0_dp, min(1.0_dp, dot_product(tb, tc))))
  end function sphere_angle

  !> Cartesian -> geographic, for output and for Exercises 10/11.
  subroutine to_lonlat(g)
    type(grid_t), intent(inout) :: g
    integer :: i
    if (allocated(g%vlon)) deallocate(g%vlon)
    if (allocated(g%vlat)) deallocate(g%vlat)
    allocate(g%vlon(g%nverts), g%vlat(g%nverts))
    do i = 1, g%nverts
       g%vlon(i) = atan2(g%vxyz(2,i), g%vxyz(1,i))
       g%vlat(i) = asin(max(-1.0_dp, min(1.0_dp, g%vxyz(3,i))))
    end do
  end subroutine to_lonlat

end module icogrid_mod


program icogrid
  use icogrid_mod
  implicit none

  type(grid_t) :: g
  integer :: nrefine = 3
  integer :: k, npass, ntest

  call read_cli()
  npass = 0; ntest = 0

  print '(a)', '=== Exercise 05: icosahedral grid construction ==='
  print '(a,i0,a)', 'refining ', nrefine, ' times (ICON calls this R2B<k>)'
  print '(a)', ''
  print '(a)', '  level    cells    edges  verts   V-E+F   area/4pi     status'
  print '(a)', '  ---------------------------------------------------------------'

  call base_icosahedron(g)
  call build_edges(g)
  call build_neighbours(g)
  call compute_areas(g)
  call report_level(0)

  do k = 1, nrefine
     call refine(g)
     call build_neighbours(g)
     call compute_areas(g)
     call report_level(k)
  end do

  print '(a)', ''
  call grid_quality()

  print '(a)', ''
  print '(a,i0,a,i0,a)', '  ', npass, ' / ', ntest, ' checks passed'
  if (npass < ntest) then
     print '(a)', '  Work the TODOs in order: build_edges -> refine ->'
     print '(a)', '  build_neighbours -> compute_areas.'
  end if

contains

  subroutine report_level(lvl)
    integer, intent(in) :: lvl
    integer  :: euler, want_c, want_e, want_v
    real(dp) :: atot
    logical  :: ok

    euler  = g%nverts - g%nedges + g%ncells
    atot   = sum(g%area)
    want_c = 20 * 4**lvl
    want_e = 30 * 4**lvl
    want_v = 10 * 4**lvl + 2

    ok = (euler == 2) .and. (g%ncells == want_c) .and. &
         (g%nedges == want_e) .and. (g%nverts == want_v) .and. &
         (abs(atot - 4.0_dp*PI) < 1.0e-10_dp)

    ntest = ntest + 1
    if (ok) npass = npass + 1

    print '(a,i5,i9,i9,i7,i8,f12.8,a)', '  ', lvl, g%ncells, g%nedges, &
         g%nverts, euler, atot / (4.0_dp*PI), merge('     PASS', '     FAIL', ok)

    if (.not. ok) then
       print '(a,i0,a,i0,a,i0)', '        expected cells=', want_c, &
            ' edges=', want_e, ' verts=', want_v
    end if
  end subroutine report_level

  ! ------------------------------------------------------------------------
  ! TODO 5: grid quality metrics.
  !
  ! Two numbers every ICON user quotes:
  !   (a) area ratio max/min — bisect-and-project does NOT produce equal
  !       cells; the ratio grows with refinement and is why ICON applies a
  !       spring-dynamics optimisation to the raw grid.
  !   (b) vertex degree histogram — exactly 12 vertices have 5 neighbours
  !       (the original icosahedron corners); every other vertex has 6.
  !       Those 12 pentagons break the regular 6-neighbour assumption that
  !       tempting fast-path code likes to make.
  ! ------------------------------------------------------------------------
  subroutine grid_quality()
    integer, allocatable :: degree(:)
    integer :: c, i, v, n5, n6, nother
    real(dp) :: amin, amax

    print '(a)', '  --- grid quality ---'

    amin = minval(g%area); amax = maxval(g%area)
    if (amin > 0.0_dp) then
       print '(a,f10.6)', '    cell area ratio max/min : ', amax / amin
    end if
    print '(a)', '      (TODO 5: this reads 1.0 while compute_areas is a stub)'

    ! TODO 5: build the vertex degree histogram from cell_vert.
    allocate(degree(g%nverts)); degree = 0
    do c = 1, g%ncells
       do i = 1, 3
          v = g%cell_vert(i, c)
          if (v >= 1 .and. v <= g%nverts) degree(v) = degree(v) + 1
       end do
    end do
    ! Each vertex is counted once per incident triangle; for a closed
    ! triangulation the number of incident triangles == the vertex degree.
    n5 = count(degree == 5); n6 = count(degree == 6)
    nother = g%nverts - n5 - n6

    print '(a,i0)', '    vertices of degree 5    : ', n5
    print '(a,i0)', '    vertices of degree 6    : ', n6
    print '(a,i0)', '    other degrees           : ', nother

    ntest = ntest + 1
    if (n5 == 12 .and. nother == 0) then
       npass = npass + 1
       print '(a)', '    pentagon check            : PASS (exactly 12)'
    else
       print '(a)', '    pentagon check            : FAIL (expect exactly 12 degree-5)'
    end if
    deallocate(degree)
  end subroutine grid_quality

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) nrefine
    end if
  end subroutine read_cli

end program icogrid

! ===========================================================================
! TODO 6 — write your answers in notes/day2.md
!
! (a) Run with nrefine = 0..6 and tabulate cells vs mean cell edge length in
!     km (Earth radius 6371 km). Which level is closest to the ~5 km used by
!     nextGEMS-class simulations? How many cells is that globally, and how
!     much memory for ONE 3D prognostic field at 90 levels in double
!     precision?
!
! (b) What is the area ratio max/min at nrefine=6? ICON does not ship the
!     raw bisected grid — look up why "spring dynamics" grid optimisation
!     exists and what it costs you.
!
! (c) The 12 pentagon points: name two numerical problems they cause in a
!     finite-volume dycore. (Hint: think about what a 3-neighbour stencil
!     assumes, and about local truncation error.)
!
! (d) Compare this grid to a regular lat-lon grid at the same nominal
!     resolution. How many cells does lat-lon need to resolve 5 km at the
!     EQUATOR, and what is the smallest zonal cell width near the pole?
!     What does that do to the CFL-limited timestep? This is THE argument
!     for icosahedral grids — be able to give the numbers.
!
! (e) Every neighbour access here is cell_neigh(i, c) — an indirection. On
!     a structured grid it would be c+1 or c+nx. Estimate the cost
!     difference: how many extra cache lines does an indirect gather touch
!     for a 3-neighbour stencil over 1e6 cells if the neighbour indices are
!     (i) sorted and local, (ii) randomly permuted? You will measure exactly
!     this in Exercise 08.
! ===========================================================================
