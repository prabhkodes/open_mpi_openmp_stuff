! ===========================================================================
! EXERCISE 10 — Conservative remapping between mismatched grids
! ===========================================================================
!
! GOAL
!   Build a first-order conservative remapping operator between two grids
!   that share a domain but nothing else — different cell counts, different
!   boundaries, no nesting. Verify that it preserves the global integral to
!   machine precision. Then show what happens when you use interpolation
!   instead, by remapping back and forth a thousand times and watching the
!   mass drift away.
!
! WHY (DKRZ)
!   This is why couplers exist.
!
!   The atmosphere gives the ocean a heat flux. If the remapping is not
!   conservative, the coupled system gains or loses energy at every single
!   coupling step. A 0.1% error per step, 48 steps a day, over a 100-year
!   run, is not a rounding error — it is a fake climate trend. Reviewers of
!   climate papers ask about this specifically.
!
!   YAC's headline feature is conservative remapping on the sphere, and the
!   posting's "extending models with new components without compromising the
!   efficient execution time" means adding fields to that machinery. You
!   cannot have that conversation without having built one.
!
!   Two DIFFERENT properties, often confused. Learn to say them separately:
!     CONSERVATION — sum(F_tgt * area_tgt) == sum(f_src * area_src)
!                    (needed for fluxes: heat, freshwater, momentum)
!     CONSISTENCY  — a constant field remaps to the same constant
!                    (i.e. weights per target cell sum to 1; needed for
!                    state variables: temperature, salinity)
!   First-order conservative remapping gives you both. Most interpolation
!   schemes give you the second and not the first.
!
! TASKS
!   TODO 1  build_overlap   — the sparse weight matrix
!   TODO 2  apply_remap     — apply it
!   TODO 3  apply_linear    — a non-conservative alternative, for contrast
!   TODO 4  round_trip      — the drift experiment
!   TODO 5  answer the questions at the bottom
!
! ACCEPTANCE
!   - conservation error < 1e-15 relative for the conservative operator
!   - consistency: a constant field survives exactly (error < 1e-15)
!   - after 1000 round trips the conservative mass is unchanged to ~1e-13,
!     while the linear one has visibly drifted
!   - you can explain why first-order conservative remapping SMOOTHS the
!     field, and why that matters over a long run
!
! HINTS
!   - The overlap of source cell i = [a_i, b_i] with target cell j =
!     [c_j, d_j] is max(0, min(b_i,d_j) - max(a_i,c_j)). That one line is
!     the entire algorithm in 1D.
!   - Weight for a conservative AVERAGE onto the target:
!         w_ij = overlap_ij / area_j
!     so that F_j = sum_i w_ij * f_i. Check sum_i w_ij == 1 for every j —
!     that is the consistency test, and it will catch a domain-coverage bug
!     instantly.
!   - Do NOT do the O(n_src * n_tgt) double loop. Both boundary lists are
!     sorted, so walk them together with two pointers: O(n_src + n_tgt).
!     At ICON resolutions the naive version is genuinely unusable.
!   - Store the result as a sparse triplet list (row, col, weight). That is
!     exactly what SCRIP/YAC weight files contain.
! ===========================================================================

module remap_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  !> Sparse remap operator in triplet form: F(row) += weight * f(col).
  !> This is the same structure a SCRIP weight file stores.
  type :: remap_t
     integer :: nnz = 0
     integer,  allocatable :: row(:), col(:)
     real(dp), allocatable :: wgt(:)
  end type remap_t

contains

  !> Cell boundaries for a grid of n cells over [0,1]. `stretch` /= 0 makes
  !> the cells non-uniform, so the two grids share no boundaries at all and
  !> every target cell genuinely straddles several source cells.
  subroutine make_grid(n, stretch, bnd)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: stretch
    real(dp), allocatable, intent(out) :: bnd(:)     ! bnd(0:n)
    integer  :: i
    real(dp) :: s
    allocate(bnd(0:n))
    do i = 0, n
       s = real(i, dp) / real(n, dp)
       ! Monotone stretching that still maps [0,1] onto [0,1] exactly.
       bnd(i) = s + stretch * sin(2.0_dp * 3.14159265358979_dp * s) / 8.0_dp
    end do
    bnd(0) = 0.0_dp
    bnd(n) = 1.0_dp
  end subroutine make_grid

  pure function cell_area(bnd, i) result(a)
    real(dp), intent(in) :: bnd(0:)
    integer,  intent(in) :: i
    real(dp) :: a
    a = bnd(i) - bnd(i-1)
  end function cell_area

  ! ------------------------------------------------------------------------
  ! TODO 1: build the conservative weight matrix.
  !
  ! For every (source i, target j) pair that overlaps:
  !     ov   = max(0, min(sb(i), tb(j)) - max(sb(i-1), tb(j-1)))
  !     w_ij = ov / area_target_j
  ! and push (j, i, w_ij) onto the triplet list.
  !
  ! Use the two-pointer sweep, not the double loop:
  !     i = 1; j = 1
  !     do while (i <= n_src .and. j <= n_tgt)
  !        compute the overlap of cell i and cell j
  !        if (it is positive) record the weight
  !        advance whichever cell ENDS first
  !     end do
  !
  ! Each cell is visited once, so this is O(n_src + n_tgt) and produces the
  ! triplets already sorted by target — which is what you want for the
  ! apply step's memory access pattern.
  ! ------------------------------------------------------------------------
  subroutine build_overlap(sb, tb, n_src, n_tgt, R)
    real(dp),     intent(in)  :: sb(0:), tb(0:)
    integer,      intent(in)  :: n_src, n_tgt
    type(remap_t), intent(out) :: R
    integer :: cap

    cap = n_src + n_tgt + 4          ! two-pointer sweep cannot exceed this
    allocate(R%row(cap), R%col(cap), R%wgt(cap))
    R%nnz = 0
    R%row = 0; R%col = 0; R%wgt = 0.0_dp

    ! TODO 1: the two-pointer sweep. Leave R%nnz = 0 and every check below
    ! will fail loudly, which is the point.
    if (size(sb) < 0 .or. size(tb) < 0) continue
  end subroutine build_overlap

  ! ------------------------------------------------------------------------
  ! TODO 2: apply the operator.
  !
  !     F = 0
  !     do k = 1, nnz
  !        F(row(k)) = F(row(k)) + wgt(k) * f(col(k))
  !     end do
  !
  ! Three lines. Note the scatter into F(row(k)) — if you ever OpenMP this,
  ! that is a race unless the triplets are grouped by row (which the
  ! two-pointer sweep gives you for free). Worth saying out loud in an
  ! interview: the data structure choice at setup decided whether the hot
  ! loop is parallelisable.
  ! ------------------------------------------------------------------------
  subroutine apply_remap(R, f, F_out)
    type(remap_t), intent(in)  :: R
    real(dp),      intent(in)  :: f(:)
    real(dp),      intent(out) :: F_out(:)
    F_out = 0.0_dp
    ! TODO 2: the three-line loop.
    if (R%nnz < 0 .or. size(f) < 0) continue
  end subroutine apply_remap

  ! ------------------------------------------------------------------------
  ! TODO 3: a NON-conservative alternative, so you can measure the contrast.
  !
  ! Piecewise-linear interpolation of the source cell-centre values onto the
  ! target cell centres. Perfectly reasonable-looking. Smooth. Second-order
  ! accurate. And it does not conserve, because nothing in it knows about
  ! cell areas.
  !
  ! Find the two source centres bracketing each target centre and linearly
  ! interpolate. Clamp at the domain ends.
  ! ------------------------------------------------------------------------
  subroutine apply_linear(sb, tb, f, F_out)
    real(dp), intent(in)  :: sb(0:), tb(0:)
    real(dp), intent(in)  :: f(:)
    real(dp), intent(out) :: F_out(:)
    integer :: j
    F_out = 0.0_dp
    do j = 1, size(F_out)
       ! TODO 3: locate and interpolate at centre 0.5*(tb(j-1)+tb(j))
       F_out(j) = 0.0_dp
    end do
    if (size(sb) < 0 .or. size(f) < 0) continue
  end subroutine apply_linear

  !> Global integral of a cell-averaged field: sum(f_i * area_i).
  !> The quantity that must not change.
  pure function integral(bnd, f) result(s)
    real(dp), intent(in) :: bnd(0:), f(:)
    real(dp) :: s
    integer  :: i
    s = 0.0_dp
    do i = 1, size(f)
       s = s + f(i) * (bnd(i) - bnd(i-1))
    end do
  end function integral

end module remap_mod


program remap
  use remap_mod
  implicit none

  integer  :: n_src = 90, n_tgt = 143      ! coprime on purpose: no alignment
  integer  :: n_trip = 1000
  real(dp), allocatable :: sb(:), tb(:)
  real(dp), allocatable :: f_src(:), f_tgt(:), f_back(:)
  type(remap_t) :: S2T, T2S
  integer  :: npass, ntest, i

  call read_cli()
  npass = 0; ntest = 0

  call make_grid(n_src,  0.6_dp, sb)
  call make_grid(n_tgt, -0.4_dp, tb)

  allocate(f_src(n_src), f_tgt(n_tgt), f_back(n_src))

  print '(a)', '=== Exercise 10: conservative remapping ==='
  print '(a,i0,a,i0,a)', 'source ', n_src, ' cells  ->  target ', n_tgt, ' cells'
  print '(a,f10.6,a,f10.6)', 'src cell width min/max: ', &
       minval([(cell_area(sb,i), i=1,n_src)]), ' / ', &
       maxval([(cell_area(sb,i), i=1,n_src)])
  print '(a)', ''

  call build_overlap(sb, tb, n_src, n_tgt, S2T)
  call build_overlap(tb, sb, n_tgt, n_src, T2S)
  print '(a,i0,a,i0)', '  nonzeros src->tgt: ', S2T%nnz, '   tgt->src: ', T2S%nnz
  print '(a,i0,a)',    '  (a correct 1D sweep gives about n_src + n_tgt = ', &
       n_src + n_tgt, ')'
  print '(a)', ''

  ! ---- test 1: consistency (a constant must survive exactly) -------------
  f_src = 1.0_dp
  call apply_remap(S2T, f_src, f_tgt)
  call check('consistency: constant field', &
       maxval(abs(f_tgt - 1.0_dp)), 1.0e-15_dp)

  ! ---- test 2: conservation ----------------------------------------------
  do i = 1, n_src
     f_src(i) = 1.0_dp + sin(12.0_dp * 0.5_dp*(sb(i-1)+sb(i)))
  end do
  call apply_remap(S2T, f_src, f_tgt)
  call check('conservation: integral src=tgt', &
       abs(integral(tb, f_tgt) - integral(sb, f_src)) / &
       max(abs(integral(sb, f_src)), 1.0e-30_dp), 1.0e-14_dp)

  ! ---- test 3: the same field through the linear operator ----------------
  call apply_linear(sb, tb, f_src, f_tgt)
  print '(a,es12.4)', '  linear interp conservation error : ', &
       abs(integral(tb, f_tgt) - integral(sb, f_src)) / &
       max(abs(integral(sb, f_src)), 1.0e-30_dp)
  print '(a)', '     (this one is EXPECTED to be nonzero -- that is the lesson)'

  ! ---- test 4: the drift experiment --------------------------------------
  print '(a)', ''
  call round_trip()

  print '(a)', ''
  print '(a,i0,a,i0,a)', '  ', npass, ' / ', ntest, ' checks passed'

  deallocate(sb, tb, f_src, f_tgt, f_back)

contains

  ! ------------------------------------------------------------------------
  ! TODO 4: the drift experiment. This is the exercise's punchline.
  !
  ! Remap src -> tgt -> src, n_trip times, tracking the global integral of
  ! the source field after each round trip. Do it for BOTH operators.
  !
  ! Expected result:
  !   conservative : integral unchanged to ~1e-13 after 1000 round trips
  !                  (but the FIELD is heavily smoothed -- print min/max too)
  !   linear       : integral drifts monotonically
  !
  ! Print the integral every 100 trips so you can see the trend, and report
  ! the drift as a percentage. Then answer TODO 5b: over a 100-year run at
  ! 48 coupling steps per day, what does that percentage become?
  ! ------------------------------------------------------------------------
  subroutine round_trip()
    real(dp) :: m0, m_cons, m_lin
    integer  :: t

    print '(a,i0,a)', '  --- ', n_trip, ' round trips src->tgt->src ---'

    do i = 1, n_src
       f_src(i) = 1.0_dp + sin(12.0_dp * 0.5_dp*(sb(i-1)+sb(i)))
    end do
    m0 = integral(sb, f_src)
    print '(a,es20.12)', '    initial integral : ', m0

    ! TODO 4a: conservative round trips
    f_back = f_src
    do t = 1, n_trip
       ! call apply_remap(S2T, f_back, f_tgt)
       ! call apply_remap(T2S, f_tgt, f_back)
       ! if (mod(t, 100) == 0) print the integral
    end do
    m_cons = integral(sb, f_back)

    ! TODO 4b: linear round trips
    f_back = f_src
    do t = 1, n_trip
       ! call apply_linear(sb, tb, f_back, f_tgt)
       ! call apply_linear(tb, sb, f_tgt, f_back)
    end do
    m_lin = integral(sb, f_back)

    print '(a,es20.12,a,es10.2,a)', '    conservative     : ', m_cons, &
         '   drift ', abs(m_cons - m0)/max(abs(m0),1.0e-30_dp) * 100.0_dp, ' %'
    print '(a,es20.12,a,es10.2,a)', '    linear           : ', m_lin, &
         '   drift ', abs(m_lin  - m0)/max(abs(m0),1.0e-30_dp) * 100.0_dp, ' %'

    ! Requiring nnz > 0 matters: without it, a stubbed-out round trip leaves
    ! the field untouched and "conserves" perfectly, which would be a very
    ! misleading PASS.
    ntest = ntest + 1
    if (S2T%nnz > 0 .and. T2S%nnz > 0 .and. &
        abs(m_cons - m0) <= 1.0e-13_dp * max(abs(m0), 1.0_dp)) then
       npass = npass + 1
       print '(a)', '    round-trip conservation : PASS'
    else if (S2T%nnz == 0 .or. T2S%nnz == 0) then
       print '(a)', '    round-trip conservation : FAIL (build_overlap is still a stub)'
    else
       print '(a)', '    round-trip conservation : FAIL (mass is not conserved)'
    end if
  end subroutine round_trip

  subroutine check(label, err, tol)
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: err, tol
    logical :: ok
    ok = (err <= tol) .and. (S2T%nnz > 0)
    ntest = ntest + 1
    if (ok) npass = npass + 1
    print '(a,a32,a,es12.4,a)', '  ', label, '  err ', err, &
         merge('   PASS', '   FAIL', ok)
  end subroutine check

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) n_src
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) n_tgt
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg); read(arg,*) n_trip
    end if
  end subroutine read_cli

end program remap

! ===========================================================================
! TODO 5 — write your answers in notes/day3.md
!
! (a) State CONSERVATION and CONSISTENCY in one sentence each, and give a
!     coupled-model field that needs each one. Can a scheme be conservative
!     but not consistent? Construct one.
!
! (b) Take your measured linear-interpolation drift per round trip. A
!     coupled model exchanges fields every 30 minutes for 100 years. What is
!     the accumulated error? Express it as a fraction of the ocean's total
!     heat content and say whether a reviewer would notice.
!
! (c) First-order conservative remapping is DIFFUSIVE — print min/max of the
!     field after 1000 round trips and you will see it flattening toward the
!     mean. Explain the mechanism. Then look up second-order conservative
!     remapping (gradient reconstruction) and say what it costs and what new
!     problem it introduces. (Hint: monotonicity. What is a limiter for?)
!
! (d) Extend the overlap idea from 1D to the sphere. Two unstructured
!     spherical polygons overlap in a spherical polygon — sketch the
!     algorithm. Why is this the expensive part of coupler setup, and why is
!     it done ONCE and written to a weight file rather than every run?
!
! (e) You built the weights on one rank. In a real coupler the source field
!     lives on the atmosphere's ranks and the target on the ocean's, and no
!     single rank holds either grid. What does that do to the two-pointer
!     sweep? Sketch the parallel algorithm. (This is the hardest part of
!     writing a coupler, and it is what YAC's setup phase does.)
!
! (f) The weight matrix is fixed for the whole run. What changes if a
!     component has a MOVING grid — a regional model that follows a storm,
!     or an ice sheet whose margin advances? What would you have to rebuild,
!     and how often?
! ===========================================================================
