! ===========================================================================
! EXERCISE 04 — Fortran/C interoperability with iso_c_binding
! ===========================================================================
!
! GOAL
!   Cross the Fortran/C boundary correctly in all five ways a coupled model
!   needs: plain arrays, multidimensional arrays (index order!), MPI
!   communicator handles, strings, and structs — plus a callback from C
!   back into Fortran.
!
! WHY (DKRZ)
!   This is not an academic exercise for this job. ICON is Fortran. YAC —
!   DKRZ's coupler, the thing you would be extending with new components —
!   is written in C and exposes a C API that ICON calls through exactly
!   these bindings:
!
!       yac_cdef_comp(comp_name, comp_id)         <- string + int
!       yac_cdef_points_unstruct(...)             <- arrays + sizes
!       yac_cget(field_id, ..., recv_field, info) <- array out
!
!   "Extending models with new components without compromising the efficient
!   execution time" means writing this glue and making sure it does not
!   copy an array it did not need to copy. The posting asks for "C/C++
!   and/or Fortran" — the honest answer for this job is both, at the seam.
!
! TASKS  (C-side tasks C1..C5 live in ckernels.c — do them together)
!   TODO 1  interface block for c_axpy, and call it
!   TODO 2  c_column_sums — get the index order right
!   TODO 3  pass MPI_COMM_WORLD to C via comm%MPI_VAL
!   TODO 4  pass a NUL-terminated string to C
!   TODO 5  bind(C) derived type matching grid_desc_t
!   TODO 6  f_scale — the Fortran routine C calls back into
!   TODO 7  answer the questions at the bottom
!
! ACCEPTANCE
!   - `make run NP=2` prints PASS for all six checks
!   - `c_sizeof_grid_desc()` equals `c_sizeof(grid_desc)` in Fortran
!   - you can explain, without looking it up, why passing a non-contiguous
!     Fortran array slice to a bind(C) routine is a trap
!
! HINTS
!   - Every interoperable interface needs `bind(C, name="...")` and dummy
!     arguments declared with c_int / c_double / c_char kinds.
!   - Scalars that C takes BY VALUE must be declared `value` in Fortran.
!     Scalars C takes by pointer must not be. Getting this backwards
!     produces a garbage value or a segfault, not a compile error.
!   - `character(kind=c_char, len=1), dimension(*)` is the portable way to
!     receive a C string; to SEND one, build it as
!     trim(name)//c_null_char and pass that.
!   - Contiguity: passing `a(1:n:2)` (a strided slice) to a bind(C) dummy
!     forces the compiler to make a temporary copy. In an ESM inner loop
!     that copy can dominate. Declare dummies `contiguous` to make the
!     compiler tell you.
! ===========================================================================

program interop
  use, intrinsic :: iso_c_binding
  use mpi_f08
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: N  = 1000

  ! -----------------------------------------------------------------------
  ! TODO 5: make this derived type interoperable with grid_desc_t in
  ! ckernels.c. Add `bind(C)` and use the c_* kinds. Member ORDER must
  ! match the C struct exactly.
  ! -----------------------------------------------------------------------
  type :: grid_desc_t          ! TODO 5: -> type, bind(C) :: grid_desc_t
     integer  :: grid_id       ! TODO 5: -> integer(c_int)
     integer  :: n_cells       ! TODO 5: -> integer(c_int)
     real(dp) :: dx            ! TODO 5: -> real(c_double)
     real(dp) :: dy            ! TODO 5: -> real(c_double)
  end type grid_desc_t

  ! -----------------------------------------------------------------------
  ! TODO 1-5: interface blocks for the C routines.
  !
  ! One is written out for you as a worked example — copy the pattern.
  ! Note `value` on the scalars: c_sizeof_grid_desc takes none, but
  ! c_axpy takes n and a BY VALUE, so both need `value`.
  ! -----------------------------------------------------------------------
  interface

     ! --- worked example: no arguments, integer return --------------------
     integer(c_int) function c_sizeof_grid_desc() bind(C, name="c_sizeof_grid_desc")
       import :: c_int
     end function c_sizeof_grid_desc

     ! --- TODO 1: void c_axpy(int n, double a, const double *x, double *y)
     ! subroutine c_axpy(n, a, x, y) bind(C, name="c_axpy")
     !   import :: c_int, c_double
     !   integer(c_int), value :: n
     !   real(c_double), value :: a
     !   real(c_double), intent(in)    :: x(*)
     !   real(c_double), intent(inout) :: y(*)
     ! end subroutine c_axpy

     ! --- TODO 2: void c_column_sums(int nrow, int ncol,
     !                                const double *m, double *colsum)

     ! --- TODO 3: double c_sum_over_comm(int fortran_comm, double local)

     ! --- TODO 4: int c_register_field(const char *name)
     !             int c_check_last_field(const char *expect)

     ! --- TODO 5: double c_grid_area(const grid_desc_t *g)
     !     (passed by reference, so NO `value` on g)

     ! --- read-only: void c_calls_fortran(int n, double factor, double *v)
     subroutine c_calls_fortran(n, factor, v) bind(C, name="c_calls_fortran")
       import :: c_int, c_double
       integer(c_int), value :: n
       real(c_double), value :: factor
       real(c_double), intent(inout) :: v(*)
     end subroutine c_calls_fortran

  end interface

  type(MPI_Comm) :: comm
  integer :: rank, nprocs, ierr
  integer :: npass, ntest

  real(c_double), allocatable, target :: x(:), y(:)
  real(c_double), allocatable, target :: m(:,:), colsum(:)
  type(grid_desc_t), target :: grid
  integer :: i, j

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)

  npass = 0; ntest = 0

  if (rank == 0) then
     print '(a)', '=== Exercise 04: Fortran/C interoperability ==='
     print '(a)', '(this is the ICON <-> YAC boundary in miniature)'
     print '(a)', ''
  end if

  ! ---- check 1: plain array kernel ---------------------------------------
  allocate(x(N), y(N))
  x = 2.0_c_double
  y = 1.0_c_double
  ! TODO 1: call c_axpy(N, 3.0_c_double, x, y)   -- expect y == 7.0 everywhere
  call check('C1 c_axpy', all(abs(y - 7.0_c_double) < 1.0e-14_c_double))

  ! ---- check 2: index order ----------------------------------------------
  ! m is 4 rows x 3 cols in FORTRAN terms; m(i,j) = i + 10*j
  allocate(m(4,3), colsum(3))
  do j = 1, 3
     do i = 1, 4
        m(i,j) = real(i + 10*j, c_double)
     end do
  end do
  colsum = 0.0_c_double
  ! TODO 2: call c_column_sums(4, 3, m, colsum)
  !         Expected column sums: (1+2+3+4) + 4*10*j = 10 + 40j
  !            j=1 -> 50,  j=2 -> 90,  j=3 -> 130
  !         If you get 33, 36, 39, ... you indexed it row-major. Fix the C.
  call check('C2 column sums (index order)', &
       abs(colsum(1) -  50.0_c_double) < 1.0e-12_c_double .and. &
       abs(colsum(2) -  90.0_c_double) < 1.0e-12_c_double .and. &
       abs(colsum(3) - 130.0_c_double) < 1.0e-12_c_double)

  ! ---- check 3: MPI communicator across the boundary ---------------------
  block
    real(c_double) :: total, expect
    total  = real(rank + 1, c_double)     ! placeholder until TODO 3
    expect = real(nprocs * (nprocs + 1) / 2, c_double)
    ! TODO 3: total = c_sum_over_comm(comm%MPI_VAL, real(rank+1, c_double))
    !
    ! comm%MPI_VAL is the Fortran integer handle inside the mpi_f08 derived
    ! type. With the older `use mpi` module the handle IS the integer, which
    ! is why so much legacy glue code passes a bare integer around.
    call check('C3 MPI_Comm_f2c across boundary', &
         abs(total - expect) < 1.0e-12_c_double)
  end block

  ! ---- check 4: strings --------------------------------------------------
  block
    character(len=*), parameter :: fname = 'sea_surface_temperature'
    integer :: slen, ok
    slen = -1; ok = 0
    ! TODO 4a: slen = c_register_field(fname//c_null_char)
    ! TODO 4b: ok   = c_check_last_field(fname//c_null_char)
    call check('C4 string crossing (len + content)', &
         slen == len(fname) .and. ok == 1)
  end block

  ! ---- check 5: struct interop -------------------------------------------
  block
    real(c_double) :: area
    logical :: size_ok
    grid%grid_id = 7
    grid%n_cells = 100
    grid%dx      = 0.5_c_double
    grid%dy      = 0.25_c_double
    area = 0.0_c_double
    ! TODO 5: area = c_grid_area(grid)      -- expect 100 * 0.5 * 0.25 = 12.5

    ! storage_size works on ANY type, so this compiles before you add
    ! bind(C). Once the type IS interoperable, switch to the proper tool:
    !     size_ok = (c_sizeof(grid) == int(c_sizeof_grid_desc(), c_size_t))
    ! c_sizeof refuses to accept a non-interoperable type at all, which is
    ! a compile-time guarantee worth having in real glue code.
    size_ok = (storage_size(grid) / 8 == c_sizeof_grid_desc())
    if (rank == 0 .and. .not. size_ok) then
       print '(a,i0,a,i0)', '     struct size mismatch: Fortran ', &
            storage_size(grid) / 8, ' vs C ', c_sizeof_grid_desc()
    end if
    call check('C5 struct interop (area + sizeof)', &
         size_ok .and. abs(area - 12.5_c_double) < 1.0e-12_c_double)
  end block

  ! ---- check 6: C calls back into Fortran --------------------------------
  block
    real(c_double), target :: v(5)
    v = 2.0_c_double
    if (rank == 0) call c_calls_fortran(5, 4.0_c_double, v)
    call MPI_Bcast(v, 5, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call check('C6 Fortran callback from C', &
         all(abs(v - 8.0_c_double) < 1.0e-14_c_double))
  end block

  ! ---- summary -----------------------------------------------------------
  if (rank == 0) then
     print '(a)', ''
     print '(a,i0,a,i0,a)', '  ', npass, ' / ', ntest, ' checks passed'
     if (npass < ntest) print '(a)', &
          '  Work through the TODOs in interop.f90 AND ckernels.c together.'
  end if

  deallocate(x, y, m, colsum)
  call MPI_Finalize(ierr)

contains

  subroutine check(label, ok)
    character(len=*), intent(in) :: label
    logical,          intent(in) :: ok
    ntest = ntest + 1
    if (ok) npass = npass + 1
    if (rank == 0) print '(a,a34,a)', '  ', label, &
         merge('  PASS', '  FAIL', ok)
  end subroutine check

end program interop

! ===========================================================================
! TODO 6 — the Fortran routine that C calls back into.
!
! It must be a MODULE procedure or an external subroutine (NOT a contained
! procedure of the program — those are not interoperable), declared
! bind(C, name="f_scale") so the C linker can find the symbol.
!
! Signature on the C side:   void f_scale(int n, double factor, double *v)
!
! Uncomment and complete:
! ---------------------------------------------------------------------------
! subroutine f_scale(n, factor, v) bind(C, name="f_scale")
!   use, intrinsic :: iso_c_binding
!   implicit none
!   integer(c_int), value :: n
!   real(c_double), value :: factor
!   real(c_double), intent(inout) :: v(*)
!   integer :: i
!   do i = 1, n
!      v(i) = v(i) * factor
!   end do
! end subroutine f_scale
! ===========================================================================

! Stub so the skeleton links before you do TODO 6. Delete it once you write
! the real one above (leaving both is a duplicate-symbol link error).
subroutine f_scale(n, factor, v) bind(C, name="f_scale")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), value :: n
  real(c_double), value :: factor
  real(c_double), intent(inout) :: v(*)
  integer :: i
  do i = 1, n
     v(i) = v(i)          ! TODO 6: multiply by factor
  end do
  i = int(factor)
end subroutine f_scale

! ===========================================================================
! TODO 7 — write your answers in notes/day1.md
!
! (a) Change c_column_sums to index m[i*ncol + j] (row-major). It still
!     compiles and still runs. What do the sums become, and why is a bug
!     that produces plausible-looking numbers worse than one that crashes?
!
! (b) Declare the c_axpy dummy `real(c_double), intent(in), contiguous :: x(:)`
!     and pass a strided slice x(1:N:2). What does gfortran do — copy, or
!     refuse? Use -Wall and check. Why does this matter for a coupler that
!     receives a field defined on every third grid cell?
!
! (c) MPI handles: with `use mpi` (old module) a communicator IS an integer,
!     with `use mpi_f08` it is a derived type. Which one can C receive
!     directly, and what does MPI_Comm_f2c actually do on OpenMPI? (Look at
!     the header if you want — it is a lookup, not a cast.)
!
! (d) Add a `character(len=32) :: name` member to grid_desc_t and the
!     matching `char name[32]` to the C struct. Does c_sizeof still agree?
!     What if you put it between the two ints and the two doubles?
!
! (e) DKRZ context: YAC's C API takes grid definitions as flat arrays of
!     cell-to-vertex connectivity. ICON stores that as a Fortran array
!     dimensioned (nproma, nblks, 3). Sketch what has to happen at the
!     boundary, and say whether it can be zero-copy. (You will meet this
!     layout for real in Exercise 08.)
! ===========================================================================
