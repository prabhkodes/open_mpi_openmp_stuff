!> A minimal grid description. Depends on `kinds`, so kinds.f90 must be
!> compiled first -- CMake works that out from this `use` statement alone.
module grid_mod
  use kinds, only: wp
  implicit none
  private
  public :: grid_t, grid_init, grid_area

  type :: grid_t
     integer  :: ncells = 0
     real(wp) :: dx = 1.0_wp, dy = 1.0_wp
     real(wp), allocatable :: field(:)
  end type grid_t

contains

  subroutine grid_init(g, ncells, dx, dy)
    type(grid_t), intent(out) :: g
    integer,      intent(in)  :: ncells
    real(wp),     intent(in)  :: dx, dy
    integer :: i
    g%ncells = ncells
    g%dx = dx; g%dy = dy
    allocate(g%field(ncells))
    do i = 1, ncells
       g%field(i) = real(i, wp)
    end do
  end subroutine grid_init

  pure function grid_area(g) result(a)
    type(grid_t), intent(in) :: g
    real(wp) :: a
    a = real(g%ncells, wp) * g%dx * g%dy
  end function grid_area

end module grid_mod
