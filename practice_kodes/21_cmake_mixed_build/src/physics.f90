!> Depends on BOTH kinds and grid_mod, and calls into C. Two levels down the
!> dependency chain -- the file CMake must compile last of the three.
module physics_mod
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use kinds,    only: wp
  use grid_mod, only: grid_t
  implicit none
  private
  public :: physics_step, physics_checksum

  interface
     subroutine c_scale_and_sum(n, factor, x, out) bind(C, name="c_scale_and_sum")
       import :: c_int, c_double
       integer(c_int), value :: n
       real(c_double), value :: factor
       real(c_double), intent(in)  :: x(*)
       real(c_double), intent(out) :: out
     end subroutine c_scale_and_sum
  end interface

contains

  subroutine physics_step(g, factor)
    type(grid_t), intent(inout) :: g
    real(wp),     intent(in)    :: factor
    integer :: i
    !$omp parallel do
    do i = 1, g%ncells
       g%field(i) = g%field(i) * factor + 1.0_wp
    end do
    !$omp end parallel do
  end subroutine physics_step

  !> Routed through C so a broken mixed-language link fails the tests.
  function physics_checksum(g) result(s)
    type(grid_t), intent(in) :: g
    real(wp) :: s
    call c_scale_and_sum(g%ncells, 1.0_c_double, g%field, s)
  end function physics_checksum

end module physics_mod
