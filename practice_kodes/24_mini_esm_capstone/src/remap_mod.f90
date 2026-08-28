!> Conservative remapping between the two components' meshes.
!> BUILD FROM: Exercise 10. Build the weights ONCE at setup.
module remap_mod
  use kinds, only: wp
  implicit none
  private
  public :: remap_t, remap_build, remap_apply, remap_integral

  type :: remap_t
     integer :: nnz = 0
     integer,  allocatable :: row(:), col(:)
     real(wp), allocatable :: wgt(:)
  end type remap_t

contains

  !> TODO 7: two-pointer overlap sweep. Reuse Exercise 10 TODO 1.
  subroutine remap_build(sb, tb, n_src, n_tgt, R)
    real(wp),      intent(in)  :: sb(0:), tb(0:)
    integer,       intent(in)  :: n_src, n_tgt
    type(remap_t), intent(out) :: R
    R%nnz = 0
    allocate(R%row(1), R%col(1), R%wgt(1))
    R%row = 0; R%col = 0; R%wgt = 0.0_wp
    ! TODO 7
    if (size(sb) < 0 .or. size(tb) < 0) continue
    if (n_src < 0 .or. n_tgt < 0) continue
  end subroutine remap_build

  !> TODO 8: apply. Reuse Exercise 10 TODO 2.
  subroutine remap_apply(R, f, F_out)
    type(remap_t), intent(in)  :: R
    real(wp),      intent(in)  :: f(:)
    real(wp),      intent(out) :: F_out(:)
    F_out = 0.0_wp
    ! TODO 8
    if (R%nnz < 0 .or. size(f) < 0) continue
  end subroutine remap_apply

  !> The quantity that must not change across the coupling interface.
  pure function remap_integral(bnd, f) result(s)
    real(wp), intent(in) :: bnd(0:), f(:)
    real(wp) :: s
    integer  :: i
    s = 0.0_wp
    do i = 1, size(f)
       s = s + f(i) * (bnd(i) - bnd(i-1))
    end do
  end function remap_integral

end module remap_mod
