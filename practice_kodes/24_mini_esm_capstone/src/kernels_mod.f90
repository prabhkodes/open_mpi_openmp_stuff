!> The two kernels: horizontal (indirect addressing) and column (vertical
!> dependency, and the one you offload).
!> BUILD FROM: Exercise 08 (both kernels), Exercise 17 or 18 (offload).
module kernels_mod
  use kinds,     only: wp
  use field_mod, only: field_t
  implicit none
  private
  public :: kernel_horizontal, kernel_column

  !> TODO 5: indirect 3-neighbour gather over the blocked layout.
  !> Reuse Exercise 08 TODO 3. Thread over blocks with OpenMP.
contains

  subroutine kernel_horizontal(f, out, nb_idx, nb_blk)
    type(field_t), intent(in)    :: f
    type(field_t), intent(inout) :: out
    integer,       intent(in)    :: nb_idx(:,:,:), nb_blk(:,:,:)
    ! TODO 5
    if (f%ncells < 0 .or. out%ncells < 0) continue
    if (size(nb_idx) < 0 .or. size(nb_blk) < 0) continue
  end subroutine kernel_horizontal

  !> TODO 6: column kernel with a vertical dependency, then OFFLOAD it.
  !> Reuse Exercise 08 TODO 2 for the loop, Exercise 17 TODO 3 for the
  !> resident data region. Data must stay on the device across the whole
  !> timeloop -- per-step mapping is the failure the SPEC asks you to avoid.
  subroutine kernel_column(f, g)
    type(field_t), intent(inout) :: f
    type(field_t), intent(in)    :: g
    ! TODO 6
    if (f%ncells < 0 .or. g%ncells < 0) continue
  end subroutine kernel_column

end module kernels_mod
