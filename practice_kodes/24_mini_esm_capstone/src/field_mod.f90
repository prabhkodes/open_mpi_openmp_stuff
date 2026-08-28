!> nproma-blocked field storage: f(nproma, nlev, nblks).
!> BUILD FROM: Exercise 08.
module field_mod
  use kinds, only: wp
  implicit none
  private
  public :: field_t, field_alloc, field_free, field_set, field_checksum, &
            idx_to_blk

  type :: field_t
     integer :: nproma = 0, nlev = 0, nblks = 0, npromz = 0, ncells = 0
     real(wp), allocatable :: v(:,:,:)
  end type field_t

contains

  subroutine field_alloc(f, ncells, nlev, nproma)
    type(field_t), intent(out) :: f
    integer,       intent(in)  :: ncells, nlev, nproma
    f%ncells = ncells; f%nlev = nlev; f%nproma = nproma
    f%nblks  = (ncells + nproma - 1) / nproma
    f%npromz = ncells - (f%nblks - 1) * nproma
    allocate(f%v(nproma, nlev, f%nblks))
    f%v = 0.0_wp
  end subroutine field_alloc

  subroutine field_free(f)
    type(field_t), intent(inout) :: f
    if (allocated(f%v)) deallocate(f%v)
  end subroutine field_free

  pure subroutine idx_to_blk(c, nproma, jc, jb)
    integer, intent(in)  :: c, nproma
    integer, intent(out) :: jc, jb
    jb = (c - 1) / nproma + 1
    jc = c - (jb - 1) * nproma
  end subroutine idx_to_blk

  subroutine field_set(f, c, k, val)
    type(field_t), intent(inout) :: f
    integer,       intent(in)    :: c, k
    real(wp),      intent(in)    :: val
    integer :: jc, jb
    call idx_to_blk(c, f%nproma, jc, jb)
    f%v(jc, k, jb) = val
  end subroutine field_set

  !> Sum over REAL cells only -- padding in the last block is garbage and
  !> including it makes every comparison meaningless.
  real(wp) function field_checksum(f) result(s)
    type(field_t), intent(in) :: f
    integer :: c, k, jc, jb
    s = 0.0_wp
    do c = 1, f%ncells
       call idx_to_blk(c, f%nproma, jc, jb)
       do k = 1, f%nlev
          s = s + f%v(jc, k, jb)
       end do
    end do
  end function field_checksum

end module field_mod
