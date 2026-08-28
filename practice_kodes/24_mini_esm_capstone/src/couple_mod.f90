!> Component split, time accumulation, and the coupled exchange.
!> BUILD FROM: Exercise 09 (split), Exercise 12 (accumulator).
module couple_mod
  use kinds, only: wp
  use mpi_f08
  implicit none
  private
  public :: comp_t, accum_t, comp_init, accumulate, window_average, &
            couple_exchange

  integer, parameter, public :: COMP_ATM = 0, COMP_OCE = 1

  type :: comp_t
     integer :: id = -1
     character(len=12) :: name = ''
     type(MPI_Comm) :: comm = MPI_COMM_NULL
     integer :: rank = -1, size = 0, world_rank = -1, world_size = 0
     integer :: remote_root = -1
  end type comp_t

  type :: accum_t
     real(wp) :: sum = 0.0_wp, dt_total = 0.0_wp
  end type accum_t

contains

  !> TODO 9: split the world. Reuse Exercise 09 TODO 1.
  subroutine comp_init(c, n_atm)
    type(comp_t), intent(out) :: c
    integer,      intent(in)  :: n_atm
    integer :: ierr
    call MPI_Comm_rank(MPI_COMM_WORLD, c%world_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, c%world_size, ierr)
    if (c%world_rank < n_atm) then
       c%id = COMP_ATM; c%name = 'atmosphere'; c%remote_root = n_atm
    else
       c%id = COMP_OCE; c%name = 'ocean';      c%remote_root = 0
    end if
    ! TODO 9: MPI_Comm_split, then fill c%rank and c%size.
    c%comm = MPI_COMM_WORLD
    call MPI_Comm_rank(c%comm, c%rank, ierr)
    call MPI_Comm_size(c%comm, c%size, ierr)
  end subroutine comp_init

  subroutine accumulate(a, flux, dt)
    type(accum_t), intent(inout) :: a
    real(wp),      intent(in)    :: flux, dt
    a%sum = a%sum + flux * dt
    a%dt_total = a%dt_total + dt
  end subroutine accumulate

  function window_average(a) result(avg)
    type(accum_t), intent(inout) :: a
    real(wp) :: avg
    if (a%dt_total > 0.0_wp) then
       avg = a%sum / a%dt_total
    else
       avg = 0.0_wp
    end if
    a%sum = 0.0_wp; a%dt_total = 0.0_wp
  end function window_average

  !> TODO 10: cross-component exchange. Reuse Exercise 09 TODO 3.
  subroutine couple_exchange(c, send, recv)
    type(comp_t), intent(inout) :: c
    real(wp),     intent(in)    :: send(:)
    real(wp),     intent(out)   :: recv(:)
    recv = 0.0_wp
    ! TODO 10
    if (size(send) < 0 .or. c%size < 0) continue
  end subroutine couple_exchange

end module couple_mod
