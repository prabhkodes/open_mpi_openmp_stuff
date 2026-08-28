!> Unstructured mesh, space-filling-curve partition, and halo schedule.
!>
!> BUILD FROM: Exercise 05 (connectivity), 06 (halo schedule), 07 (SFC).
!> Most of this you can lift almost verbatim from those exercises -- that is
!> intentional. The capstone is about integration, not about re-deriving the
!> Hilbert curve at hour one.
module mesh_mod
  use kinds, only: wp
  use mpi_f08
  implicit none
  private
  public :: mesh_t, mesh_build, mesh_partition, mesh_halo_schedule, &
            mesh_exchange, mesh_stats

  type :: mesh_t
     integer :: nglobal = 0, nowned = 0, nhalo = 0
     integer, allocatable :: cell_neigh(:,:)     ! (3, nglobal) global indices
     integer, allocatable :: owner(:)            ! (nglobal)
     integer, allocatable :: owned_gidx(:), halo_gidx(:), g2l(:)
     ! communication schedule -- built once, replayed every step
     integer :: nneigh = 0
     integer, allocatable :: neigh_rank(:)
     integer, allocatable :: send_count(:), send_displ(:), send_lidx(:)
     integer, allocatable :: recv_count(:), recv_displ(:), recv_slot(:)
     real(wp), allocatable :: sendbuf(:), recvbuf(:)
  end type mesh_t

contains

  !> TODO 1: build the mesh topology. Reuse Exercise 06's generator, or
  !> Exercise 05's icosahedral one if you completed it.
  subroutine mesh_build(m, ncells)
    type(mesh_t), intent(out) :: m
    integer,      intent(in)  :: ncells
    m%nglobal = ncells
    allocate(m%cell_neigh(3, ncells))
    m%cell_neigh = 0
    ! TODO 1
  end subroutine mesh_build

  !> TODO 2: space-filling-curve partition. Reuse Exercise 07.
  !> Report the edge cut -- SPEC requires it.
  subroutine mesh_partition(m, comm)
    type(mesh_t),   intent(inout) :: m
    type(MPI_Comm), intent(in)    :: comm
    ! TODO 2
    if (m%nglobal < 0) continue
    if (comm == MPI_COMM_NULL) continue
  end subroutine mesh_partition

  !> TODO 3: discover the halo and invert it into a send schedule.
  !> Reuse Exercise 06 TODO 1 and TODO 2 -- the hard part is already solved.
  subroutine mesh_halo_schedule(m, comm)
    type(mesh_t),   intent(inout) :: m
    type(MPI_Comm), intent(in)    :: comm
    ! TODO 3
    if (m%nhalo < 0) continue
    if (comm == MPI_COMM_NULL) continue
  end subroutine mesh_halo_schedule

  !> TODO 4: pack, exchange, unpack. Reuse Exercise 06 TODO 3.
  subroutine mesh_exchange(m, f, comm)
    type(mesh_t),   intent(inout) :: m
    real(wp),       intent(inout) :: f(:)
    type(MPI_Comm), intent(in)    :: comm
    ! TODO 4
    if (size(f) < 0) continue
    if (comm == MPI_COMM_NULL) continue
  end subroutine mesh_exchange

  !> Edge cut and load imbalance -- the SPEC report table needs both.
  subroutine mesh_stats(m, comm, label)
    type(mesh_t),     intent(in) :: m
    type(MPI_Comm),   intent(in) :: comm
    character(len=*), intent(in) :: label
    integer :: rank, e, maxhalo, maxneigh
    call MPI_Comm_rank(comm, rank, e)
    call MPI_Reduce(m%nhalo,  maxhalo,  1, MPI_INTEGER, MPI_MAX, 0, comm, e)
    call MPI_Reduce(m%nneigh, maxneigh, 1, MPI_INTEGER, MPI_MAX, 0, comm, e)
    if (rank == 0) print '(a,a,a,i0,a,i0,a,i0)', '    ', label, &
         ': cells ', m%nglobal, '  max halo ', maxhalo, '  max nbrs ', maxneigh
  end subroutine mesh_stats

end module mesh_mod
