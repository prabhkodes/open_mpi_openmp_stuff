!> Exercises every layer: MPI, OpenMP, the Fortran module chain, and the
!> C kernel. The CMake test looks for "ALL CHECKS PASSED" in this output.
program driver
  use mpi
  use omp_lib
  use kinds,      only: wp
  use grid_mod,   only: grid_t, grid_init, grid_area
  use physics_mod, only: physics_step, physics_checksum
  implicit none

  type(grid_t) :: g
  integer  :: rank, nprocs, ierr, i
  real(wp) :: chk, total, expect
  logical  :: ok

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  call grid_init(g, 1000, 0.5_wp, 0.25_wp)
  do i = 1, 10
     call physics_step(g, 0.9_wp)
  end do
  chk = physics_checksum(g)

  call MPI_Allreduce(chk, total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                     MPI_COMM_WORLD, ierr)

  ! Every rank does identical work, so the sum is nprocs * chk. If the C
  ! kernel did not link, chk is garbage and this fails.
  expect = chk * real(nprocs, wp)
  ok = abs(total - expect) < 1.0e-9_wp * abs(expect) .and. chk /= 0.0_wp

  if (rank == 0) then
     print '(a,i0,a,i0)', 'ranks = ', nprocs, '   threads = ', omp_get_max_threads()
     print '(a,f14.4)',   'grid area      = ', grid_area(g)
     print '(a,es20.12)', 'checksum (C)   = ', chk
     print '(a,es20.12)', 'allreduce      = ', total
     if (ok) then
        print '(a)', 'ALL CHECKS PASSED'
     else
        print '(a)', 'CHECKS FAILED'
     end if
  end if

  call MPI_Finalize(ierr)
  if (.not. ok) stop 1
end program driver
