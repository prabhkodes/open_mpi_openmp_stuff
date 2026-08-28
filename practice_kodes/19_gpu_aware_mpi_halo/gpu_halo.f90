! ===========================================================================
! EXERCISE 19 — GPU-aware MPI: halo exchange without a round trip to the host
! ===========================================================================
!
! GOAL
!   Exchange halos between MPI ranks whose data lives on the GPU. Three
!   ways, timed:
!     A  staged   — copy device->host, MPI, copy host->device (the obvious way)
!     B  direct   — pass DEVICE pointers straight to MPI (GPU-aware MPI)
!     C  overlap  — direct, plus computing the interior while it flies
!
! WHY (DKRZ)
!   This is where a GPU port of a distributed model succeeds or fails.
!
!   The kernels are the easy part. Once the dycore runs on the GPU, the
!   halo exchange becomes the bottleneck, and the naive version does four
!   bus crossings per exchange (down, over, back, and the same in reverse).
!   On a machine where each node has 4 GPUs and the network is attached to
!   the CPU, that traffic can cost more than the computation it enables.
!
!   GPU-aware MPI lets the library take a device pointer directly, and on
!   the right hardware the data goes GPU->NIC->GPU without touching host
!   memory at all (GPUDirect RDMA). The application-side change is small;
!   knowing whether your MPI actually supports it, and proving it does, is
!   the real skill.
!
!   Combine with Exercise 03's overlap idea and you have the complete
!   picture: compute the interior on the GPU while the halo is in flight.
!   That is the state of the art for an ESM dycore, and it is variant C.
!
! RUNNING WITHOUT A GPU
!   gfortran runs the OpenACC regions on the host, so all three variants
!   work and produce correct answers — but A and B will be nearly identical,
!   because there is no bus to cross. Develop and verify correctness here;
!   the performance story only appears on real hardware.
!
!   Checking for GPU-aware MPI on a cluster:
!       ompi_info --parsable --all | grep mpi_built_with_cuda_support
!       # or for UCX:
!       ucx_info -d | grep -i cuda
!   With OpenMPI you often also need:  export OMPI_MCA_opal_cuda_support=1
!
! TASKS
!   TODO 1  exchange_staged     — device->host->MPI->host->device
!   TODO 2  exchange_direct     — host_data use_device, MPI on device pointers
!   TODO 3  exchange_overlap    — interior compute during the exchange
!   TODO 4  measure the bus traffic in each case
!   TODO 5  answer the questions at the bottom
!
! ACCEPTANCE
!   - all three produce identical halo contents (exact global-index check)
!   - on a GPU cluster: B beats A, and you can quote the factor
!   - you can state whether your MPI is GPU-aware and how you PROVED it
!     (not "the docs say so" -- how you demonstrated it)
!
! HINTS
!   - The magic clause is OpenACC's:
!         !$acc host_data use_device(sendbuf, recvbuf)
!            call MPI_Isend(sendbuf, ...)
!         !$acc end host_data
!     Inside that block, the Fortran array name resolves to the DEVICE
!     address. The OpenMP equivalent is
!         !$omp target data use_device_ptr(sendbuf)
!   - If your MPI is NOT GPU-aware, passing a device pointer usually
!     segfaults rather than failing gracefully. Test on a tiny message first.
!   - Pack the halo into a contiguous device buffer with a kernel; do not
!     rely on MPI to gather a strided device array. Strided device access
!     from the NIC is either unsupported or very slow.
!   - For variant C, remember the ordering from Exercise 03: post receives,
!     pack and send, compute interior, wait, compute boundary.
! ===========================================================================

program gpu_halo
  use mpi_f08
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer :: nx = 512, ny = 512, nsteps = 100
  type(MPI_Comm) :: comm
  integer :: rank, nprocs, ierr, up, down
  real(dp), allocatable :: u(:,:), unew(:,:)
  real(dp), allocatable :: sbuf_up(:), sbuf_dn(:), rbuf_up(:), rbuf_dn(:)
  real(dp) :: t_staged, t_direct, t_overlap
  logical  :: ok_s, ok_d, ok_o

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)
  call read_cli()

  ! 1D decomposition in y; each rank owns ny rows with one halo row each side.
  up   = modulo(rank + 1, nprocs)
  down = modulo(rank - 1, nprocs)

  allocate(u(nx, 0:ny+1), unew(nx, 0:ny+1))
  allocate(sbuf_up(nx), sbuf_dn(nx), rbuf_up(nx), rbuf_dn(nx))

  if (rank == 0) then
     print '(a)', '=== Exercise 19: GPU-aware MPI halo exchange ==='
     print '(a,i0,a,i0,a,i0)', 'local tile ', nx, ' x ', ny, '   ranks ', nprocs
     print '(a,f8.2,a)', 'halo message size: ', real(nx,dp)*8.0_dp/1024.0_dp, ' KiB'
     print '(a,i0)', 'steps ', nsteps
     print '(a)', ''
  end if

  call run_variant('A: staged via host    ', 1, t_staged,  ok_s)
  call run_variant('B: direct device ptr  ', 2, t_direct,  ok_d)
  call run_variant('C: direct + overlap   ', 3, t_overlap, ok_o)

  if (rank == 0) then
     print '(a)', ''
     print '(a,f10.2,a)', '  bus bytes per step, variant A: ', &
          4.0_dp * real(nx,dp) * 8.0_dp / 1024.0_dp, ' KiB (4 crossings)'
     print '(a)', '  variant B: TODO 4 -- how many crossings, and why?'
     if (t_direct > 0.0_dp) print '(a,f8.2,a)', &
          '  speedup B over A : ', t_staged / t_direct, 'x'
     if (t_overlap > 0.0_dp) print '(a,f8.2,a)', &
          '  speedup C over A : ', t_staged / t_overlap, 'x'
     print '(a)', ''
     print '(a)', '  On a host-only build these will be nearly equal --'
     print '(a)', '  there is no bus to avoid. Run on a GPU cluster.'
  end if

  deallocate(u, unew, sbuf_up, sbuf_dn, rbuf_up, rbuf_dn)
  call MPI_Finalize(ierr)

contains

  subroutine run_variant(label, which, t, ok)
    character(len=*), intent(in)  :: label
    integer,          intent(in)  :: which
    real(dp),         intent(out) :: t
    logical,          intent(out) :: ok
    integer :: s, e

    call init_field()

    ! !$acc data copy(u) create(unew, sbuf_up, sbuf_dn, rbuf_up, rbuf_dn)
    ! TODO: wrap the whole timeloop in a data region -- exactly as in
    ! Exercise 18. Without it you are measuring transfers, not exchange.

    call MPI_Barrier(comm, e)
    t = MPI_Wtime()
    do s = 1, nsteps
       select case (which)
       case (1); call exchange_staged()
       case (2); call exchange_direct()
       case (3); call exchange_overlap()
       end select
    end do
    t = MPI_Wtime() - t
    ! !$acc end data

    ok = halo_is_correct()
    call finish_report(label, t, ok)
  end subroutine run_variant

  ! ------------------------------------------------------------------------
  ! TODO 1: variant A — stage through the host.
  !
  !   1. pack the edge rows into sbuf_* with an !$acc parallel loop
  !   2. !$acc update host(sbuf_up, sbuf_dn)        <- bus crossing 1 and 2
  !   3. MPI_Isend/Irecv on the HOST buffers, MPI_Waitall
  !   4. !$acc update device(rbuf_up, rbuf_dn)      <- crossing 3 and 4
  !   5. unpack into the halo rows with an !$acc parallel loop
  !
  ! This is what you get if you port the kernels and leave the communication
  ! alone. It is correct, it is simple, and on a real machine it is slow.
  ! ------------------------------------------------------------------------
  subroutine exchange_staged()
    type(MPI_Request) :: req(4)
    integer :: e
    ! TODO 1
    req = MPI_REQUEST_NULL; e = 0
  end subroutine exchange_staged

  ! ------------------------------------------------------------------------
  ! TODO 2: variant B — hand MPI the device pointer.
  !
  !   1. pack into sbuf_* on the device (same kernel as variant A)
  !   2. !$acc host_data use_device(sbuf_up, sbuf_dn, rbuf_up, rbuf_dn)
  !         MPI_Irecv / MPI_Isend / MPI_Waitall
  !      !$acc end host_data
  !   3. unpack on the device
  !
  ! No `update` calls at all. On hardware with GPUDirect RDMA the bytes go
  ! straight from GPU memory to the NIC.
  !
  ! If this segfaults, your MPI is not GPU-aware. That is a legitimate
  ! finding — record it and how you diagnosed it.
  ! ------------------------------------------------------------------------
  subroutine exchange_direct()
    type(MPI_Request) :: req(4)
    integer :: e
    ! TODO 2
    req = MPI_REQUEST_NULL; e = 0
  end subroutine exchange_direct

  ! ------------------------------------------------------------------------
  ! TODO 3: variant C — overlap the exchange with interior compute.
  !
  !   1. pack and post everything (as variant B), using !$acc async(1)
  !   2. run the stencil over the INTERIOR rows 2..ny-1 -- these need no
  !      halo, so they can proceed immediately
  !   3. MPI_Waitall, unpack
  !   4. run the stencil over rows 1 and ny
  !
  ! Same idea as Exercise 03, one level up the stack. On a GPU there is an
  ! extra subtlety: the interior kernel and the pack kernel are on different
  ! streams, so you need async() and wait() to express which may run
  ! concurrently. Get that wrong and you serialise anyway, silently.
  ! ------------------------------------------------------------------------
  subroutine exchange_overlap()
    type(MPI_Request) :: req(4)
    integer :: e
    ! TODO 3
    req = MPI_REQUEST_NULL; e = 0
  end subroutine exchange_overlap

  ! ---- harness -----------------------------------------------------------

  !> Interior rows carry their global row index, so a correct halo has an
  !> exactly predictable value and the check is not statistical.
  subroutine init_field()
    integer :: i, j
    u = -1.0_dp
    do j = 1, ny
       do i = 1, nx
          u(i,j) = real(rank * ny + j, dp)
       end do
    end do
    unew = u
  end subroutine init_field

  logical function halo_is_correct() result(ok)
    integer  :: i
    real(dp) :: want_dn, want_up
    ok = .true.
    want_dn = real(modulo(rank*ny - 1, nprocs*ny) + 1, dp)
    want_up = real(modulo(rank*ny + ny, nprocs*ny) + 1, dp)
    do i = 1, nx
       if (u(i, 0)    /= want_dn) ok = .false.
       if (u(i, ny+1) /= want_up) ok = .false.
    end do
  end function halo_is_correct

  subroutine finish_report(label, t, ok)
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: t
    logical,          intent(in) :: ok
    real(dp) :: tmax
    logical  :: all_ok
    integer  :: e
    call MPI_Reduce(t, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, e)
    call MPI_Allreduce(ok, all_ok, 1, MPI_LOGICAL, MPI_LAND, comm, e)
    if (rank == 0) then
       print '(a,a,a,f9.4,a,f9.2,a,a)', '  ', label, '  total ', tmax, &
            ' s   per-exchange ', tmax/real(nsteps,dp)*1.0e6_dp, ' us', &
            merge('   HALO OK  ', '  HALO FAIL ', all_ok)
    end if
  end subroutine finish_report

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) nx
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) ny
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg); read(arg,*) nsteps
    end if
  end subroutine read_cli

end program gpu_halo

! ===========================================================================
! TODO 5 — write your answers in notes/day5.md
!
! (a) Count the bus crossings per exchange for A and for B. Then, at your
!     halo size, compute the PCIe time for A and compare to the measured
!     difference. Do they agree? If B is not faster on your cluster, GPU
!     awareness may be silently falling back to staging inside the MPI
!     library -- how would you tell?
!
! (b) How did you PROVE your MPI is GPU-aware? "The docs say so" is not a
!     proof. Give a measurement or a tool output that settles it.
!
! (c) Variant C: what fraction of the exchange did you hide? Sweep the tile
!     size -- at what interior:halo ratio does overlap stop helping? Compare
!     the shape of that curve with the one you measured in Exercise 03 on
!     CPU. Are they the same shape?
!
! (d) Multi-GPU per node: with 4 GPUs and 4 ranks per node, two of your
!     neighbours are on the SAME node. Should those go over MPI at all?
!     Look up CUDA IPC / NVLink peer-to-peer and say what an MPI library
!     does with an intra-node device-to-device transfer.
!
! (e) The pack kernel: it gathers a strided row into a contiguous buffer.
!     Estimate its cost relative to the transfer. At what message size does
!     packing dominate, and what would you do about it? (Consider whether
!     the layout could make the halo contiguous in the first place -- and
!     connect this to the nproma layout question from Exercise 08.)
!
! (f) DKRZ context: an ICON GPU run on a machine like JUPITER has thousands
!     of GPUs. Halo exchange is 30+ fields per timestep, each a small
!     message. Given what you measured, which matters more at that scale --
!     bandwidth or message rate? Design the exchange you would recommend,
!     and say what it costs in memory.
! ===========================================================================
