! ===========================================================================
! EXERCISE 20 — Parallel I/O: writing a decomposed field
! ===========================================================================
!
! GOAL
!   Write one distributed 2D field to one file, four ways, and measure:
!     A  file-per-rank         — trivially parallel, operationally awful
!     B  MPI-IO independent    — one file, everyone seeks and writes
!     C  MPI-IO collective     — one file, subarray view, write_all
!     D  gather-to-root        — the accidental serial bottleneck
!   Then tune the collective-buffering hints and see what they buy.
!
! WHY (DKRZ)
!   DKRZ is a data centre. A high-resolution coupled run produces petabytes,
!   and I/O is routinely 20-40% of wall time in a production ESM. It is also
!   the part most often left un-tuned, because it is nobody's specialty.
!
!   Variant D is the one to internalise. Gathering to rank 0 and writing
!   from there is what almost every model does first, because it is easy and
!   it works. It also does not scale at all: one rank's memory, one rank's
!   bandwidth, and every other rank blocked waiting. Recognising that
!   pattern in someone else's model, and knowing the fix, is a concrete
!   thing you can offer.
!
!   Variant A is the trap that looks like a solution: it is genuinely fast,
!   and it is why so many models produce one file per rank. Then someone has
!   to post-process 100000 files, and the filesystem metadata server falls
!   over, and the scientist cannot open the output in anything. Fast for the
!   writer, catastrophic for everyone downstream — a real engineering
!   tradeoff you should be able to argue both sides of.
!
! TASKS
!   TODO 1  write_file_per_rank
!   TODO 2  write_mpiio_independent  — MPI_File_write_at
!   TODO 3  write_mpiio_collective   — subarray datatype + set_view + write_all
!   TODO 4  write_gather_root        — the anti-pattern, measured
!   TODO 5  MPI_Info hints, and verification that the file is correct
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - all four produce a byte-identical file (verified by reading it back
!     and checking every element against its global index)
!   - a bandwidth table in MiB/s for all four
!   - you can explain why collective beats independent even though they
!     write exactly the same bytes
!
! HINTS
!   - MPI_Type_create_subarray is the whole trick: it describes "my tile
!     inside the global array" as a datatype. Then MPI_File_set_view makes
!     the file look like it contains only your tile, and you write your
!     local buffer with no offset arithmetic at all.
!   - Collective I/O wins through two-phase I/O: the library reorganises
!     scattered small writes into a few large contiguous ones by a small set
!     of aggregator ranks. That is why it beats independent writes of
!     identical data.
!   - Hints worth setting:
!         cb_nodes            number of aggregators
!         cb_buffer_size      aggregation buffer
!         striping_factor     Lustre OSTs (cluster only)
!         striping_unit       Lustre stripe size
!   - Always MPI_File_sync or close before timing the end, or you time the
!     page cache and get a beautiful, meaningless number.
! ===========================================================================

program pario
  use mpi_f08
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer :: nx_g = 2048, ny_g = 2048       ! global field
  integer :: nrep = 3

  type(MPI_Comm) :: comm
  integer :: rank, nprocs, ierr
  integer :: px, py, cx, cy                 ! process grid and my position
  integer :: nx_l, ny_l, x0, y0             ! my tile and its global origin
  real(dp), allocatable :: tile(:,:)
  real(dp) :: bw_a, bw_b, bw_c, bw_d
  logical  :: ok_a, ok_b, ok_c, ok_d

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)
  call read_cli()

  call make_process_grid()

  nx_l = nx_g / px
  ny_l = ny_g / py
  x0   = cx * nx_l
  y0   = cy * ny_l

  allocate(tile(nx_l, ny_l))
  call fill_tile()

  if (rank == 0) then
     print '(a)', '=== Exercise 20: parallel I/O ==='
     print '(a,i0,a,i0,a,f10.2,a)', 'global field ', nx_g, ' x ', ny_g, &
          ' = ', real(nx_g,dp)*real(ny_g,dp)*8.0_dp/1048576.0_dp, ' MiB'
     print '(a,i0,a,i0,a,i0,a,i0)', 'process grid ', px, ' x ', py, &
          '   tile ', nx_l, ' x ', ny_l
     print '(a)', ''
     print '(a)', '  method                        time (s)    MiB/s   verified'
     print '(a)', '  ------------------------------------------------------------'
  end if

  call timed('A: file per rank            ', 1, bw_a, ok_a)
  call timed('B: MPI-IO independent       ', 2, bw_b, ok_b)
  call timed('C: MPI-IO collective        ', 3, bw_c, ok_c)
  call timed('D: gather to root           ', 4, bw_d, ok_d)

  if (rank == 0) then
     print '(a)', ''
     print '(a)', '  Expect C > B (two-phase aggregation) and D to be'
     print '(a)', '  catastrophically slow at high rank counts.'
     print '(a)', '  Then try: make hints   and   make ranks'
  end if

  deallocate(tile)
  call MPI_Finalize(ierr)

contains

  !> Factor nprocs into the squarest px x py that divides the global grid.
  subroutine make_process_grid()
    integer :: p
    px = int(sqrt(real(nprocs)))
    do while (px > 1)
       if (mod(nprocs, px) == 0 .and. mod(nx_g, px) == 0 .and. &
           mod(ny_g, nprocs/px) == 0) exit
       px = px - 1
    end do
    if (px < 1) px = 1
    py = nprocs / px
    cx = mod(rank, px)
    cy = rank / px
    p = px * py
    if (p /= nprocs .and. rank == 0) &
         print '(a)', '  WARNING: process grid does not use all ranks'
  end subroutine make_process_grid

  !> Each element holds its GLOBAL linear index, so verification is exact:
  !> read the file back and element k must equal k.
  subroutine fill_tile()
    integer :: i, j, gi, gj
    do j = 1, ny_l
       do i = 1, nx_l
          gi = x0 + i
          gj = y0 + j
          tile(i,j) = real((gj - 1) * nx_g + gi, dp)
       end do
    end do
  end subroutine fill_tile

  subroutine timed(label, which, bw, ok)
    character(len=*), intent(in)  :: label
    integer,          intent(in)  :: which
    real(dp),         intent(out) :: bw
    logical,          intent(out) :: ok
    real(dp) :: t, tmax, mib
    integer  :: r, e

    call MPI_Barrier(comm, e)
    t = MPI_Wtime()
    do r = 1, nrep
       select case (which)
       case (1); call write_file_per_rank()
       case (2); call write_mpiio_independent()
       case (3); call write_mpiio_collective()
       case (4); call write_gather_root()
       end select
    end do
    call MPI_Barrier(comm, e)
    t = (MPI_Wtime() - t) / real(nrep, dp)

    call MPI_Reduce(t, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, e)
    ok = verify_output(which)

    mib = real(nx_g,dp) * real(ny_g,dp) * 8.0_dp / 1048576.0_dp
    if (rank == 0) then
       bw = mib / max(tmax, 1.0e-9_dp)
       ! 1e-4 s to write tens of MiB would be hundreds of GB/s -- no real
       ! filesystem does that, so a time this small means the routine is
       ! still a stub. Printing the "bandwidth" would be pure noise.
       if (tmax < 1.0e-4_dp) then
          print '(a,a,a)', '  ', label, '     -- not implemented --'
          bw = 0.0_dp
       else
          print '(a,a,f10.4,f11.1,a)', '  ', label, tmax, bw, &
               merge('       yes', '        NO', ok)
       end if
    end if
  end subroutine timed

  ! ------------------------------------------------------------------------
  ! TODO 1: one file per rank.
  !
  ! Plain Fortran unformatted stream I/O to 'out_rank_NNNN.bin'. No MPI
  ! involved. This will probably be the FASTEST option here, which is
  ! exactly why the pattern is so common and so damaging — think about who
  ! pays the cost, and when.
  ! ------------------------------------------------------------------------
  subroutine write_file_per_rank()
    character(len=64) :: fname
    integer :: u
    write(fname, '(a,i4.4,a)') 'out_rank_', rank, '.bin'
    ! TODO 1: open stream, write tile, close
    u = 0
    if (len_trim(fname) < 0) continue
  end subroutine write_file_per_rank

  ! ------------------------------------------------------------------------
  ! TODO 2: MPI-IO, independent writes.
  !
  !   MPI_File_open(comm, 'out_indep.bin', MPI_MODE_CREATE+MPI_MODE_WRONLY,
  !                 MPI_INFO_NULL, fh)
  !   then for each of my rows j, compute the byte offset of that row's
  !   segment in the global array and MPI_File_write_at it.
  !
  ! One call per ROW, because my tile is not contiguous in the global file:
  ! my row j occupies bytes [(y0+j-1)*nx_g + x0 ... + nx_l] and the next row
  ! is nx_g away, not nx_l. That fragmentation is the whole problem, and
  ! TODO 3 is how you stop doing it by hand.
  ! ------------------------------------------------------------------------
  subroutine write_mpiio_independent()
    type(MPI_File) :: fh
    integer :: e
    ! TODO 2
    e = 0
    if (.false.) call MPI_File_close(fh, e)
  end subroutine write_mpiio_independent

  ! ------------------------------------------------------------------------
  ! TODO 3: MPI-IO, collective, with a subarray view. The right answer.
  !
  !   integer :: sizes(2), subsizes(2), starts(2)
  !   sizes    = [nx_g, ny_g]
  !   subsizes = [nx_l, ny_l]
  !   starts   = [x0,   y0]          ! ZERO-based, even in Fortran
  !
  !   call MPI_Type_create_subarray(2, sizes, subsizes, starts,          &
  !          MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, e)
  !   call MPI_Type_commit(filetype, e)
  !   call MPI_File_open(comm, 'out_coll.bin', ..., fh, e)
  !   call MPI_File_set_view(fh, 0_MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION,  &
  !                          filetype, 'native', info, e)
  !   call MPI_File_write_all(fh, tile, nx_l*ny_l, MPI_DOUBLE_PRECISION,   &
  !                           MPI_STATUS_IGNORE, e)
  !   call MPI_File_close(fh, e); call MPI_Type_free(filetype, e)
  !
  ! ONE call, no offset arithmetic. The `_all` suffix is what lets the
  ! library aggregate across ranks — MPI_File_write (no _all) with the same
  ! view would be correct but would not aggregate.
  !
  ! TODO 5: build an MPI_Info with cb_nodes / cb_buffer_size and pass it to
  ! set_view instead of MPI_INFO_NULL. Sweep the values.
  ! ------------------------------------------------------------------------
  subroutine write_mpiio_collective()
    type(MPI_File)     :: fh
    type(MPI_Datatype) :: filetype
    type(MPI_Info)     :: info
    integer :: e
    ! TODO 3
    e = 0
    if (.false.) then
       call MPI_File_close(fh, e)
       call MPI_Type_free(filetype, e)
       call MPI_Info_free(info, e)
    end if
  end subroutine write_mpiio_collective

  ! ------------------------------------------------------------------------
  ! TODO 4: the anti-pattern, measured honestly.
  !
  ! MPI_Gather every tile to rank 0, reassemble the global array there, and
  ! write it with plain Fortran I/O.
  !
  ! Note what this costs before you run it: rank 0 needs the WHOLE global
  ! array in memory (at 2048^2 that is 32 MiB, at a realistic ESM resolution
  ! it is hundreds of GB and simply impossible), and every other rank sits
  ! idle. Measure it anyway — you need the number to argue against it.
  ! ------------------------------------------------------------------------
  subroutine write_gather_root()
    real(dp), allocatable :: global(:,:)
    integer :: e
    ! TODO 4
    e = 0
    if (allocated(global)) deallocate(global)
  end subroutine write_gather_root

  ! ------------------------------------------------------------------------
  ! TODO 5b: verify. Read the file back on rank 0 and check that element k
  ! holds the value k. A wrong subarray `starts` produces a file of exactly
  ! the right SIZE with the tiles in the wrong places — size checks will not
  ! catch it, and neither will a plot unless you look carefully.
  ! ------------------------------------------------------------------------
  logical function verify_output(which) result(ok)
    integer, intent(in) :: which
    ok = .false.
    ! TODO 5b: read back the file written by `which` and check every element.
    ! File-per-rank needs a different check (many files) -- decide what
    ! "verified" means for it and be explicit about the difference.
    if (which < 0) continue
  end function verify_output

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) nx_g
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) ny_g
    end if
  end subroutine read_cli

end program pario

! ===========================================================================
! TODO 6 — write your answers in notes/day5.md
!
! (a) Bandwidth table for all four at 4, 8, 16 ranks. Which scales and which
!     does not? Plot D's time against rank count and state its asymptotic
!     behaviour.
!
! (b) B and C write exactly the same bytes to exactly the same file, yet C
!     is faster. Explain two-phase I/O: who aggregates, what buffer is used,
!     and why fewer larger writes beat many smaller ones on a parallel
!     filesystem.
!
! (c) `make hints`: sweep cb_nodes. Is more always better? Relate the
!     optimum to the number of OSTs (Lustre) or the number of nodes.
!
! (d) File-per-rank is fastest here. Write the argument AGAINST it that you
!     would give a scientist who likes it: think about 100000 files, the
!     metadata server, post-processing, archival, and what happens when
!     someone reruns on a different rank count.
!
! (e) Real ESM output is netCDF/HDF5, not raw binary, because it is
!     self-describing. What does the parallel HDF5 layer add on top of what
!     you built, and what does it cost? Try it: the toolchain here has
!     hdf5-mpi (`h5pcc -show`). Compare against variant C.
!
! (f) DKRZ context: production ESMs often use dedicated I/O SERVERS -- a
!     subset of ranks that do nothing but receive fields and write them,
!     so the compute ranks never block on I/O. Sketch that design. What is
!     the buffering requirement, and what happens when the writers fall
!     behind? (Search terms: XIOS, CDI-PIO. CDI-PIO is the DKRZ one.)
! ===========================================================================
