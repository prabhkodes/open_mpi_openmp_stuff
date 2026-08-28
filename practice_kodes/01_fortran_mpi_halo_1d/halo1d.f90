! ===========================================================================
! EXERCISE 01 — 1D halo exchange in modern Fortran with mpi_f08
! ===========================================================================
!
! GOAL
!   Decompose a 1D field across MPI ranks and keep one layer of ghost cells
!   up to date, three ways: blocking Sendrecv, non-blocking Isend/Irecv, and
!   persistent requests. Then measure which is fastest and at what message
!   size the ordering changes.
!
! WHY (DKRZ)
!   This is the single most executed communication pattern in any Earth
!   System Model. ICON does it every timestep on every prognostic variable.
!   The job asks for "parallel programming with MPI/OpenMP" and "Fortran
!   under UNIX/LINUX" — this exercise is the intersection of both.
!   Note the `use mpi_f08` interface: handles are derived TYPES
!   (type(MPI_Comm), type(MPI_Request)) rather than bare integers, so the
!   compiler catches argument-order mistakes that the old `mpi` module
!   silently accepted. Modern ICON code uses mpi_f08.
!
! TASKS
!   TODO 1  decompose_1d      — split n_global over nprocs with the remainder
!                               spread over the first mod(n_global,nprocs) ranks
!   TODO 2  exchange_blocking — one MPI_Sendrecv per direction
!   TODO 3  exchange_nonblock — 2 x MPI_Irecv + 2 x MPI_Isend + MPI_Waitall
!   TODO 4  exchange_persist  — MPI_Send_init/Recv_init once, MPI_Startall in
!                               the loop (this is what ICON actually does)
!   TODO 5  answer the question at the bottom of this file
!
! ACCEPTANCE
!   - `make run NP=4` prints "HALO CHECK: PASS" for all three variants
!   - `make sweep` produces a table of time vs halo width for all three
!   - You can state the crossover point where non-blocking beats blocking,
!     and explain why persistent mode wins at small messages
!
! HINTS
!   - MPI_PROC_NULL is a valid rank: send/recv to it is a no-op that returns
!     immediately. Use it for the physical boundaries instead of if-branching.
!   - Post ALL receives before ANY sends. If you send first you rely on
!     eager-protocol buffering, which silently deadlocks above the eager
!     threshold (typically ~64 KB with OpenMPI).
!   - Persistent requests amortise the request-object setup, which is a real
!     cost when you exchange 30 small fields per timestep.
! ===========================================================================


program halo1d
  use mpi_f08
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! -- configuration (override on the command line: ./halo1d 1000000 4 200)
  integer :: n_global = 1000000   ! total interior points
  integer :: halo     = 1         ! ghost layers on each side
  integer :: n_iter   = 200       ! exchanges to time

  type(MPI_Comm) :: comm
  integer :: rank, nprocs, ierr
  integer :: n_local, i_start          ! this rank's slice
  integer :: left, right               ! neighbour ranks (or MPI_PROC_NULL)
  real(dp), allocatable :: u(:)        ! u(1-halo : n_local+halo)
  real(dp) :: t_block, t_nb, t_pers

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)

  call read_cli()

  ! ---- decomposition -----------------------------------------------------
  call decompose_1d(n_global, nprocs, rank, n_local, i_start)

  ! Periodic domain: rank 0's left neighbour is the last rank.
  ! Switch these two lines to MPI_PROC_NULL at the ends for a non-periodic
  ! domain and confirm the halo check still passes.
  left  = modulo(rank - 1, nprocs)
  right = modulo(rank + 1, nprocs)

  allocate(u(1-halo : n_local+halo))

  if (rank == 0) then
     print '(a)',        '=== Exercise 01: 1D halo exchange ==='
     print '(a,i0,a,i0,a,i0)', 'n_global=', n_global, '  nprocs=', nprocs, &
                               '  halo=', halo
     print '(a,i0)',     'bytes per exchange per side: ', halo * 8
     print '(a)',        ''
  end if

  ! ---- variant 1: blocking Sendrecv --------------------------------------
  call fill_with_global_index()
  call MPI_Barrier(comm, ierr)
  t_block = MPI_Wtime()
  block
    integer :: it
    do it = 1, n_iter
       call exchange_blocking()
    end do
  end block
  t_block = MPI_Wtime() - t_block
  call report('blocking Sendrecv', t_block)

  ! ---- variant 2: non-blocking Isend/Irecv -------------------------------
  call fill_with_global_index()
  call MPI_Barrier(comm, ierr)
  t_nb = MPI_Wtime()
  block
    integer :: it
    do it = 1, n_iter
       call exchange_nonblock()
    end do
  end block
  t_nb = MPI_Wtime() - t_nb
  call report('non-blocking Isend/Irecv', t_nb)

  ! ---- variant 3: persistent requests ------------------------------------
  call fill_with_global_index()
  call MPI_Barrier(comm, ierr)
  t_pers = MPI_Wtime()
  call exchange_persistent(n_iter)
  t_pers = MPI_Wtime() - t_pers
  call report('persistent Startall', t_pers)

  if (rank == 0) then
     print '(a)', ''
     print '(a)', 'Record these three timings in notes/day1.md, then re-run'
     print '(a)', 'with:  make sweep     (halo = 1,4,16,64,256,1024,4096)'
  end if

  deallocate(u)
  call MPI_Finalize(ierr)

contains

  ! ------------------------------------------------------------------------
  ! TODO 1: block decomposition with remainder distribution.
  !
  ! With n=10 and p=4 the slices must be 3,3,2,2 — NOT 2,2,2,4. Getting this
  ! wrong is the classic source of load imbalance in ESMs, because the last
  ! rank ends up doing double work and every collective waits for it.
  !
  ! Set:
  !   nloc   = points owned by `r`
  !   istart = 1-based global index of this rank's first point
  ! and make sure sum(nloc) over all ranks == n.
  ! ------------------------------------------------------------------------
  subroutine decompose_1d(n, p, r, nloc, istart)
    integer, intent(in)  :: n, p, r
    integer, intent(out) :: nloc, istart

    nloc   = n / p        ! TODO 1: add the remainder correction
    istart = r * nloc + 1 ! TODO 1: fix istart to match

    ! Leave this check in — it will fail loudly until TODO 1 is correct.
    block
      integer :: total, e
      call MPI_Allreduce(nloc, total, 1, MPI_INTEGER, MPI_SUM, comm, e)
      if (total /= n .and. rank == 0) then
         print '(a,i0,a,i0)', '  [decompose_1d] WARNING: slices sum to ', &
              total, ' but n_global = ', n
      end if
    end block
  end subroutine decompose_1d

  ! ------------------------------------------------------------------------
  ! TODO 2: blocking exchange with MPI_Sendrecv.
  !
  ! Two calls total. Send your rightmost `halo` interior cells to `right`
  ! while receiving into your left ghost cells from `left`, then mirror it.
  !
  !   interior right edge : u(n_local-halo+1 : n_local)
  !   left ghost cells    : u(1-halo : 0)
  !   interior left edge  : u(1 : halo)
  !   right ghost cells   : u(n_local+1 : n_local+halo)
  !
  ! MPI_Sendrecv(sendbuf, sendcount, sendtype, dest,   sendtag, &
  !              recvbuf, recvcount, recvtype, source, recvtag, &
  !              comm, status, ierror)
  ! ------------------------------------------------------------------------
  subroutine exchange_blocking()
    integer :: e
    type(MPI_Status) :: st

    ! TODO 2a: send right edge -> right neighbour, recv into left ghosts
    ! TODO 2b: send left  edge -> left  neighbour, recv into right ghosts

    e = 0; st = MPI_STATUS_IGNORE   ! placeholder so the skeleton compiles
  end subroutine exchange_blocking

  ! ------------------------------------------------------------------------
  ! TODO 3: non-blocking exchange.
  !
  ! Order matters: post BOTH MPI_Irecv calls first, then both MPI_Isend,
  ! then a single MPI_Waitall over all four requests. Posting receives first
  ! lets the MPI progress engine land incoming data directly in your buffer
  ! (rendezvous protocol) instead of staging it through an internal one.
  !
  ! Once this works, look at exercise 03 — the whole point of non-blocking
  ! is that you can compute the interior between the Isend and the Waitall.
  ! ------------------------------------------------------------------------
  subroutine exchange_nonblock()
    type(MPI_Request) :: req(4)
    integer :: e

    ! TODO 3a: req(1) = Irecv into left  ghosts from left
    ! TODO 3b: req(2) = Irecv into right ghosts from right
    ! TODO 3c: req(3) = Isend right edge to right
    ! TODO 3d: req(4) = Isend left  edge to left
    ! TODO 3e: MPI_Waitall(4, req, MPI_STATUSES_IGNORE, e)

    req = MPI_REQUEST_NULL; e = 0   ! placeholder
  end subroutine exchange_nonblock

  ! ------------------------------------------------------------------------
  ! TODO 4: persistent requests.
  !
  ! MPI_Send_init / MPI_Recv_init build the request objects ONCE, outside the
  ! timestep loop. Inside the loop you only call MPI_Startall + MPI_Waitall.
  ! This removes per-call argument marshalling and lets the MPI library
  ! pre-register the buffers with the NIC.
  !
  ! ICON uses exactly this pattern, which is why its halo buffers are
  ! allocated once at setup and reused for the whole run.
  !
  ! Remember MPI_Request_free on all four at the end.
  ! ------------------------------------------------------------------------
  subroutine exchange_persistent(iters)
    integer, intent(in) :: iters
    type(MPI_Request) :: req(4)
    integer :: it, e

    ! TODO 4a: MPI_Recv_init x2, MPI_Send_init x2  (same buffers as TODO 3)

    req = MPI_REQUEST_NULL
    do it = 1, iters
       ! TODO 4b: MPI_Startall(4, req, e)
       ! TODO 4c: MPI_Waitall(4, req, MPI_STATUSES_IGNORE, e)
    end do

    ! TODO 4d: MPI_Request_free on each request
    e = 0   ! placeholder
  end subroutine exchange_persistent

  ! ------------------------------------------------------------------------
  ! Test harness below this line — you should not need to modify it.
  ! ------------------------------------------------------------------------

  !> Fill interior with the global index so the correct halo contents are
  !> known analytically. This makes the check exact, not approximate.
  subroutine fill_with_global_index()
    integer :: i
    u = -huge(1.0_dp)                      ! poison the ghosts
    do i = 1, n_local
       u(i) = real(i_start + i - 1, dp)
    end do
  end subroutine fill_with_global_index

  !> After a correct exchange, ghost cell j must hold the global index of the
  !> point it shadows, wrapped periodically into [1, n_global].
  logical function halo_is_correct() result(ok)
    integer :: k, gidx
    real(dp) :: want
    ok = .true.
    do k = 1, halo
       gidx = modulo(i_start - 1 - k, n_global) + 1     ! left ghost u(1-k)
       want = real(gidx, dp)
       if (abs(u(1-k) - want) > 0.0_dp) ok = .false.

       gidx = modulo(i_start - 1 + n_local + k - 1, n_global) + 1
       want = real(gidx, dp)
       if (abs(u(n_local+k) - want) > 0.0_dp) ok = .false.
    end do
  end function halo_is_correct

  subroutine report(label, t)
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: t
    logical :: ok, all_ok
    real(dp) :: tmax
    integer :: e

    ok = halo_is_correct()
    call MPI_Allreduce(ok, all_ok, 1, MPI_LOGICAL, MPI_LAND, comm, e)
    call MPI_Reduce(t, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, e)

    if (rank == 0) then
       print '(a,a24,a,f9.4,a,f9.3,a)', '  ', label, &
            '  total ', tmax, ' s   per-exchange ', &
            tmax / real(n_iter, dp) * 1.0e6_dp, ' us'
       if (all_ok) then
          print '(a)', '     HALO CHECK: PASS'
       else
          print '(a)', '     HALO CHECK: FAIL  <-- ghost cells are wrong'
       end if
    end if
  end subroutine report

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) n_global
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) halo
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg); read(arg,*) n_iter
    end if
  end subroutine read_cli

end program halo1d

! ===========================================================================
! TODO 5 — write your answers in notes/day1.md
!
! (a) At halo=1 (8 bytes) which variant wins, and at halo=4096 (32 KB)?
!     Where is the crossover, and what does that tell you about the fixed
!     cost of an MPI call versus the bandwidth cost?
!
! (b) OpenMPI switches from the eager to the rendezvous protocol at
!     btl_*_eager_limit. Find it with:
!         ompi_info --param btl all --level 9 | grep eager
!     Does your measured crossover line up with it?
!
! (c) Why must you post receives before sends in TODO 3? Construct the
!     deadlock: make exchange_blocking use MPI_Send + MPI_Recv in that order
!     and run with halo=100000. Explain what you see.
!
! (d) ICON exchanges ~30 fields per timestep. If each is a separate
!     Isend/Irecv pair, what optimisation would you propose, and what does
!     it cost you in memory?
! ===========================================================================
