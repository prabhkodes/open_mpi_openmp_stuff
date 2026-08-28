! ===========================================================================
! EXERCISE 06 — Halo exchange on an UNSTRUCTURED mesh
! ===========================================================================
!
! GOAL
!   Exchange ghost data when there is no regular neighbour pattern: build
!   the send/receive index lists from scratch, pack into contiguous buffers,
!   exchange, unpack. Then do the same thing with a distributed graph
!   communicator and MPI_Neighbor_alltoallv, and compare.
!
! WHY (DKRZ)
!   In Exercise 01 the neighbour of rank r was r+1. On ICON's icosahedral
!   grid a rank's neighbours are whoever happens to own the cells adjacent
!   to its partition boundary — could be 4 ranks, could be 11, and it is
!   different for every rank. There is no formula; you must DISCOVER the
!   communication pattern at setup time and then replay it every timestep.
!
!   The hard part is not the send. It is that YOU know what you need to
!   RECEIVE, but the sender does not know what to send you. Inverting that
!   — turning "here is what I need" into "here is what I must send" — is
!   the core of every unstructured-mesh code, and it is TODO 2.
!
!   This is the exercise that most directly separates "I have used MPI"
!   from "I have worked on an unstructured HPC code". Do not skip it.
!
! MESH
!   To keep the exercise self-contained, the mesh is a 2D nx-by-ny grid of
!   cells, but stored the way an unstructured code stores it: an explicit
!   cell_neigh(4, nglobal) table of GLOBAL indices, with no (i,j) arithmetic
!   anywhere in the solver. That gives exact verification while keeping the
!   irregular communication structure real — with a linear partition the
!   north/south neighbours land in a different rank's block, so the pattern
!   is genuinely discovered, not assumed.
!
! TASKS
!   TODO 1  build_halo         — which remote cells do I need?
!   TODO 2  invert_to_send     — turn my recv lists into everyone's send lists
!   TODO 3  exchange_p2p       — pack, Isend/Irecv, Waitall, unpack
!   TODO 4  build_graph_comm   — MPI_Dist_graph_create_adjacent
!   TODO 5  exchange_neighbor  — MPI_Neighbor_alltoallv
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - "HALO CHECK: PASS" for both exchange implementations
!   - a printed table of halo cells per rank, neighbour ranks per rank,
!     bytes exchanged, and time per exchange for both methods
!   - you can explain when Neighbor_alltoallv beats hand-rolled p2p
!
! HINTS
!   - TODO 2 has a standard solution: MPI_Alltoall of the per-rank COUNTS
!     (so everyone learns how many items each peer wants), then
!     MPI_Alltoallv of the actual global indices. Two collectives at setup,
!     zero collectives per timestep.
!   - Sort each rank's request list by global index. Then the sender packs
!     in a deterministic order and the receiver unpacks in the same order,
!     with no tag matching needed.
!   - The global->local map here is a full-size array for clarity. Note in
!     your answer why a real code at 1e8 cells cannot afford that, and what
!     it uses instead.
! ===========================================================================

module uhalo_mod
  use mpi_f08
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  type :: mesh_t
     integer :: nx = 0, ny = 0, nglobal = 0
     integer, allocatable :: cell_neigh(:,:)   ! (4, nglobal), global indices

     ! --- this rank's partition -------------------------------------------
     integer :: nowned = 0, nhalo = 0
     integer, allocatable :: owned_gidx(:)     ! (nowned) local -> global
     integer, allocatable :: halo_gidx(:)      ! (nhalo)  halo slot -> global
     integer, allocatable :: g2l(:)            ! (nglobal) global -> local, 0=absent
     integer, allocatable :: owner(:)          ! (nglobal) owning rank

     ! --- communication schedule (built once, replayed every timestep) ------
     integer :: nneigh = 0
     integer, allocatable :: neigh_rank(:)     ! (nneigh)
     integer, allocatable :: recv_count(:)     ! (nneigh)
     integer, allocatable :: recv_displ(:)     ! (nneigh)
     integer, allocatable :: recv_slot(:)      ! flat: halo slot per recv item
     integer, allocatable :: send_count(:)     ! (nneigh)
     integer, allocatable :: send_displ(:)     ! (nneigh)
     integer, allocatable :: send_lidx(:)      ! flat: owned local idx per item

     real(dp), allocatable :: sendbuf(:), recvbuf(:)
     type(MPI_Comm) :: graph_comm
     logical :: graph_ready = .false.
  end type mesh_t

contains

  !> Build the global mesh topology. Every rank builds the whole table --
  !> fine at exercise scale, and it keeps the exercise focused on the
  !> COMMUNICATION rather than on distributed mesh reading. Periodic in both
  !> directions, so every cell has exactly 4 neighbours and there are no
  !> boundary special cases to muddy the verification.
  subroutine build_mesh(m, nx, ny)
    type(mesh_t), intent(out) :: m
    integer, intent(in) :: nx, ny
    integer :: i, j, c

    m%nx = nx; m%ny = ny; m%nglobal = nx * ny
    allocate(m%cell_neigh(4, m%nglobal))

    do j = 1, ny
       do i = 1, nx
          c = (j-1)*nx + i
          m%cell_neigh(1, c) = (j-1)*nx + modulo(i-2, nx) + 1     ! west
          m%cell_neigh(2, c) = (j-1)*nx + modulo(i,   nx) + 1     ! east
          m%cell_neigh(3, c) = modulo(j-2, ny)*nx + i             ! south
          m%cell_neigh(4, c) = modulo(j,   ny)*nx + i             ! north
       end do
    end do
  end subroutine build_mesh

  !> Linear partition: rank r owns a contiguous block of global indices.
  !> Exercise 07 replaces this with something much better.
  subroutine partition_linear(m, rank, nprocs)
    type(mesh_t), intent(inout) :: m
    integer, intent(in) :: rank, nprocs
    integer :: base, rem, r, lo, hi, i

    allocate(m%owner(m%nglobal))
    base = m%nglobal / nprocs
    rem  = mod(m%nglobal, nprocs)

    lo = 1
    do r = 0, nprocs - 1
       hi = lo + base - 1
       if (r < rem) hi = hi + 1
       m%owner(lo:hi) = r
       if (r == rank) then
          m%nowned = hi - lo + 1
          allocate(m%owned_gidx(m%nowned))
          do i = 1, m%nowned
             m%owned_gidx(i) = lo + i - 1
          end do
       end if
       lo = hi + 1
    end do

    allocate(m%g2l(m%nglobal))
    m%g2l = 0
    do i = 1, m%nowned
       m%g2l(m%owned_gidx(i)) = i
    end do
  end subroutine partition_linear

  ! ------------------------------------------------------------------------
  ! TODO 1: discover the halo.
  !
  ! Walk every owned cell's neighbour list. Any neighbour whose owner is not
  ! me is a halo cell. Collect the UNIQUE set of such global indices, sort
  ! it, and assign each one a halo slot.
  !
  ! Fill in:
  !   m%nhalo
  !   m%halo_gidx(1:nhalo)   -- sorted ascending
  !   m%g2l(gidx) = m%nowned + slot     for each halo cell
  !
  ! Why sorted: it makes the recv order deterministic and lets the sender
  ! reproduce it without extra metadata. Every unstructured code does this.
  ! ------------------------------------------------------------------------
  subroutine build_halo(m, rank)
    type(mesh_t), intent(inout) :: m
    integer, intent(in) :: rank

    ! TODO 1: replace this stub.
    !   - loop c = 1..nowned, k = 1..4 over m%cell_neigh(k, m%owned_gidx(c))
    !   - if m%owner(gn) /= rank, remember gn
    !   - unique + sort (a mark array over nglobal is the easy way here)
    !   - allocate m%halo_gidx, fill it, and extend m%g2l
    m%nhalo = 0
    allocate(m%halo_gidx(max(m%nhalo,1)))
    m%halo_gidx = 0
  end subroutine build_halo

  ! ------------------------------------------------------------------------
  ! TODO 2: invert the recv schedule into a send schedule.
  !
  ! THIS IS THE HEART OF THE EXERCISE.
  !
  ! You know: "I need cells 400,401,999 from rank 2 and 1500 from rank 5."
  ! Rank 2 knows nothing about that. It must learn what to pack.
  !
  ! Standard two-step, both at setup time only:
  !
  !   step A  every rank counts how many items it wants from every other
  !           rank -> want(0:nprocs-1). One MPI_Alltoall turns that into
  !           give(0:nprocs-1): how many items each peer wants FROM ME.
  !
  !   step B  MPI_Alltoallv of the actual global indices, using want/give as
  !           the counts. Now I hold the list of global indices each peer
  !           needs from me; map them through g2l to get local indices ->
  !           m%send_lidx.
  !
  ! Then compress to the peers with nonzero traffic:
  !   m%nneigh, m%neigh_rank, m%send_count/displ, m%recv_count/displ
  !
  ! Also fill m%recv_slot: for recv item i, which halo slot it lands in.
  ! Because halo_gidx is sorted and grouped by owner, this is just a
  ! running index — but write it explicitly, you will need it for TODO 5.
  ! ------------------------------------------------------------------------
  subroutine invert_to_send(m, comm, rank, nprocs)
    type(mesh_t),   intent(inout) :: m
    type(MPI_Comm), intent(in)    :: comm
    integer,        intent(in)    :: rank, nprocs

    ! TODO 2: replace this stub with the Alltoall + Alltoallv inversion.
    m%nneigh = 0
    allocate(m%neigh_rank(1), m%send_count(1), m%send_displ(1))
    allocate(m%recv_count(1), m%recv_displ(1))
    allocate(m%send_lidx(1), m%recv_slot(1))
    m%neigh_rank = 0; m%send_count = 0; m%send_displ = 0
    m%recv_count = 0; m%recv_displ = 0; m%send_lidx = 1; m%recv_slot = 1
    allocate(m%sendbuf(1), m%recvbuf(1))
    m%sendbuf = 0.0_dp; m%recvbuf = 0.0_dp
    if (rank < 0 .or. nprocs < 0) continue
    if (comm == MPI_COMM_NULL) continue
  end subroutine invert_to_send

  ! ------------------------------------------------------------------------
  ! TODO 3: the point-to-point exchange, replayed every timestep.
  !
  !   pack:     sendbuf(i) = f(send_lidx(i))          for all i
  !   exchange: one Irecv + one Isend per neighbour, then Waitall
  !   unpack:   f(nowned + recv_slot(i)) = recvbuf(i) for all i
  !
  ! Keep the pack/unpack loops separate from the MPI calls — that is what
  ! lets you OpenMP-parallelise the packing (it is a gather over 1e5+ items
  ! in a real model) and what lets you overlap, as in Exercise 03.
  ! ------------------------------------------------------------------------
  subroutine exchange_p2p(m, f, comm)
    type(mesh_t),   intent(inout) :: m
    real(dp),       intent(inout) :: f(:)
    type(MPI_Comm), intent(in)    :: comm
    type(MPI_Request), allocatable :: req(:)
    integer :: n, e

    n = max(m%nneigh, 1)
    allocate(req(2*n))
    req = MPI_REQUEST_NULL

    ! TODO 3a: pack   -- sendbuf from f via send_lidx
    ! TODO 3b: post Irecv for every neighbour into recvbuf at recv_displ
    ! TODO 3c: post Isend for every neighbour from sendbuf at send_displ
    ! TODO 3d: MPI_Waitall(2*m%nneigh, req, MPI_STATUSES_IGNORE, e)
    ! TODO 3e: unpack -- f(nowned + recv_slot) from recvbuf

    e = 0
    deallocate(req)
    if (size(f) < 0) continue
    if (comm == MPI_COMM_NULL) continue
  end subroutine exchange_p2p

  ! ------------------------------------------------------------------------
  ! TODO 4: build a distributed graph communicator.
  !
  !   call MPI_Dist_graph_create_adjacent(comm,                        &
  !          indegree,  sources,      sourceweights,                   &
  !          outdegree, destinations, destweights,                     &
  !          MPI_INFO_NULL, .false., m%graph_comm, ierr)
  !
  ! Here sources == destinations == m%neigh_rank (the pattern is symmetric:
  ! if I need data from you, you need data from me — true for a shared-face
  ! halo, NOT true in general, so do not assume it in other codes).
  !
  ! Passing `reorder = .true.` lets the MPI library renumber ranks to match
  ! the machine topology. On a real cluster that is free performance; try
  ! both and see whether your MPI actually implements it.
  ! ------------------------------------------------------------------------
  subroutine build_graph_comm(m, comm)
    type(mesh_t),   intent(inout) :: m
    type(MPI_Comm), intent(in)    :: comm
    ! TODO 4: create m%graph_comm, then set m%graph_ready = .true.
    m%graph_comm = comm
    m%graph_ready = .false.
  end subroutine build_graph_comm

  ! ------------------------------------------------------------------------
  ! TODO 5: the same exchange as one collective call.
  !
  !   call MPI_Neighbor_alltoallv(m%sendbuf, m%send_count, m%send_displ,   &
  !                               MPI_DOUBLE_PRECISION,                    &
  !                               m%recvbuf, m%recv_count, m%recv_displ,   &
  !                               MPI_DOUBLE_PRECISION, m%graph_comm, ierr)
  !
  ! Same pack and unpack as TODO 3 — only the middle changes. The MPI
  ! library now knows the whole pattern up front and may schedule it far
  ! better than your loop of Isends (message aggregation, topology-aware
  ! ordering, hardware collectives on some fabrics).
  !
  ! Whether it actually IS faster is an empirical question. Measure it.
  ! ------------------------------------------------------------------------
  subroutine exchange_neighbor(m, f)
    type(mesh_t), intent(inout) :: m
    real(dp),     intent(inout) :: f(:)
    integer :: e
    ! TODO 5a: pack
    ! TODO 5b: MPI_Neighbor_alltoallv
    ! TODO 5c: unpack
    e = 0
    if (size(f) < 0) continue
  end subroutine exchange_neighbor

end module uhalo_mod


program uhalo
  use uhalo_mod
  use mpi_f08
  implicit none

  type(mesh_t) :: m
  type(MPI_Comm) :: comm
  integer :: rank, nprocs, ierr
  integer :: nx = 400, ny = 400, n_iter = 200
  real(dp), allocatable :: f(:)
  real(dp) :: t_p2p, t_nbr
  integer :: max_halo, max_neigh, tot_bytes

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)
  call read_cli()

  ! ---- setup: done ONCE, cost amortised over the whole run ---------------
  call build_mesh(m, nx, ny)
  call partition_linear(m, rank, nprocs)
  call build_halo(m, rank)
  call invert_to_send(m, comm, rank, nprocs)
  call build_graph_comm(m, comm)

  allocate(f(m%nowned + m%nhalo))

  call print_header()

  ! ---- variant 1: point-to-point -----------------------------------------
  call fill_with_global_index()
  call MPI_Barrier(comm, ierr)
  t_p2p = MPI_Wtime()
  block
    integer :: it
    do it = 1, n_iter
       call exchange_p2p(m, f, comm)
    end do
  end block
  t_p2p = MPI_Wtime() - t_p2p
  call report('point-to-point Isend/Irecv', t_p2p)

  ! ---- variant 2: neighbourhood collective -------------------------------
  call fill_with_global_index()
  call MPI_Barrier(comm, ierr)
  t_nbr = MPI_Wtime()
  block
    integer :: it
    do it = 1, n_iter
       call exchange_neighbor(m, f)
    end do
  end block
  t_nbr = MPI_Wtime() - t_nbr
  call report('MPI_Neighbor_alltoallv    ', t_nbr)

  if (rank == 0) then
     print '(a)', ''
     print '(a)', '  Record both timings and the halo/neighbour stats in'
     print '(a)', '  notes/day2.md, then try:  make strong   and   make shape'
  end if

  deallocate(f)
  call MPI_Finalize(ierr)

contains

  !> Poison everything, then write the exact global index into owned cells.
  !> After a correct exchange each halo slot must hold its own global index,
  !> which makes the check exact rather than statistical.
  subroutine fill_with_global_index()
    integer :: i
    f = -1.0_dp
    do i = 1, m%nowned
       f(i) = real(m%owned_gidx(i), dp)
    end do
  end subroutine fill_with_global_index

  logical function halo_is_correct() result(ok)
    integer :: s
    ok = .true.
    do s = 1, m%nhalo
       if (abs(f(m%nowned + s) - real(m%halo_gidx(s), dp)) > 0.0_dp) then
          ok = .false.
          return
       end if
    end do
    ! A zero-length halo is not a pass — it means TODO 1 is still a stub.
    if (m%nhalo == 0) ok = .false.
  end function halo_is_correct

  subroutine print_header()
    integer :: e
    call MPI_Reduce(m%nhalo,  max_halo,  1, MPI_INTEGER, MPI_MAX, 0, comm, e)
    call MPI_Reduce(m%nneigh, max_neigh, 1, MPI_INTEGER, MPI_MAX, 0, comm, e)
    call MPI_Reduce(m%nhalo*8, tot_bytes, 1, MPI_INTEGER, MPI_SUM, 0, comm, e)
    if (rank == 0) then
       print '(a)', '=== Exercise 06: unstructured halo exchange ==='
       print '(a,i0,a,i0,a,i0)', 'mesh ', nx, ' x ', ny, ' = ', m%nglobal
       print '(a,i0,a,i0)', 'ranks ', nprocs, '   cells/rank ~ ', m%nglobal/nprocs
       print '(a,i0)', 'max halo cells on any rank : ', max_halo
       print '(a,i0)', 'max neighbour ranks        : ', max_neigh
       print '(a,i0)', 'bytes exchanged per step   : ', tot_bytes
       print '(a)', ''
    end if
  end subroutine print_header

  subroutine report(label, t)
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: t
    logical  :: ok, all_ok
    real(dp) :: tmax
    integer  :: e
    ok = halo_is_correct()
    call MPI_Allreduce(ok, all_ok, 1, MPI_LOGICAL, MPI_LAND, comm, e)
    call MPI_Reduce(t, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, e)
    if (rank == 0) then
       print '(a,a,a,f9.4,a,f9.2,a)', '  ', label, '  total ', tmax, &
            ' s   per-exchange ', tmax/real(n_iter,dp)*1.0e6_dp, ' us'
       if (all_ok) then
          print '(a)', '     HALO CHECK: PASS'
       else
          print '(a)', '     HALO CHECK: FAIL  <-- halo contents are wrong'
       end if
    end if
  end subroutine report

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) nx
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) ny
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg); read(arg,*) n_iter
    end if
  end subroutine read_cli

end program uhalo

! ===========================================================================
! TODO 6 — write your answers in notes/day2.md
!
! (a) With the linear partition, how does the halo size scale with rank
!     count for a fixed 400x400 mesh? Derive it: a rank owns a horizontal
!     strip, so its halo is 2*nx regardless of how many rows it owns.
!     What is the surface-to-volume ratio, and what does that predict about
!     strong scaling? Verify against `make strong`.
!
! (b) Run `make shape`: 400x400 vs 1600x100 vs 100x1600 at fixed rank count.
!     The cell count is identical; the halo is not. Explain the difference
!     and say what it implies about how ICON should partition the sphere.
!
! (c) Which exchange won, p2p or Neighbor_alltoallv? Now raise the
!     neighbour count (many ranks, small mesh) and re-measure. At what
!     neighbour count does the collective start to win, and why?
!
! (d) m%g2l and m%owner are both O(nglobal) on EVERY rank. At 1e8 cells and
!     1000 ranks that is fatal. Name two replacements and give their
!     lookup complexity. (Hint: ICON keeps a sorted local list; think about
!     what a distributed directory buys you.)
!
! (e) The exchange currently moves one field. A dycore moves ~30 per step.
!     Restructure the schedule to exchange all 30 in ONE set of messages.
!     How much does that cut total time, and what is the tradeoff? Relate
!     your answer to the message-rate limit of a real interconnect
!     (~1e6-1e7 messages/s/node).
!
! (f) invert_to_send uses MPI_Alltoall over ALL ranks, which is O(P) memory
!     and O(P) time per rank at setup. At P = 100000 that alone is a
!     problem. Look up how a "sparse data exchange" / NBX algorithm avoids
!     it, using MPI_Ibarrier. Sketch it in three lines.
! ===========================================================================
