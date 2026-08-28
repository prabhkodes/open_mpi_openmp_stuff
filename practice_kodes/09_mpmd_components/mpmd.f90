! ===========================================================================
! EXERCISE 09 — Two components, one MPI_COMM_WORLD: the coupled-model skeleton
! ===========================================================================
!
! GOAL
!   Split the world communicator into an "atmosphere" and an "ocean" running
!   CONCURRENTLY on disjoint rank sets, build an intercommunicator between
!   them, and exchange fields every coupling step. Then measure the thing
!   that actually matters: how long each component spends WAITING for the
!   other.
!
! WHY (DKRZ)
!   This is the architecture of every coupled Earth System Model, and it is
!   the structure the job description is describing when it says "coupled
!   Earth System models" and "extending models with new components".
!
!   Concurrent (not sequential) execution is the whole point: atmosphere and
!   ocean advance at the same wall-clock time on different cores, and only
!   synchronise at coupling points. That is also the whole problem — the
!   coupled timeline runs at the speed of whichever component is slower, so
!   a component that takes 10% longer wastes 10% of the OTHER component's
!   cores for the entire run.
!
!   The number you produce here (idle time per component) is exactly the
!   number you would put in front of a scientist to argue for a different
!   rank split. Exercise 23 makes you fix it.
!
! TASKS
!   TODO 1  split the world by colour into component communicators
!   TODO 2  build the intercommunicator
!   TODO 3  exchange fields across it
!   TODO 4  instrument the wait time
!   TODO 5  find the rank split that minimises total idle time
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - `make run NP=8` prints PASS for the field exchange in both directions
!   - a table of: component, ranks, compute time, wait time, % idle
!   - you can state the optimal atmosphere:ocean rank ratio for the built-in
!     cost model, and you derived it before measuring it
!
! HINTS
!   - MPI_Comm_split(comm, colour, key, newcomm) — same colour lands in the
!     same communicator. This is the ONLY thing an MPMD launch really needs.
!   - MPI_Intercomm_create(local_comm, local_leader, peer_comm, remote_leader,
!     tag, newintercomm) — note local_leader is a rank in LOCAL_COMM, and
!     remote_leader is a rank in PEER_COMM (usually MPI_COMM_WORLD). Mixing
!     up which communicator each rank number refers to is the classic bug.
!   - On an intercommunicator, rank N means "rank N in the REMOTE group".
!     MPI_Comm_size gives your local size; MPI_Comm_remote_size gives theirs.
!   - Real couplers do not exchange whole fields root-to-root; they build a
!     rank-to-rank schedule so every pair that shares grid overlap talks
!     directly. Start root-to-root, then read TODO 6c.
! ===========================================================================

module component_mod
  use mpi_f08
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: COMP_ATM = 0
  integer, parameter :: COMP_OCE = 1

  type :: component_t
     integer :: id = -1
     character(len=16) :: name = ''
     type(MPI_Comm) :: comm      = MPI_COMM_NULL   ! my component
     type(MPI_Comm) :: inter     = MPI_COMM_NULL   ! to the other component
     integer :: rank = -1, size = 0, remote_size = 0
     integer :: world_rank = -1, world_size = 0
     integer :: local_leader = 0, remote_leader = -1
     real(dp) :: t_compute = 0.0_dp
     real(dp) :: t_wait    = 0.0_dp
  end type component_t

contains

  ! ------------------------------------------------------------------------
  ! TODO 1 + TODO 2: build the component and intercommunicators.
  !
  ! Layout: world ranks [0, n_atm) are atmosphere, [n_atm, world_size) are
  ! ocean. So:
  !     colour = merge(COMP_ATM, COMP_OCE, world_rank < n_atm)
  !
  ! TODO 1: c%comm = MPI_Comm_split(MPI_COMM_WORLD, colour, world_rank)
  !         then fill c%rank and c%size from it.
  !
  ! TODO 2: the intercommunicator.
  !         Atmosphere's local leader is world rank 0        -> local_leader = 0
  !         Ocean's local leader is world rank n_atm         -> local_leader = 0
  !         Each side's REMOTE leader, expressed in MPI_COMM_WORLD:
  !             atmosphere's remote leader = n_atm
  !             ocean's      remote leader = 0
  !
  !         call MPI_Intercomm_create(c%comm, 0, MPI_COMM_WORLD,           &
  !                                   c%remote_leader, 99, c%inter, ierr)
  !
  !         Then MPI_Comm_remote_size(c%inter, c%remote_size, ierr).
  ! ------------------------------------------------------------------------
  subroutine component_init(c, n_atm)
    type(component_t), intent(out) :: c
    integer,           intent(in)  :: n_atm
    integer :: ierr, colour

    call MPI_Comm_rank(MPI_COMM_WORLD, c%world_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, c%world_size, ierr)

    if (c%world_rank < n_atm) then
       c%id = COMP_ATM; c%name = 'atmosphere'; c%remote_leader = n_atm
    else
       c%id = COMP_OCE; c%name = 'ocean';      c%remote_leader = 0
    end if
    colour = c%id

    ! TODO 1: split the world.
    c%comm = MPI_COMM_WORLD                 ! placeholder -- wrong on purpose
    call MPI_Comm_rank(c%comm, c%rank, ierr)
    call MPI_Comm_size(c%comm, c%size, ierr)

    ! TODO 2: create the intercommunicator and set c%remote_size.
    c%inter       = MPI_COMM_NULL
    c%remote_size = 0

    if (colour < 0) continue
  end subroutine component_init

  ! ------------------------------------------------------------------------
  ! TODO 3: exchange a field across the intercommunicator.
  !
  ! Start simple, root-to-root:
  !   - each component gathers its field onto its own local root
  !   - the two roots MPI_Sendrecv across the intercommunicator
  !     (on an intercomm, destination 0 means "rank 0 of the OTHER group")
  !   - each root scatters/broadcasts the received field to its own ranks
  !
  ! `send` is what I give the other component, `recv` is what I get back.
  ! For this exercise both are single scalars per rank, gathered into a
  ! vector — enough to verify correctness without drowning in remapping,
  ! which is Exercise 10's job.
  !
  ! Time this whole routine into c%t_wait: from the caller's point of view,
  ! everything here that is not local computation is coupling overhead, and
  ! most of it will turn out to be waiting for the other component to arrive.
  ! ------------------------------------------------------------------------
  subroutine couple_exchange(c, send, recv)
    type(component_t), intent(inout) :: c
    real(dp),          intent(in)    :: send(:)
    real(dp),          intent(out)   :: recv(:)
    integer :: ierr

    ! TODO 3: implement the exchange. Suggested structure:
    !   real(dp), allocatable :: gsend(:), grecv(:)
    !   allocate(gsend(c%size), grecv(c%remote_size))
    !   call MPI_Gather(send, 1, MPI_DOUBLE_PRECISION, gsend, 1, ..., 0, c%comm)
    !   if (c%rank == 0) then
    !      call MPI_Sendrecv(gsend, c%size,        MPI_DOUBLE_PRECISION, 0, 7, &
    !                        grecv, c%remote_size, MPI_DOUBLE_PRECISION, 0, 7, &
    !                        c%inter, MPI_STATUS_IGNORE, ierr)
    !   end if
    !   call MPI_Bcast(grecv, c%remote_size, MPI_DOUBLE_PRECISION, 0, c%comm)
    !   then map grecv into recv (see the verification below for the rule)

    recv = -1.0_dp     ! placeholder: nothing received
    ierr = 0
    if (size(send) < 0) continue
  end subroutine couple_exchange

  !> Synthetic per-step cost. The atmosphere is deliberately more expensive
  !> per rank than the ocean — as it is in reality, where the atmosphere
  !> carries radiation and microphysics. This is what makes the rank split
  !> in TODO 5 a real optimisation problem.
  subroutine do_work(c, work_units)
    type(component_t), intent(inout) :: c
    integer,           intent(in)    :: work_units
    real(dp) :: t0, acc
    integer  :: i, n
    t0 = MPI_Wtime()
    n = work_units
    if (c%id == COMP_ATM) n = n * 3       ! atmosphere costs 3x per rank
    acc = 0.0_dp
    do i = 1, n * 20000
       acc = acc + sqrt(real(i, dp))
    end do
    if (acc < 0.0_dp) print *, acc        ! keep the compiler honest
    c%t_compute = c%t_compute + (MPI_Wtime() - t0)
  end subroutine do_work

end module component_mod


program mpmd
  use component_mod
  use mpi_f08
  implicit none

  type(component_t) :: c
  integer :: ierr, n_atm = -1, nsteps = 20, work = 5
  integer :: i, npass, ntest
  real(dp), allocatable :: send(:), recv(:)
  real(dp) :: t0, t_total

  call MPI_Init(ierr)
  call read_cli()

  block
    integer :: wsize
    call MPI_Comm_size(MPI_COMM_WORLD, wsize, ierr)
    if (n_atm < 0) n_atm = max(1, wsize / 2)
    if (n_atm >= wsize) n_atm = max(1, wsize - 1)
  end block

  call component_init(c, n_atm)

  allocate(send(1), recv(max(c%remote_size, 1)))
  npass = 0; ntest = 0

  call print_header()

  ! ---- the coupled timeloop ----------------------------------------------
  t0 = MPI_Wtime()
  do i = 1, nsteps
     call do_work(c, work)

     ! Each rank contributes a value the other side can verify exactly:
     !   atmosphere rank r sends  1000 + r
     !   ocean      rank r sends  2000 + r
     send(1) = real(merge(1000, 2000, c%id == COMP_ATM) + c%rank, dp)

     block
       real(dp) :: tw
       tw = MPI_Wtime()
       call couple_exchange(c, send, recv)
       c%t_wait = c%t_wait + (MPI_Wtime() - tw)
     end block
  end do
  t_total = MPI_Wtime() - t0

  ! ---- verification ------------------------------------------------------
  call check_exchange()
  call report(t_total)

  deallocate(send, recv)
  call MPI_Finalize(ierr)

contains

  subroutine print_header()
    integer :: e
    call MPI_Barrier(MPI_COMM_WORLD, e)
    if (c%world_rank == 0) then
       print '(a)', '=== Exercise 09: concurrent coupled components ==='
       print '(a,i0,a,i0,a,i0)', 'world ', c%world_size, ' ranks -> atm ', &
            n_atm, ' + oce ', c%world_size - n_atm
       print '(a,i0,a,i0)', 'coupling steps ', nsteps, '   work units ', work
       print '(a)', '(atmosphere is 3x more expensive per rank -- as in reality)'
       print '(a)', ''
    end if
    call MPI_Barrier(MPI_COMM_WORLD, e)
  end subroutine print_header

  !> Every rank must have received the OTHER component's values.
  subroutine check_exchange()
    integer :: r, expect_base
    logical :: ok
    expect_base = merge(2000, 1000, c%id == COMP_ATM)
    ok = (c%remote_size > 0)
    do r = 1, c%remote_size
       if (abs(recv(r) - real(expect_base + r - 1, dp)) > 0.0_dp) ok = .false.
    end do
    call tally('field exchange', ok)
  end subroutine check_exchange

  subroutine tally(label, ok)
    character(len=*), intent(in) :: label
    logical,          intent(in) :: ok
    logical :: all_ok
    integer :: e
    ntest = ntest + 1
    call MPI_Allreduce(ok, all_ok, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, e)
    if (all_ok) npass = npass + 1
    if (c%world_rank == 0) print '(a,a24,a)', '  ', label, &
         merge('  PASS', '  FAIL', all_ok)
  end subroutine tally

  ! ------------------------------------------------------------------------
  ! TODO 4: the idle-time report. This is the deliverable.
  !
  ! For each component print: ranks, mean compute time, mean wait time, and
  ! wait as a percentage of total. Then print the wasted core-seconds:
  !
  !     wasted = sum over components of (ranks * wait_time)
  !
  ! That number is what you would show a scientist to justify changing the
  ! rank split. On a real machine it is money.
  ! ------------------------------------------------------------------------
  subroutine report(t_total)
    real(dp), intent(in) :: t_total
    real(dp) :: comp_mean, wait_mean, tmax
    integer  :: e

    call MPI_Reduce(c%t_compute, comp_mean, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, c%comm, e)
    call MPI_Reduce(c%t_wait,    wait_mean, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, c%comm, e)
    if (c%rank == 0 .and. c%size > 0) then
       comp_mean = comp_mean / real(c%size, dp)
       wait_mean = wait_mean / real(c%size, dp)
    end if
    call MPI_Reduce(t_total, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                    MPI_COMM_WORLD, e)

    ! Both component roots print; serialise them so the output is readable.
    block
      integer :: turn
      do turn = 0, 1
         call MPI_Barrier(MPI_COMM_WORLD, e)
         if (c%rank == 0 .and. c%id == turn) then
            print '(a,a12,i6,a,f9.4,a,f9.4,a,f6.1,a)',                  &
                 '  ', trim(c%name), c%size, ' ranks   compute ',       &
                 comp_mean, ' s   wait ', wait_mean, ' s   idle ',       &
                 100.0_dp * wait_mean / max(comp_mean + wait_mean, 1.0e-12_dp), &
                 ' %'
         end if
      end do
      call MPI_Barrier(MPI_COMM_WORLD, e)
    end block

    if (c%world_rank == 0) then
       print '(a)', ''
       print '(a,f9.4,a)', '  coupled wall time : ', tmax, ' s'
       print '(a)',        '  wasted core-sec   : TODO 4'
       print '(a,i0,a,i0,a)', '  ', npass, ' / ', ntest, ' checks passed'
       print '(a)', ''
       print '(a)', '  Now run  make split  and find the rank ratio that'
       print '(a)', '  minimises wasted core-seconds. Predict it first.'
    end if
  end subroutine report

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) n_atm
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) nsteps
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg); read(arg,*) work
    end if
  end subroutine read_cli

end program mpmd

! ===========================================================================
! TODO 6 — write your answers in notes/day3.md
!
! (a) The atmosphere costs 3 work units per rank per step, the ocean 1. With
!     W total ranks, derive the split n_atm that equalises the two
!     components' step times. Then run `make split` and check. Does the
!     measured optimum match? If not, what is the extra cost you ignored?
!
! (b) At the optimal split, what fraction of the machine is still idle?
!     (It is not zero. Explain where it goes.)
!
! (c) Root-to-root exchange serialises all coupling traffic through two
!     ranks. Sketch what a real coupler does instead, and estimate the
!     speedup for a field of 1e6 values across 100 + 100 ranks. What extra
!     information does the coupler need at setup time to do that? (You build
!     that information in Exercise 10.)
!
! (d) SEQUENTIAL vs CONCURRENT coupling: some models run atmosphere and
!     ocean one after the other on ALL ranks instead of side by side.
!     Give one advantage of each. Which one makes the ice-sheet problem of
!     Exercise 23 easier, and which makes it harder?
!
! (e) MPI_Intercomm_create needs a "peer communicator" containing both
!     leaders. Why can that not just be the local communicator? What
!     happens if two components try to create intercommunicators with the
!     same tag at the same time?
!
! (f) DKRZ context: YAC does not use intercommunicators at all — it builds
!     its own rank-to-rank schedule over MPI_COMM_WORLD from the grid
!     overlap. Give one reason to prefer that over the intercommunicator
!     approach you just built.
! ===========================================================================
