! ===========================================================================
! EXERCISE 23 — CAPSTONE A
! Add an ice sheet component without compromising the execution time
! ===========================================================================
!
! THE BRIEF
!   You have a working coupled model: atmosphere and ocean, running
!   concurrently on disjoint rank sets, exchanging fields every step
!   (Exercise 09). A glaciologist wants to add an ice sheet model.
!
!   The ice sheet:
!     - is expensive: one call costs ~8x an atmosphere step
!     - is infrequent: it only needs to run every 100 atmosphere steps
!     - consumes the TIME-AVERAGED surface mass balance over that window
!       (Exercise 12 -- instantaneous sampling is not acceptable)
!     - returns an updated ice-sheet geometry the atmosphere needs
!
!   Naively bolted on, this creates a periodic stall: every 100 steps the
!   whole coupled system stops while the ice sheet runs. The average cost
!   looks small (8/100 = 8%) but the machine is idle for all of it, so what
!   you actually pay is 8% of EVERY core, and the model's throughput --
!   simulated years per wall day, the number scientists care about -- drops
!   by more than that once you include the synchronisation.
!
!   Your job: get the ice sheet in without paying that.
!
! WHY (DKRZ)
!   This is the posting's headline duty, made concrete:
!
!     "extending [coupled Earth System models] with new components
!      (e.g. ice sheet models) without compromising the efficient
!      execution time"
!
!   If you prepare one thing for this interview, prepare this. Be able to
!   draw the timeline, name the three strategies, and quote the numbers you
!   measured for each.
!
! THE THREE STRATEGIES
!
!   A  SYNCHRONOUS (the baseline, and the problem)
!      Ice sheet lives on the atmosphere's ranks. Every 100 steps the
!      atmosphere stops, runs the ice sheet, and resumes. Ocean waits too,
!      because the next coupling point cannot happen until the atmosphere
!      gets there. Everything stalls.
!
!   B  CONCURRENT, LAGGED (usually the right answer)
!      Ice sheet gets its OWN ranks and runs continuously alongside the
!      others. At each coupling window it receives the accumulated forcing
!      non-blockingly, computes while the atmosphere and ocean carry on, and
!      delivers its result at the NEXT window. The atmosphere consumes an
!      ice geometry that is one window old.
!
!      Is that acceptable? Physically, yes: ice responds over centuries, so
!      a lag of one coupling window is far below the timescale of anything
!      it does. Say that explicitly -- the argument is what makes the
!      engineering choice legitimate rather than a corner cut.
!
!   C  CONCURRENT + REBALANCED
!      As B, but tune how many ranks the ice sheet gets. Too few and it
!      cannot finish inside a window (it falls behind, and the lag grows
!      without bound -- detect this!). Too many and you have stolen ranks
!      from the atmosphere for a component that then sits idle.
!
! TASKS
!   TODO 1  strategy A: synchronous ice sheet, measure the stall
!   TODO 2  strategy B: give it its own ranks; non-blocking, lagged exchange
!   TODO 3  detect the failure mode: ice sheet cannot keep up
!   TODO 4  strategy C: sweep the rank split, find the optimum
!   TODO 5  report simulated-years-per-wall-day for all three
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - all three strategies transfer the SAME total accumulated forcing
!     (conservation across the coupling interface -- verified, not assumed)
!   - B measurably beats A on throughput; you can quote the factor
!   - you can draw the coupled timeline for A and for B from memory
!   - you have a rank split from C and a one-line justification for it
!   - your code DETECTS and reports the "ice sheet fell behind" case rather
!     than silently growing the lag
!
! HINTS
!   - Reuse: Exercise 09 for the communicator split, Exercise 12 for the
!     accumulator, Exercise 03 for the non-blocking overlap idea.
!   - The key primitive for B is MPI_Irecv posted EARLY: the ice sheet posts
!     its receive for the next window's forcing before it starts computing
!     on the current one. Then the atmosphere's send never blocks.
!   - MPI_Test (not MPI_Wait) is how the atmosphere checks whether the ice
!     result has arrived without stalling. If it has not, use last window's
!     value and carry on -- and count that event.
!   - Throughput, not step time, is the metric. Report simulated years per
!     wall-clock day. That is what appears in a project proposal.
! ===========================================================================

module coupled_mod
  use mpi_f08
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: COMP_ATM = 0
  integer, parameter :: COMP_OCE = 1
  integer, parameter :: COMP_ICE = 2

  type :: comp_t
     integer :: id = -1
     character(len=12) :: name = ''
     type(MPI_Comm) :: comm = MPI_COMM_NULL
     integer :: rank = -1, size = 0
     integer :: world_rank = -1, world_size = 0
     ! world ranks of each component's root, for cross-component messages
     integer :: root_atm = 0, root_oce = 0, root_ice = 0
     real(dp) :: t_compute = 0.0_dp, t_wait = 0.0_dp, t_stall = 0.0_dp
  end type comp_t

  !> Time-accumulator, as Exercise 12. The ice sheet must consume the
  !> WINDOW AVERAGE of the surface mass balance, never an instantaneous
  !> sample -- otherwise it misses every melt event between coupling points.
  type :: accum_t
     real(dp) :: sum = 0.0_dp, dt_total = 0.0_dp
  end type accum_t

contains

  !> Split the world three ways. Layout:
  !>   [0, n_atm)                  atmosphere
  !>   [n_atm, n_atm+n_oce)        ocean
  !>   [n_atm+n_oce, world_size)   ice sheet   (0 ranks => strategy A)
  subroutine comp_init(c, n_atm, n_oce)
    type(comp_t), intent(out) :: c
    integer,      intent(in)  :: n_atm, n_oce
    integer :: ierr, colour

    call MPI_Comm_rank(MPI_COMM_WORLD, c%world_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, c%world_size, ierr)

    if (c%world_rank < n_atm) then
       c%id = COMP_ATM; c%name = 'atmosphere'
    else if (c%world_rank < n_atm + n_oce) then
       c%id = COMP_OCE; c%name = 'ocean'
    else
       c%id = COMP_ICE; c%name = 'ice_sheet'
    end if
    colour = c%id

    call MPI_Comm_split(MPI_COMM_WORLD, colour, c%world_rank, c%comm, ierr)
    call MPI_Comm_rank(c%comm, c%rank, ierr)
    call MPI_Comm_size(c%comm, c%size, ierr)

    c%root_atm = 0
    c%root_oce = n_atm
    c%root_ice = n_atm + n_oce
  end subroutine comp_init

  subroutine accumulate(a, flux, dt)
    type(accum_t), intent(inout) :: a
    real(dp),      intent(in)    :: flux, dt
    a%sum      = a%sum + flux * dt
    a%dt_total = a%dt_total + dt
  end subroutine accumulate

  function window_average(a) result(avg)
    type(accum_t), intent(inout) :: a
    real(dp) :: avg
    if (a%dt_total > 0.0_dp) then
       avg = a%sum / a%dt_total
    else
       avg = 0.0_dp
    end if
    a%sum = 0.0_dp; a%dt_total = 0.0_dp
  end function window_average

  !> Synthetic component cost. Relative weights chosen to make the problem
  !> real: the ice sheet is 8x an atmosphere step, but runs 100x less often.
  subroutine burn(units)
    real(dp), intent(in) :: units
    real(dp) :: acc
    integer  :: i, n
    n = int(units * 40000.0_dp)
    acc = 0.0_dp
    do i = 1, n
       acc = acc + sqrt(real(i, dp))
    end do
    if (acc < 0.0_dp) print *, acc
  end subroutine burn

end module coupled_mod


program icesheet
  use coupled_mod
  use mpi_f08
  implicit none

  type(comp_t)  :: c
  type(accum_t) :: acc_ice
  integer :: ierr
  integer :: nsteps   = 400        ! atmosphere steps
  integer :: ice_every = 100       ! ice sheet coupling window
  integer :: n_atm = -1, n_oce = -1, n_ice = -1
  character(len=12) :: strategy = 'A'

  real(dp) :: t_total, forcing_sent, forcing_check
  integer  :: n_late                ! times the ice result was not ready

  call MPI_Init(ierr)
  call read_cli()
  call decide_split()
  call comp_init(c, n_atm, n_oce)

  call banner()

  select case (trim(strategy))
  case ('A'); call run_synchronous(t_total)
  case ('B'); call run_concurrent(t_total)
  case ('C'); call run_concurrent(t_total)     ! same code, different split
  case default
     if (c%world_rank == 0) print '(a)', 'unknown strategy (use A, B or C)'
     call MPI_Finalize(ierr); stop
  end select

  call report(t_total)

  call MPI_Finalize(ierr)

contains

  subroutine decide_split()
    integer :: w, e
    call MPI_Comm_size(MPI_COMM_WORLD, w, e)
    if (trim(strategy) == 'A') then
       ! Strategy A: no dedicated ice ranks. The ice sheet runs INSIDE the
       ! atmosphere's timeloop, on the atmosphere's ranks.
       if (n_atm < 0) n_atm = max(1, w / 2)
       n_oce = w - n_atm
       n_ice = 0
    else
       ! Strategies B and C: the ice sheet gets its own ranks.
       if (n_ice < 0) n_ice = max(1, w / 8)
       if (n_atm < 0) n_atm = max(1, (w - n_ice) / 2)
       n_oce = w - n_atm - n_ice
    end if
    if (n_oce < 1) then
       n_oce = 1; n_atm = max(1, w - n_oce - n_ice)
    end if
  end subroutine decide_split

  ! ------------------------------------------------------------------------
  ! TODO 1: STRATEGY A — synchronous. The baseline you must beat.
  !
  ! Structure (all components loop together):
  !
  !   do step = 1, nsteps
  !      atmosphere: burn(3.0)      ocean: burn(1.0)
  !      accumulate the surface mass balance into acc_ice
  !      exchange atm <-> oce  (as Exercise 09)
  !
  !      if (mod(step, ice_every) == 0) then
  !         ! EVERYTHING STOPS HERE
  !         avg = window_average(acc_ice)
  !         atmosphere ranks: burn(8.0 * 3.0)     ! the ice sheet
  !         barrier across the whole world        ! ocean waits too
  !      end if
  !   end do
  !
  ! Time the ice-sheet block separately into c%t_stall. That number is the
  ! thing you are trying to eliminate, and you need it to prove you did.
  ! ------------------------------------------------------------------------
  subroutine run_synchronous(t_out)
    real(dp), intent(out) :: t_out
    integer  :: step, e
    real(dp) :: t0, tw, smb

    forcing_sent = 0.0_dp
    n_late = 0
    call MPI_Barrier(MPI_COMM_WORLD, e)
    t0 = MPI_Wtime()

    do step = 1, nsteps
       ! TODO 1a: component compute (atm 3.0 units, oce 1.0 units)
       ! TODO 1b: accumulate the surface mass balance
       smb = surface_mass_balance(step)
       call accumulate(acc_ice, smb, 1.0_dp)

       ! TODO 1c: atm <-> oce exchange

       if (mod(step, ice_every) == 0) then
          tw = MPI_Wtime()
          ! TODO 1d: window average, run the ice sheet on the atm ranks,
          !          then a WORLD barrier so the ocean waits too.
          forcing_sent = forcing_sent + window_average(acc_ice) * real(ice_every, dp)
          call MPI_Barrier(MPI_COMM_WORLD, e)
          c%t_stall = c%t_stall + (MPI_Wtime() - tw)
       end if
    end do

    t_out = MPI_Wtime() - t0
  end subroutine run_synchronous

  ! ------------------------------------------------------------------------
  ! TODO 2: STRATEGY B — concurrent and lagged. The fix.
  !
  ! Now the three components run genuinely side by side:
  !
  !   ATMOSPHERE + OCEAN ranks:
  !     do step = 1, nsteps
  !        compute; accumulate smb; exchange atm<->oce
  !        if (window closes) then
  !           avg = window_average(acc_ice)
  !           MPI_Isend avg to the ice root        ! does NOT block
  !           MPI_Test the pending ice result:
  !              arrived  -> adopt it
  !              not yet  -> keep the previous one, n_late = n_late + 1
  !           MPI_Irecv for the NEXT result
  !        end if
  !     end do
  !
  !   ICE SHEET ranks:
  !     do window = 1, nsteps / ice_every
  !        MPI_Recv the accumulated forcing         ! blocks -- fine, it has
  !                                                 ! nothing else to do
  !        burn(8.0 * 3.0)                          ! runs CONCURRENTLY with
  !                                                 ! atm+oce advancing
  !        MPI_Isend the result back
  !     end do
  !
  ! The atmosphere never waits for the ice sheet. The ice sheet's answer
  ! arrives one window late, which is physically fine for ice.
  !
  ! TODO 3: count n_late. If it is nonzero EVERY window, the ice sheet
  ! cannot finish inside a window and the lag is growing without bound --
  ! that is a real failure, not a lag. Report it loudly.
  ! ------------------------------------------------------------------------
  subroutine run_concurrent(t_out)
    real(dp), intent(out) :: t_out
    integer  :: step, e
    real(dp) :: t0, smb

    forcing_sent = 0.0_dp
    n_late = 0
    call MPI_Barrier(MPI_COMM_WORLD, e)
    t0 = MPI_Wtime()

    if (c%id == COMP_ICE) then
       ! TODO 2b: the ice sheet's own loop -- recv, compute, isend.
       block
         integer :: w, nwin
         nwin = nsteps / ice_every
         do w = 1, nwin
            ! TODO 2b: MPI_Recv forcing from the atmosphere root,
            !          burn(24.0), MPI_Isend the result back.
            call burn(0.0_dp)
         end do
       end block
    else
       do step = 1, nsteps
          ! TODO 2a: atm/oce compute and exchange, as strategy A
          smb = surface_mass_balance(step)
          call accumulate(acc_ice, smb, 1.0_dp)

          if (mod(step, ice_every) == 0) then
             ! TODO 2a: Isend the window average to the ice root;
             !          MPI_Test the outstanding result; post the next Irecv.
             forcing_sent = forcing_sent + &
                  window_average(acc_ice) * real(ice_every, dp)
          end if
       end do
    end if

    call MPI_Barrier(MPI_COMM_WORLD, e)
    t_out = MPI_Wtime() - t0
  end subroutine run_concurrent

  !> The forcing the ice sheet consumes: a seasonal cycle plus melt events.
  !> Sharp enough that instantaneous sampling would give the wrong answer,
  !> which is why the accumulator exists.
  pure function surface_mass_balance(step) result(f)
    integer, intent(in) :: step
    real(dp) :: f, t
    t = real(step, dp)
    f = 1.0_dp + 0.5_dp * sin(2.0_dp * 3.14159265358979_dp * t / 50.0_dp)
    if (modulo(step, 137) < 8) f = f + 4.0_dp        ! a melt event
  end function surface_mass_balance

  subroutine banner()
    integer :: e
    call MPI_Barrier(MPI_COMM_WORLD, e)
    if (c%world_rank == 0) then
       print '(a)', '=== Exercise 23 (CAPSTONE A): adding an ice sheet ==='
       print '(a,a)',  'strategy      : ', trim(strategy)
       print '(a,i0,a,i0,a,i0,a,i0)', 'ranks         : world ', c%world_size, &
            '  = atm ', n_atm, ' + oce ', n_oce, ' + ice ', n_ice
       print '(a,i0,a,i0)', 'steps         : ', nsteps, &
            '   ice window every ', ice_every
       print '(a)', 'costs         : atm 3.0  oce 1.0  ice 24.0 units per call'
       print '(a)', ''
    end if
    call MPI_Barrier(MPI_COMM_WORLD, e)
  end subroutine banner

  ! ------------------------------------------------------------------------
  ! TODO 5: the report. Throughput is the metric that matters.
  ! ------------------------------------------------------------------------
  subroutine report(t)
    real(dp), intent(in) :: t
    real(dp) :: tmax, stall_max, sent_tot
    integer  :: late_tot, e

    call MPI_Reduce(t,           tmax,      1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, e)
    call MPI_Reduce(c%t_stall,   stall_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, e)
    call MPI_Reduce(forcing_sent, sent_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, e)
    call MPI_Reduce(n_late,      late_tot,  1, MPI_INTEGER,          MPI_SUM, 0, MPI_COMM_WORLD, e)

    if (c%world_rank == 0) then
       print '(a,f10.4,a)', '  coupled wall time      : ', tmax, ' s'
       print '(a,f10.4,a)', '  time stalled on ice    : ', stall_max, ' s'
       if (tmax > 0.0_dp) &
            print '(a,f10.2,a)', '  stall fraction         : ', &
            100.0_dp * stall_max / tmax, ' %'
       ! es, not f: with the kernels still stubbed the wall time is tiny and
       ! a fixed-width field just prints asterisks.
       print '(a,es12.4)',  '  steps per second       : ', &
            real(nsteps, dp) / max(tmax, 1.0e-12_dp)
       print '(a)',         '  simulated yr / wall d  : TODO 5'
       print '(a,es16.8)',  '  total forcing to ice   : ', sent_tot
       print '(a,i0,a,i0)', '  windows the ice was late: ', late_tot, &
            ' of ', nsteps / ice_every
       if (late_tot >= nsteps / ice_every .and. trim(strategy) /= 'A') then
          print '(a)', '  *** ICE SHEET NEVER KEEPS UP -- the lag is growing.'
          print '(a)', '  *** Give it more ranks, or lengthen the window.'
       end if
       print '(a)', ''
       print '(a)', '  Compare A and B on wall time AND on total forcing.'
       print '(a)', '  If the forcing differs, B is not conservative and the'
       print '(a)', '  speedup is worthless -- check that first, every time.'
    end if
  end subroutine report

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); strategy = trim(arg)
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) nsteps
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg); read(arg,*) n_ice
    end if
    if (command_argument_count() >= 4) then
       call get_command_argument(4, arg); read(arg,*) n_atm
    end if
  end subroutine read_cli

end program icesheet

! ===========================================================================
! TODO 6 — write your answers in notes/day6.md  (and INTERVIEW.md)
!
! (a) Draw the coupled timeline for strategy A and for strategy B: three
!     horizontal bars (atm, oce, ice), 300 steps, marking compute and idle.
!     This diagram is the single most useful thing you can put in front of
!     an interviewer for this job. Practise drawing it in 30 seconds.
!
! (b) Measured: what is B's speedup over A? Now predict it from first
!     principles -- ice cost 24 units every 100 steps against an atmosphere
!     step of 3 units -- and reconcile any gap between prediction and
!     measurement.
!
! (c) Strategy C: sweep n_ice with `make split`. Plot throughput against
!     ice ranks. There are two failure modes at the two ends -- name them,
!     and identify the optimum. How would you find this optimum for a REAL
!     model where you cannot run the sweep cheaply?
!
! (d) The lag. Strategy B gives the atmosphere an ice geometry one window
!     old. Justify that physically for ice. Then name a component where the
!     same trick would NOT be acceptable, and say what makes the difference.
!
! (e) The failure mode in TODO 3: if the ice sheet cannot finish inside a
!     window, the lag grows every window. How would you DETECT that in a
!     production run that lasts three months, and what should the model do
!     about it -- abort, warn, or adapt? Defend your choice.
!
! (f) Conservation: strategy B sends the forcing non-blockingly and the ice
!     sheet may still be working when the next window closes. What happens
!     to the second window's forcing? Design the buffering so nothing is
!     lost or double-counted, and say how you would TEST that.
!
! (g) Now the real-world version. A glaciologist wants to couple PISM (a
!     real ice-sheet model, ~100k lines, its own grid, its own MPI setup)
!     into ICON via YAC. List the first five things you would do, in order.
!     This is very close to what the job actually is -- have an answer ready.
! ===========================================================================
