! ===========================================================================
! EXERCISE 12 — Multi-rate coupling: the ice-sheet problem
! ===========================================================================
!
! GOAL
!   Couple three components whose timesteps differ by two orders of
!   magnitude — atmosphere dt=1, ocean dt=4, ice sheet dt=100 — and keep the
!   coupled system conservative in TIME, not just in space. Then make it
!   restartable, which turns out to be the same problem wearing a hat.
!
! WHY (DKRZ)
!   The posting names this explicitly:
!
!     "extending [models] with new components (e.g. ICE SHEET MODELS)
!      without compromising the efficient execution time"
!
!   Ice sheets are the canonical hard case. Ice responds over centuries, so
!   its timestep is enormous compared to the atmosphere's. You cannot just
!   sample the atmosphere's surface mass balance once every 100 steps and
!   call it coupled — you would miss every melt event between samples, and
!   the ice sheet would grow or shrink for entirely numerical reasons.
!
!   What you must do instead is ACCUMULATE the flux over the whole coupling
!   window and hand over the time-average. That keeps the total mass and
!   energy transferred exactly right regardless of how the timesteps line up.
!
!   And then: the accumulator is model state. If it is not in the restart
!   file, a restarted run silently differs from a continuous one — which is
!   the single most common "why doesn't my restart reproduce?" bug in
!   coupled modelling, and a thing you will absolutely be asked to debug.
!
! TASKS
!   TODO 1  accumulate       — flux-weighted accumulation over a window
!   TODO 2  window_average   — hand over the average, reset the accumulator
!   TODO 3  naive_sample     — instantaneous sampling, for contrast
!   TODO 4  write_restart / read_restart — including the accumulators
!   TODO 5  lagged coupling for the ice sheet
!   TODO 6  answer the questions at the bottom
!
! ACCEPTANCE
!   - accumulated coupling transfers EXACTLY the same total flux as a
!     reference that exchanges every atmosphere step (error < 1e-13)
!   - naive instantaneous sampling shows a clear, non-vanishing error
!   - a run restarted at step N is BIT-IDENTICAL to an uninterrupted run
!   - you can state why the ice sheet must use the previous window's average
!     rather than the current one, and what that costs physically
!
! HINTS
!   - "Conservative in time" means: integral over the window of the flux,
!     divided by the window length. If all steps are the same length that is
!     a plain mean; if they are not, it is dt-weighted. Write the weighted
!     version — real models have variable timesteps.
!   - The accumulator must be reset at exactly the right moment. Off-by-one
!     here gives an error of one step out of the window, which for a 100-step
!     window is 1% — small enough to look like physics, big enough to ruin a
!     century-long run.
!   - For the restart test, dump every piece of state you think matters, then
!     diff. When it does not reproduce, the thing you forgot IS the answer.
!   - The lag is not a bug. The ice sheet cannot use the average of a window
!     that has not finished yet, so it necessarily runs one window behind.
!     Understand the difference between that and an accidental lag.
! ===========================================================================

module multirate_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  !> A time-accumulator for one coupled field. This little type is the whole
  !> idea of the exercise: it is what turns a fast component's rapidly
  !> varying flux into something a slow component can consume without losing
  !> mass or energy.
  type :: accumulator_t
     real(dp) :: sum      = 0.0_dp     ! integral of flux*dt over the window
     real(dp) :: dt_total = 0.0_dp     ! window length so far
     integer  :: nsamples = 0
  end type accumulator_t

  type :: component_t
     character(len=16) :: name = ''
     real(dp) :: dt      = 1.0_dp
     real(dp) :: time    = 0.0_dp
     real(dp) :: state   = 0.0_dp      ! whatever this component integrates
     real(dp) :: received = 0.0_dp     ! last value handed to it
     real(dp) :: total_received = 0.0_dp  ! time-integral of what it received
  end type component_t

contains

  !> The atmosphere's surface flux: strongly time-varying, with a diurnal
  !> cycle and a sharp melt event. If you sample this instantaneously every
  !> 100 steps you will catch the event or miss it depending purely on phase
  !> — which is exactly the failure this exercise is about.
  pure function surface_flux(t) result(f)
    real(dp), intent(in) :: t
    real(dp) :: f
    f = 1.0_dp + 0.8_dp * sin(2.0_dp * 3.14159265358979_dp * t / 24.0_dp)
    ! A melt event: brief, large, and easy to miss.
    if (t > 300.0_dp .and. t < 340.0_dp) f = f + 5.0_dp
  end function surface_flux

  ! ------------------------------------------------------------------------
  ! TODO 1: accumulate one sample into the window.
  !
  !   acc%sum      = acc%sum + flux * dt
  !   acc%dt_total = acc%dt_total + dt
  !   acc%nsamples = acc%nsamples + 1
  !
  ! Multiplying by dt is the entire point. If you accumulate the flux without
  ! weighting and then divide by nsamples, you get the right answer ONLY when
  ! every timestep is identical. Real models adapt their timestep, so write
  ! the weighted form from the start.
  ! ------------------------------------------------------------------------
  subroutine accumulate(acc, flux, dt)
    type(accumulator_t), intent(inout) :: acc
    real(dp),            intent(in)    :: flux, dt
    ! TODO 1
    if (flux < -huge(1.0_dp) .or. dt < 0.0_dp) continue
    if (acc%nsamples < 0) continue
  end subroutine accumulate

  ! ------------------------------------------------------------------------
  ! TODO 2: close the window — return the time-average and reset.
  !
  !   avg = acc%sum / acc%dt_total
  !   then zero the accumulator
  !
  ! Guard dt_total == 0 (a window in which nothing was accumulated). Decide
  ! what that should mean and be able to defend it: returning 0 and returning
  ! the previous value are both defensible, and they give different physics.
  ! ------------------------------------------------------------------------
  function window_average(acc) result(avg)
    type(accumulator_t), intent(inout) :: acc
    real(dp) :: avg
    avg = 0.0_dp
    ! TODO 2
    if (acc%nsamples < 0) continue
  end function window_average

  ! ------------------------------------------------------------------------
  ! TODO 3: the naive alternative — instantaneous sampling.
  !
  ! Just read the flux at the coupling instant and hand it over, ignoring
  ! everything that happened in between. This is what you get if you couple
  ! by "call the exchange every N steps" without thinking about time
  ! averaging, and it is a very easy mistake to make because it LOOKS right
  ! and runs fine.
  !
  ! Implement it so you can measure how wrong it is.
  ! ------------------------------------------------------------------------
  function naive_sample(t) result(f)
    real(dp), intent(in) :: t
    real(dp) :: f
    f = 0.0_dp
    ! TODO 3: return surface_flux(t)
    if (t < -1.0_dp) continue
  end function naive_sample

  ! ------------------------------------------------------------------------
  ! TODO 4: restart I/O.
  !
  ! Write EVERY piece of state needed to continue: each component's time,
  ! state, received value and running total, AND every accumulator's sum,
  ! dt_total and nsamples.
  !
  ! The accumulators are the ones people forget. They are not "the model
  ! state" in the physics sense, so they get left out of the restart file,
  ! and then a restarted run differs from a continuous one by up to one
  ! coupling window of flux. Small. Plausible. Extremely annoying to track
  ! down six months later.
  !
  ! Use unformatted stream I/O so the file is exact — a formatted write would
  ! round, and then "bit-identical" is impossible by construction.
  ! ------------------------------------------------------------------------
  subroutine write_restart(fname, comps, accs)
    character(len=*),    intent(in) :: fname
    type(component_t),   intent(in) :: comps(:)
    type(accumulator_t), intent(in) :: accs(:)
    integer :: u
    open(newunit=u, file=fname, form='unformatted', access='stream', &
         status='replace', action='write')
    ! TODO 4a: write comps and accs. Writing the derived types directly is
    ! not portable across compilers -- write the members explicitly.
    close(u)
    if (size(comps) < 0 .or. size(accs) < 0) continue
  end subroutine write_restart

  subroutine read_restart(fname, comps, accs)
    character(len=*),    intent(in)    :: fname
    type(component_t),   intent(inout) :: comps(:)
    type(accumulator_t), intent(inout) :: accs(:)
    integer :: u, ios
    open(newunit=u, file=fname, form='unformatted', access='stream', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) return
    ! TODO 4b: read back in EXACTLY the order written.
    close(u)
    if (size(comps) < 0 .or. size(accs) < 0) continue
  end subroutine read_restart

end module multirate_mod


program multirate
  use multirate_mod
  implicit none

  integer, parameter :: ATM = 1, OCE = 2, ICE = 3
  integer  :: nsteps = 1000
  integer  :: restart_at = 500
  real(dp) :: dt_atm = 1.0_dp, dt_oce = 4.0_dp, dt_ice = 100.0_dp

  type(component_t)   :: comps(3)
  type(accumulator_t) :: accs(2)          ! 1: atm->oce   2: atm->ice
  integer :: npass, ntest
  real(dp) :: ref_integral, acc_integral, naive_integral
  real(dp) :: cont_state, rst_state

  call read_cli()
  npass = 0; ntest = 0

  print '(a)', '=== Exercise 12: multi-rate coupling (the ice-sheet problem) ==='
  print '(a,f5.1,a,f5.1,a,f6.1)', 'dt: atmosphere ', dt_atm, &
       '   ocean ', dt_oce, '   ice sheet ', dt_ice
  print '(a,i0,a,f8.1)', 'steps ', nsteps, '   total time ', nsteps*dt_atm
  print '(a)', '(the flux has a diurnal cycle plus a melt event at t=300..340)'
  print '(a)', ''

  ! ---- the reference: exact time-integral of the flux --------------------
  ref_integral = reference_integral()
  print '(a,es20.12)', '  reference integral of flux : ', ref_integral

  ! ---- accumulated (correct) coupling ------------------------------------
  acc_integral = run_accumulated()
  call check('accumulated coupling', &
       abs(acc_integral - ref_integral) / abs(ref_integral), 1.0e-13_dp)

  ! ---- naive instantaneous sampling --------------------------------------
  naive_integral = run_naive()
  print '(a,es20.12,a,f8.3,a)', '  naive sampled integral     : ', &
       naive_integral, '   error ', &
       100.0_dp * abs(naive_integral - ref_integral)/abs(ref_integral), ' %'
  print '(a)', '     (expected to be WRONG -- that is the lesson)'

  ! ---- restart reproducibility -------------------------------------------
  print '(a)', ''
  print '(a)', '  --- restart reproducibility ---'
  cont_state = run_to(nsteps, .false., 0)
  rst_state  = run_to(nsteps, .true., restart_at)
  print '(a,es22.14)', '    continuous run final state : ', cont_state
  print '(a,es22.14)', '    restarted  run final state : ', rst_state
  if (cont_state == 0.0_dp) then
     ! Two zeros match perfectly, which would be a very misleading PASS.
     ! Until the ice sheet actually receives something there is nothing to
     ! reproduce, so this is a failure, not a success.
     ntest = ntest + 1
     print '(a)', '    restart is bit-identical      : FAIL (nothing was'
     print '(a)', '                                    transferred -- TODO 5)'
  else
     call check('restart is bit-identical', &
          abs(cont_state - rst_state), 0.0_dp)
  end if

  print '(a)', ''
  print '(a,i0,a,i0,a)', '  ', npass, ' / ', ntest, ' checks passed'
  if (npass < ntest) then
     print '(a)', '  Order: TODO 1 -> 2 -> 3 -> 4. The restart test will keep'
     print '(a)', '  failing until the ACCUMULATORS are in the restart file.'
  end if

contains

  !> Trapezoidal integral of the flux over the whole run, at atmosphere
  !> resolution. This is the answer any correct coupling must reproduce.
  real(dp) function reference_integral() result(s)
    integer :: i
    real(dp) :: t
    s = 0.0_dp
    do i = 1, nsteps
       t = real(i-1, dp) * dt_atm
       s = s + surface_flux(t) * dt_atm
    end do
  end function reference_integral

  !> Correct coupling: accumulate every atmosphere step, hand over the
  !> time-average when the ocean/ice window closes.
  real(dp) function run_accumulated() result(total)
    integer  :: i
    real(dp) :: t, f, avg
    type(accumulator_t) :: a_oce, a_ice

    a_oce = accumulator_t(); a_ice = accumulator_t()
    total = 0.0_dp

    do i = 1, nsteps
       t = real(i-1, dp) * dt_atm
       f = surface_flux(t)
       call accumulate(a_oce, f, dt_atm)
       call accumulate(a_ice, f, dt_atm)

       ! Ocean window closes every dt_oce.
       if (mod(i * int(dt_atm), int(dt_oce)) == 0) then
          avg = window_average(a_oce)
          total = total + avg * dt_oce
       end if

       ! Ice-sheet window closes every dt_ice. Note we do NOT add this to
       ! `total` -- the ocean already counts the flux once. Accumulating the
       ! same flux for two consumers is fine; double-counting it in a
       ! conservation budget is not. Be able to explain that distinction.
       if (mod(i * int(dt_atm), int(dt_ice)) == 0) then
          avg = window_average(a_ice)
       end if
    end do
  end function run_accumulated

  !> Naive coupling: sample the instantaneous flux at coupling time only.
  real(dp) function run_naive() result(total)
    integer  :: i
    real(dp) :: t
    total = 0.0_dp
    do i = 1, nsteps
       t = real(i-1, dp) * dt_atm
       if (mod(i * int(dt_atm), int(dt_oce)) == 0) then
          total = total + naive_sample(t) * dt_oce
       end if
    end do
  end function run_naive

  ! ------------------------------------------------------------------------
  ! TODO 5: the run driver, with an optional restart in the middle.
  !
  ! When `do_restart` is true: run to `stop_at`, write the restart, zero
  ! everything, read the restart back, and continue to nsteps. The final
  ! state must be bit-identical to the uninterrupted run.
  !
  ! The ice sheet uses the PREVIOUS window's average (TODO 5): it cannot use
  ! the current window, because that window has not finished. Hold the last
  ! completed average in a variable and apply it throughout the next window.
  ! That one-window lag is physically correct and unavoidable — make sure
  ! you can explain the difference between this and an accidental lag.
  ! ------------------------------------------------------------------------
  real(dp) function run_to(n, do_restart, stop_at) result(final_state)
    integer, intent(in) :: n, stop_at
    logical, intent(in) :: do_restart
    integer  :: i, i0
    real(dp) :: t, f, ice_lagged

    comps(ATM) = component_t('atmosphere', dt_atm, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    comps(OCE) = component_t('ocean',      dt_oce, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    comps(ICE) = component_t('ice_sheet',  dt_ice, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    accs = accumulator_t()
    ice_lagged = 0.0_dp
    i0 = 1

    if (do_restart) then
       call step_range(1, stop_at, ice_lagged)
       call write_restart('restart.bin', comps, accs)
       ! Wipe everything, exactly as a fresh process would start.
       comps(ATM) = component_t('atmosphere', dt_atm, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
       comps(OCE) = component_t('ocean',      dt_oce, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
       comps(ICE) = component_t('ice_sheet',  dt_ice, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
       accs = accumulator_t()
       ice_lagged = 0.0_dp
       call read_restart('restart.bin', comps, accs)
       ! TODO 5: ice_lagged is state too. Is it in your restart file?
       i0 = stop_at + 1
    end if

    call step_range(i0, n, ice_lagged)
    final_state = comps(ICE)%total_received
    if (t < -1.0_dp .or. f < -1.0_dp) continue
  end function run_to

  subroutine step_range(ifrom, ito, ice_lagged)
    integer,  intent(in)    :: ifrom, ito
    real(dp), intent(inout) :: ice_lagged
    integer  :: i
    real(dp) :: t, f, avg

    do i = ifrom, ito
       t = real(i-1, dp) * dt_atm
       f = surface_flux(t)

       call accumulate(accs(1), f, dt_atm)
       call accumulate(accs(2), f, dt_atm)

       if (mod(i * int(dt_atm), int(dt_oce)) == 0) then
          avg = window_average(accs(1))
          comps(OCE)%received = avg
          comps(OCE)%total_received = comps(OCE)%total_received + avg * dt_oce
       end if

       if (mod(i * int(dt_atm), int(dt_ice)) == 0) then
          ! Close the window and PARK the average for the next one.
          ice_lagged = window_average(accs(2))
       end if

       ! TODO 5: the ice sheet consumes the LAGGED average, once per its own
       ! timestep. Add ice_lagged * dt_ice to comps(ICE)%total_received at
       ! each ice timestep.
       if (mod(i * int(dt_atm), int(dt_ice)) == 0) then
          ! comps(ICE)%total_received = comps(ICE)%total_received + ice_lagged * dt_ice
          continue
       end if
    end do
  end subroutine step_range

  subroutine check(label, err, tol)
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: err, tol
    logical :: ok
    ok = (err <= tol)
    ntest = ntest + 1
    if (ok) npass = npass + 1
    print '(a,a28,a,es12.4,a)', '  ', label, '  err ', err, &
         merge('   PASS', '   FAIL', ok)
  end subroutine check

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) nsteps
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) restart_at
    end if
  end subroutine read_cli

end program multirate

! ===========================================================================
! TODO 6 — write your answers in notes/day3.md
!
! (a) How wrong was the naive instantaneous sampling? Now shift the melt
!     event by 50 time units and re-run. Does the error change sign? Explain
!     why an error that changes sign with phase is MORE dangerous in a long
!     climate run than a consistent bias would be.
!
! (b) The accumulator holds sum(flux*dt), not mean(flux). Construct a case
!     with variable timesteps where the two give different answers, and give
!     the relative error. (Real models shorten the timestep during violent
!     weather — precisely when the flux is largest.)
!
! (c) You restarted at step 500, which is not a multiple of the ice window
!     (100). What is in the accumulator at that moment, and what happens if
!     you leave it out of the restart file? Quantify the error. Now restart
!     at step 400 instead — does the bug disappear? Why is a bug that only
!     appears at some restart points worse than one that always appears?
!
! (d) The ice sheet lags by one coupling window (100 time units). At a
!     realistic ice-sheet coupling interval of one year, what does that lag
!     mean physically, and when does it stop being acceptable? (Consider a
!     rapidly retreating marine-terminating glacier.)
!
! (e) Now the performance question, which is the one the posting is really
!     asking. The ice sheet is 100x slower per call but called 100x less
!     often. Sketch the coupled timeline. Is the ice sheet on the critical
!     path? What if it takes longer than 100 atmosphere steps to run?
!     Name two ways to keep it off the critical path. (You implement one of
!     them in Exercise 23.)
!
! (f) Conservation across a coupling interface is checked in real models by
!     a "conservation diagnostic" that reports the imbalance every N steps.
!     Design one for this setup: what do you sum on each side, and what
!     tolerance would you flag on? Why can the tolerance not be zero?
! ===========================================================================
