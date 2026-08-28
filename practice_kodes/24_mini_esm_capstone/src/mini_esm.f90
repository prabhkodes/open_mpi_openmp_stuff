! ===========================================================================
! EXERCISE 24 — CAPSTONE B: the driver
!
! Read SPEC.md first. Build in the order it gives, and VERIFY EACH STEP
! before moving on -- a conservation bug introduced at step 5 and found at
! step 8 will cost you the whole session.
! ===========================================================================
program mini_esm
  use kinds
  use mpi_f08
  use mesh_mod
  use field_mod
  use kernels_mod
  use remap_mod
  use couple_mod
  implicit none

  type(comp_t)  :: c
  type(mesh_t)  :: mesh
  type(field_t) :: u, g
  type(accum_t) :: acc
  type(remap_t) :: R

  integer :: ierr, n_atm = -1
  integer :: nsteps = 100, nlev = 40, nproma = 64
  integer :: ncells_atm = 40000, ncells_oce = 25000
  integer :: ncells, step
  real(wp) :: cons_error

  call MPI_Init(ierr)
  call read_cli()

  block
    integer :: w
    call MPI_Comm_size(MPI_COMM_WORLD, w, ierr)
    if (n_atm < 0) n_atm = max(1, w / 2)
    if (n_atm >= w) n_atm = max(1, w - 1)
  end block

  call comp_init(c, n_atm)
  ncells = merge(ncells_atm, ncells_oce, c%id == COMP_ATM)

  if (c%world_rank == 0) then
     print '(a)', '=== Exercise 24 (CAPSTONE B): mini coupled ESM ==='
     print '(a,i0,a,i0,a,i0)', 'ranks ', c%world_size, ' = atm ', n_atm, &
          ' + oce ', c%world_size - n_atm
     print '(a,i0,a,i0)', 'cells: atm ', ncells_atm, '   oce ', ncells_oce
     print '(a,i0,a,i0)', 'nlev ', nlev, '   nproma ', nproma
     print '(a)', ''
     print '(a)', 'Work through SPEC.md in order. Each stage has a'
     print '(a)', 'verification -- do not skip past a failing one.'
     print '(a)', ''
  end if

  ! ---- stage 1-3: mesh, partition, halo (SPEC section 1) -----------------
  call mesh_build(mesh, ncells)
  call mesh_partition(mesh, c%comm)
  call mesh_halo_schedule(mesh, c%comm)
  if (c%world_rank == 0) print '(a)', '  --- mesh ---'
  call mesh_stats(mesh, c%comm, trim(c%name))

  ! ---- stage 2: fields (SPEC section 2) ----------------------------------
  call field_alloc(u, ncells, nlev, nproma)
  call field_alloc(g, ncells, nlev, nproma)

  ! ---- stage 5: remap weights, built ONCE (SPEC section 3) ---------------
  ! TODO: build the atm<->oce remap here, at setup, not in the timeloop.

  ! ---- stage 6: the coupled timeloop -------------------------------------
  do step = 1, nsteps
     ! TODO: halo exchange, horizontal kernel, column kernel,
     !       accumulate the coupling flux, exchange at the coupling period,
     !       and CHECK CONSERVATION every step.
     continue
  end do

  ! ---- report ------------------------------------------------------------
  cons_error = 0.0_wp
  if (c%world_rank == 0) then
     print '(a)', ''
     print '(a)', '  --- report (fill the SPEC.md table from these) ---'
     print '(a,es12.4)', '    conservation error : ', cons_error
     print '(a)',        '    (0.0 means the coupling is not implemented yet)'
     print '(a)', ''
     print '(a)', '  Next: SPEC.md section 5 (scaling) and 6 (tests).'
  end if

  call field_free(u); call field_free(g)
  call MPI_Finalize(ierr)

contains

  subroutine read_cli()
    character(len=32) :: arg
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg); read(arg,*) nsteps
    end if
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg); read(arg,*) nproma
    end if
    if (command_argument_count() >= 3) then
       call get_command_argument(3, arg); read(arg,*) n_atm
    end if
  end subroutine read_cli

end program mini_esm
