!! This module handles reading the namelist and provides access to some other flags
!! that control a specific CARMA model's behavior.
!!
!! By default the specific CARMA model does not have any unique namelist values. If
!! a CARMA model wishes to have its own namelist, then this file needs to be copied
!! from physics/cam to physics/model/<model_name> and the code needed to read in the
!! namelist values added there. This file will take the place of the one in
!! physics/cam. 
!!
!! It needs to be in its own file to resolve some circular dependencies.
!!
!! @author  Chuck Bardeen
!! @version Mar-2011
module carma_model_flags_mod

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc

  ! Flags for integration with CAM Microphysics
  public carma_model_readnl                   ! read the carma model namelist
  

  ! Namelist flags
  !
  ! Create a public definition of any new namelist variables that you wish to have,
  ! and default them to an inital value.
  character(len=256), public     :: carma_mice_file      = 'mice_warren2008.nc'   ! name of the ice refractive index file
  character(len=32), public      :: carma_sulfate_method = "fixed"                ! prescribed sulfate method
  logical, public                :: carma_do_initice     = .false.                ! If .true. than initialize carma ice bins from CLDICE
  logical, public                :: carma_do_initliq     = .false.                ! If .true. than initialize carma liquid bins from CLDLIQ
  logical, public                :: carma_do_mass_check  = .false.                ! If .true. then CARMA will check for mass loss by CARMA
  logical, public                :: carma_do_mass_check2 = .false.                ! If .true. then CARMA will check for mass loss (internal steps, e.g. detrain, diagnoseBIns, ...)
  logical, public                :: carma_do_mass_check3 = .false.                ! If .true. then CARMA will check for incoming mass loss (CAM -> CARMA)
  logical, public                :: carma_do_mass_fix    = .false.                ! If .true. then CARMA will fix for mass loss between cldice and ice bins
  logical, public                :: carma_do_print_fix   = .false.                ! If .true. then CARMA will print the value of the mass fix  
  integer, public                :: carma_dropact_bin    = 1                      ! Indicates the bin number into which activated droplets will be placed.

contains


  !! Read the CARMA model runtime options from the namelist
  !!
  !! @author  Chuck Bardeen
  !! @version Mar-2011
  subroutine carma_model_readnl(nlfile)
  
    ! Read carma namelist group.
  
    use cam_abortutils,  only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
  
    ! args
  
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  
    ! local vars
  
    integer :: unitn, ierr
  
    ! read namelist for CARMA
    namelist /carma_model_nl/ &
      carma_mice_file, carma_sulfate_method, carma_do_initice, carma_do_initliq, &
      carma_do_mass_check, carma_do_mass_check2, carma_do_mass_check3, &
      carma_do_mass_fix, carma_do_print_fix, carma_dropact_bin
  
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'carma_model_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, carma_model_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('carma_model_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if
  
#ifdef SPMD
    call mpibcast(carma_mice_file,      len(carma_mice_file),       mpichar, 0, mpicom)
    call mpibcast(carma_sulfate_method, len(carma_sulfate_method),  mpichar, 0, mpicom)
    call mpibcast(carma_do_initice,     1,  mpilog, 0, mpicom)
    call mpibcast(carma_do_initliq,     1,  mpilog, 0, mpicom)
    call mpibcast(carma_do_mass_check,  1,  mpilog, 0, mpicom)
    call mpibcast(carma_do_mass_check2, 1,  mpilog, 0, mpicom)
    call mpibcast(carma_do_mass_check3, 1,  mpilog, 0, mpicom)
    call mpibcast(carma_do_mass_fix,    1,  mpilog, 0, mpicom)
    call mpibcast(carma_do_print_fix,   1,  mpilog, 0, mpicom)
    call mpibcast(carma_dropact_bin,    1,  mpiint, 0, mpicom)
#endif
  
  end subroutine carma_model_readnl

end module carma_model_flags_mod
