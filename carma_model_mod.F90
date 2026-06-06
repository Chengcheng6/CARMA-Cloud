!!  This module is used to define a particular CARMA microphysical model. For
!! simple cases, this may be the only code that needs to be modified.
!!                                          - Chuck Bardeen circa 2010
!!
!!  This module is for the Carma Cloud model. This module is based on the
!! Carma Cirrus model written by Chuck Bardeen (CB).
!!                                          - Jamison A Smith (JAS) 2019
!! ----------------------------------------------------------------------------
!! @version 2020 October
!! @author Jamison A. Smith (JAS)
!! @version 2019 November
!! @author JAS
!!
!!  Prescribed sulfuric acid aerosol particles are used for nucleation of in
!! situ cirrus as written by Chuck Bardeen.  This prescribed sulfuric acid
!! aerosol is a diagnostic aerosol.  The aerosol is diagnosed from a seasonal
!! cycle of sulfuric acid aerosol produced in another simulation where the
!! sulfuric acid is a prognostic aerosol that is produced by the nucleation
!! of sulfuric acid molecules.  Jason English produced these simulations.
!!  Mam aerosol are activated to produce water droplets. This activation is
!! done with Subroutine dropmixnuc ( ndrop.F90 ) in the Cam. Activation is
!! occuring on aerosol particles that are outside of the Carma.
!!  Later Mam aerosol may be converted into a new diagnostic group to use the
!! Carma's activation code.
!!  Use small indicies for aerosol groups. Use large indicies for cloud
!! groups. The numerical methods act differently on particle bins depending on
!! whether the particles are moving upwards in the particle grid ( Subroutine
!! upgxfer, upgxfer.F90 ) or downwards ( Subroutine downgxfer, downgxfer.F90 ).
!! ----------------------------------------------------------------------------
!! + CB circa 2010
!!
!! This module is used to define a particular CARMA microphysical model. For 
!! simple cases, this may be the only code that needs to be modified. This module
!! defines several constants and has the following methods:
!!
!!   - CARMA_DiagnoseBins()
!!   - CARMA_DiagnoseBulk()
!!   - CARMA_DefineModel()
!!   - CARMA_Detrain()
!!   - CARMA_EmitParticle()
!!   - CARMA_InitializeModel()
!!   - CARMA_InitializeParticle()
!!
!! These methods define the microphysical model, the particle emissions and
!! the initial conditions of the particles. For diagnostic groups, there are
!! also routines that diagnose the mass in the bins of that group from the
!! parent model's state inforamtion and that calculate the tendency on the
!! parent model's state based upon changes in the bins.
!!
!! This cirrus cloud model allows CARMA bin microphysics to do the ice microphysics
!! while MG does the liquid microphysics. The MG microphysics here should not update
!! CLDICE or NUMICE, since those values will not be reflected in the CARMA ice
!! bins, which are the true state variables for ice. In this situation, CLDICE and
!! NUMICE are merely diagnostic variables available as input to the rest of CAM.
!!
!! The CARMA microphysics will run before MG and will handle:
!!   - Detrainment (liquid and ice)
!!   - Homogeneous ice nucleation (currently with prescribed sulfates)
!!   - Heterogeneous ice nucleation (future)
!!   - Bergeron process
!!   - Melting of detrained ice
!!   - Freezing of cloud drops
!!   - Autoconversion (ice -> snow)
!!   - Variable ice density (function of particle size)
!!   - In-cloud values (dividing by cloud fraction)
!!
!! Some potential issues that are not currently handled by CARMA:
!!   - collection of ice by snow
!!   - aggregation of ice
!!   - sub-grid vertical velocity for CARMA
!!   - Goff & Gratch vs. Murphy & Koop vapor pressures
!!   - Radiation using CARMA size distribution (each bin as tracer)
!!   - Hallet-Mossop Process
!!
!! The following variables will have been set by CARMA:
!!   - (S) CLDICE, (S) NUMICE
!!   - (S) CLDLIQ, (S) NUMLIQ
!!   - (S) T
!!   - (P) TNDQSNOW,  (P) TNDNSNOW
!!   - (P) REICE
!!
!! Variables with an S will be in the physics_state and variables with a P are
!! parameters passed into the MG microphysics.
!!
!! The module carma_intr defines a few flags that indicate what portion of the
!! cloud microphysics is handled by CARMA:
!! ---
!! JAS:  I apologize for interrupting, but carma_do_cldice and carma_do_cldliq
!!      are misnomers. See my discussion below.
!! --- 
!!   - carma_do_cldice  - CARMA does ice clouds
!!   - carma_do_cldliq  - CARMA does liquid clouds
!!
!! - CB circa 2010
!! ----------------------------------------------------------------------------
!!
!! + JAS:  Discussion of carma_do_cldliq and carma_do_cldice being misnomers
!!
!!  carma_do_cldliq is a misnomer. Carma is going to calculate changes in
!! cloud liquid when carma_do_cldliq is true AND when carma_do_cldliq is false.
!! In the Carma Cirrus model, Carma calculates changes in cloud liquid even
!! though carma_do_cldliq is false. The changes are from the melting of Carma
!! Ice and the freezing of Carma Liquid. Liquid-Ice collisions are possible,
!! too.
!!  carma_do_cldliq is really answering the question, "Am I accepting the Carma
!! Liquid Bins as prognostic variables from Dynamics (Finite Volume) and
!! SubGridTurb (Clubb)?" The alternative, when carma_do_cldliq = false, is to
!! diagnose the Carma liquid bins from CLDLIQ, NUMLIQ, and the gamma function.
!!
!!  In the Carma Cirrus model, carma_do_cldliq = false
!! 1) Carma diagnoses CRLIQxx Bins (xx = 01-28) from CLDLIQ, NUMLIQ, and the
!!    gamma function
!! 2) Carma allows the Carma Ice Bins (CRDICExx, CRSICExx, CRCORxx) and Liquid
!!    Bins (CRLIQxx) to exchange material via physical changes
!! 3) Carma diagnoses the new Bulk variables for cloud liquid CLDLIQ and
!!    NUMLIQ from the Carma Bins
!! 4) Carma calculates the tendencies for CLDLIQ and NUMLIQ
!!
!! - JAS:  End of Discussion of carma_do_cldice and carma_do_cldliq...
!!
!! ----------------------------------------------------------------------------
!!
!! Each realization of CARMA microphysics has its own version of this file.
!!
!! This model replaces the ice microphysics from the MG two-moment scheme with
!! a CARMA bin microphysics representation of the ice. The purpose of this
!! model is to provide a more detail description of the thin cirrus clouds that
!! form in the TTL and to investigate the impact of these clouds on radiative
!! forcing, troposphere-to-stratosphere transport, and control of water vapor
!! in the UT/LS.
!!
!! @version July-2009 
!! @author  Chuck Bardeen 
!!
!! @version 2020 February - November 2019
!! @author  Jamison A. Smith (JAS) 

module carma_model_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagas_mod
  use carmagroup_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod 

  ! ---------------------------------------------------------------------------

  use carma_flags_mod  ! JAS: The value of carma_do_cldice comes from here

  !  JAS:  physics/cam/carma_flags_mod.F90 Sub read_carma_nl reads
  ! carma_do_cldice from the Carma namelist

  ! ---------------------------------------------------------------------------

  use carma_model_flags_mod
  
  use spmd_utils,     only: masterproc
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use radconstants,   only: nswbands, nlwbands
  use cam_abortutils, only: endrun
  use physics_types,  only: physics_state, physics_ptend
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_field, pbuf_get_index
  use physconst,      only: gravit

  ! JAS:  Need consistent definition of QSmall to avoid Floating Point
  !      Exceptions (FPEs)

  Use Micro_MG_Utils,  Only: QSmall 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! LW (Sep 2022): added for homo. & hetero. nuc.
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_num, rad_cnst_get_aer_mmr, &
                            rad_cnst_get_aer_props, rad_cnst_get_mode_props, &
                            rad_cnst_get_mode_num_idx, rad_cnst_get_mam_mmr_idx
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#if ( defined SPMD )
  use mpishorthand
#endif  

  implicit none

  private

  ! Declare the public methods.
  public CARMA_DefineModel
  public CARMA_Detrain
  public CARMA_DiagnoseBins
  public CARMA_DiagnoseBulk
  public CARMA_EmitParticle
  public CARMA_InitializeModel
  public CARMA_InitializeParticle
  public CARMA_OutputBudgetDiagnostics
  public CARMA_OutputDiagnostics
  public CARMA_WetDeposition
  
  ! Declare public constants

  integer, public, parameter      :: NGROUP   = 6               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 8               !! Number of particle elements
! JAS 2020 Nov 20 integer, public, parameter      :: NBIN     = 28              !! Number of particle bins
                  integer, public, parameter      :: NBIN     = 48              !! Number of particle bins  ! JAS
! JAS 2021 Jun 29 integer, public, parameter      :: NBIN     = 60              !! Number of particle bins  ! JAS
  integer, public, parameter      :: NSOLUTE  = 1               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 1               !! Number of gases

!                                                                                                  1
!        1         2         3         4         5         6         7         8         9         0
!234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890

  !  JAS 2020 November 2: Chuck's text is too wide for me to fit two windows side by side on my
  ! larger screen.  Reduce the width.  It doesn't have to be 80 chars, but 128 is too wide for me to
  ! fit two windows on the same screen.

  ! These need to be defined, but are only used when the particles are radiatively active.

  !! Number of relative humidities for mie calculations

  integer, public, parameter      :: NMIE_RH  = 8
  real(kind=f), public            :: mie_rh(NMIE_RH)
  
  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.

  !! Should the group be transported in the first phase?

  logical, public                 :: is_convtran1(NGROUP) = .false.

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.

  integer, public, parameter      :: I_H2SO4   = 1               !! sulfate aerosol composition
  integer, public, parameter      :: I_WATER   = 2               !! water
  integer, public, parameter      :: I_ICE     = 3               !! ice
  integer, public, parameter      :: I_DUST    = 4               !! coarse mode MAM4 dust (LW: added for hetero. nuc.)
 
  ! Define group, element, solute and gas indexes.

  integer, public, parameter      :: I_GRP_CRCN     = 1             !! sulfate aerosol
  integer, public, parameter      :: I_GRP_CRLIQ    = 2             !! liquid drop
  integer, public, parameter      :: I_GRP_CRDICE   = 3             !! detrained ice
  integer, public, parameter      :: I_GRP_CRSICE   = 4             !! in-situ ice
  integer, public, parameter      :: I_GRP_CRGRP    = 5             !! graupel
  integer, public, parameter      :: I_GRP_CRDUST   = 6             !! dust aerosol (LW: added for hetero. nuc.)

  integer, public, parameter      :: I_ELEM_CRCN    = 1             !! sulfate
  integer, public, parameter      :: I_ELEM_CRLIQ   = 2             !! liquid water
  integer, public, parameter      :: I_ELEM_CRDICE  = 3             !! detrained ice
  integer, public, parameter      :: I_ELEM_CRSICE  = 4             !! in-situ ice
  integer, public, parameter      :: I_ELEM_CRCORE  = 5             !! sulfate core
  integer, public, parameter      :: I_ELEM_CRDCOR  = 6             !! dust core (LW: added for hetero. nuc.; CM:combination of all modes)
  integer, public, parameter      :: I_ELEM_CRGRP   = 7             !! graupel core
  integer, public, parameter      :: I_ELEM_CRDUST  = 8             !! dust (LW: added for hetero. nuc.)

  integer, public, parameter      :: I_SOL_CRH2SO4  = 1             !! sulfuric acid

  integer, public, parameter      :: I_GAS_H2O      = 1             !! water vapor


  ! From Morrison & Gettelman [2008] and micro_mg.F90 (formerly cldwat2m_micro.F90)
  !
  ! NOTE: In the bin model, the bin boundaries are also important for determining the threshold,
  ! since the whole bin is autoconverted if the threshold is less than the bin midpoint radius.

  real(kind=f), public, parameter :: CAM_RHOCI = 0.5_f    !! (g/cm3) MG bulk density for cloud ice
  real(kind=f), public, parameter :: CAM_RHOSN = 0.1_f    !! (g/cm3) MG bulk density for snow

  
  ! Parameters and variabls that control the detrainment process.
  integer, parameter :: NINTS_BINS = 10        !! number of steps to integrate bin fractions
    
  real(kind=f), parameter :: r_dliq_lnd  =  8e-4_f   !! detrained liquid radius (cm)
! real(kind=f), parameter :: r_dliq_lnd  = 18e-4_f   !! detrained liquid radius (cm)
  real(kind=f), parameter :: r_dliq_ocn  =  8e-4_f  !! detrained liquid radius (cm)
! real(kind=f), parameter :: r_dliq_ocn  = 14e-4_f  !! detrained liquid radius (cm)
! real(kind=f), parameter :: r_dliq_ocn  = 18e-4_f  !! detrained liquid radius (cm)

!  integer, parameter :: MIN_DTEMP      = -60            !! Miniumum detrainment temperature (C)
  integer, parameter  :: MIN_DTEMP      = -90            !! Miniumum detrainment temperature (C)
  integer, parameter  :: NDTEMP         = -MIN_DTEMP + 1 !! Number of detrainment temperature bins
  
! character(len=12), parameter    :: carma_dice_method    = "mono"
  real(kind=f), parameter :: r_dice_mono          = 25e-4_f   !! detrained ice radius, monodisperse (cm)

  ! This distribution varies the size disribution as a function of temperature, with the
  ! distribution biased toards larer particles at warm temperature and small particles at
  ! cold temperatures. This fit is from eq. 7 of Heymsfield and Schmitt [2010]. The Jensen
  ! fit used above is similar to the cold end of this range.
  character(len=12), parameter    :: carma_dice_method    = "dist_hym2010"  ! From eq 7 in Heymsefield & Schmitt [2010] (cm -1)
!  real(kind=f), parameter         :: dist_hym2010_alpha   = 14.26_f       !! alpha (stratiform)
!  real(kind=f), parameter         :: dist_hym2010_beta    = -0.0538_f     !! beta (stratiform)
  real(kind=f), parameter         :: dist_hym2010_alpha   = 2.425_f        !! alpha (convective)
  real(kind=f), parameter         :: dist_hym2010_beta    = -0.088_f       !! beta (convective)
  real(kind=f)                    :: dice_bin_fraction(NBIN, NDTEMP)       !! detrained mass fraction, ice bin

  logical, public, parameter      :: carma_do_bulk_tend  = .true.          ! If .true. then update CAM bulk tendencies

  integer                         :: ixcldice
  integer                         :: ixnumice
  integer                         :: ixcldliq
  integer                         :: ixnumliq

  ! JAS: MG 2.0 has prognostic precipitation variables

  integer                         :: ixrainqm
  integer                         :: ixnumrai
  integer                         :: ixsnowqm
  integer                         :: ixnumsno

  integer                             :: warren_nwave            ! number of wavelengths in file
  real(r8), allocatable, dimension(:) :: warren_wave             ! Warren & Brandt 2008, wavelengths
  real(r8), allocatable, dimension(:) :: warren_real             ! Warren & Brandt 2008, real part of m
  real(r8), allocatable, dimension(:) :: warren_imag             ! Warren & Brandt 2008, imag part of m

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! LW (Sep 2022): This section is added for homo. & hetero. nuc., 
!                based on C. Maloney's version for CESM1


  ! MAM 4 aerosol variables
  integer                         :: ixso4_a1
  integer                         :: ixso4_a2
  integer                         :: ixso4_a3
  integer                         :: ixnum_a1
  integer                         :: ixnum_a3
  integer                         :: ixdust_a1
  integer                         :: ixdust_a3


  integer :: ntote_amode              ! number of aerosol modes
  integer :: ncnst_tot                ! total number of mode number concentration and mode species
    type(physics_buffer_desc), pointer   :: pbuf(:)   !! physics buffer

  ! ptr2d_t is used to create arrays of pointers to 2D fields
  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t  

  ! mode aerosol indexing variables
   integer :: mode_accum_idx  = -1  ! index of sulfate accumulation mode
   integer :: mode_aitken_idx = -1  ! index of sulfate aitken mode
   integer :: mode_coarse_idx = -1  ! index of sulfate coarse mode
   integer :: so4_coarse_idx  = -1  ! index of sulfate coarse mode
   integer :: so4_accum_idx   = -1  ! index of sulfate accumulation mode
   integer :: so4_aitken_idx  = -1  ! index of sulfate aitken mode
   integer :: sslt_coarse_idx = -1  ! index of sea salt coarse mode
   integer :: sslt_accum_idx  = -1  ! index of sea salt accumulation mode
   integer :: dust_coarse_idx = -1  ! index of dust coarse mode
   integer :: dust_accum_idx  = -1  ! index of dust accumulation mode
   integer :: bc_accum_idx    = -1  ! index of black carbon accumulation mode
   integer :: pom_accum_idx   = -1  ! index of primary organics accumulaiton mode
   integer :: soa_accum_idx   = -1  ! index of secondary organics accumulation mode

  ! mode geometrical sigma variables
   real(r8)              :: sigmag_coarse                ! coarse
   real(r8)              :: sigmag_accum                 ! accumulation
   real(r8)              :: sigmag_aitken
   real(r8)              :: dgnumlo_co
   real(r8)              :: dgnumhi_co
   real(r8)              :: dgnumlo_at
   real(r8)              :: dgnumhi_at
   real(r8)              :: dgnumlo_ac
   real(r8)              :: dgnumhi_ac

    ! Indexing variables for modes and species (3 species so4, dust, sea salt)
    integer :: m, nn, ii, l

    ! Aerosol constituents info
    character(len=32)     :: str32
    integer               :: nspec                       ! number of chemical species
    integer               :: nmodes                       ! number of aerosol modes
!    integer               :: dgnum_idx   = -1
!    real(r8), pointer     :: dgnum(:,:,:)                        ! aerosol mode dry diameter (m)
    real(r8), pointer                     :: so4_mmr_co(:,:)    ! interstitial sulfate mmr, coarse mose (kg/kg)
    real(r8), pointer                     :: so4_mmr_ac(:,:)    ! interstitial sulfate mmr, accumulation mode (kg/kg)
    real(r8), pointer                     :: sslt_mmr_co(:,:)    ! interstitial sea salt mmr, coarse mode (kg/kg)
    real(r8), pointer                     :: sslt_mmr_ac(:,:)
    real(r8), pointer                     :: dust_mmr_co(:,:)   ! interstitial dust mmr, coarse mode (kg/kg)
    real(r8), pointer                     :: dust_mmr_ac(:,:)   ! interstitial dust mmr, accumulation mode (kg/kg)
    real(r8), pointer                     :: num_coarse(:,:)    ! number mixing ratio of coarse mode (#/kg)
    real(r8), pointer                     :: num_accum(:,:)
   ! For DiagnoseBulk
    real(r8), pointer                    :: so4_mmr(:,:)             ! interstitial sulfate mmr
    real(r8), pointer                    :: dust_mmr(:,:)             ! interstitial dust mmr
    real(r8), pointer                    :: sslt_mmr(:,:)             ! interstitial sea salt mmr
    real(r8), pointer                    :: mode_nmr(:,:)             ! interstitial nmr (#/kg)
    real(r8), pointer                    :: mode_nmr_co(:,:)             ! interstitial nmr (#/kg)
    real(r8), pointer                    :: mode_nmr_ac(:,:)             ! interstitial nmr (#/kg)
    real(r8), pointer                    :: bc_mmr_ac(:,:)               ! interstitial black carbon mmr (kg/kg)               
    real(r8), pointer                    :: soa_mmr_ac(:,:)              ! itnerestitial secondary organics mmr
    real(r8), pointer                    :: pom_mmr_ac(:,:)             ! interstitial primary organcis mmr 

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


contains

  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  !!  @version December 2019
  !!  @author  Jamison A. Smith (JAS)
 
  subroutine CARMA_DefineModel(carma, rc)
    use physconst,          only: latice, latvap
    use ioFileMod,          only: getfil
    use wrap_nf

    implicit none
    
    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure
        
    ! Local variables

    !real(kind=f), parameter            :: rmin_ice   = 5.e-5_f  ! min radius for ice bins (cm)
    real(kind=f), parameter            :: rmin_ice   = 3.e-5_f  ! min radius for ice bins (cm)        ! JAS 2020Nov20
    real(kind=f), parameter            :: rmin_cn    = 1.e-7_f  ! min radius for sulfate bins (cm)
    !real(kind=f), parameter            :: rmin_cn    = 5.e-8_f  ! min radius for sulfate bins (cm)    ! JAS 2020Nov20
    real(kind=f), parameter            :: RHO_CN     = 1.78_f   ! density of sulfate particles (g/cm)
    real(kind=f), parameter            :: RHO_DUST   = 2.65_f    ! dry density of dust particles (g/cm^3) -Lin Su (LW: added for hetero. nuc.)
    real(kind=f), parameter            :: rmin_dust  = 1.19e-5_f ! minimum radius (cm) (LW: added for hetero. nuc.)
    real(kind=f), parameter            :: vmrat_dust = 1.305_f   ! volume ratio (LW: 48-bin dust, added for hetero. nuc.)
    real(kind=f)                       :: rmassmin              ! mass of the first radius bin (g)
    real(kind=f)                       :: vmrat                 ! volume ratio between adjacent bin
    real(kind=f)                       :: rhoelem(NBIN)         ! element density per bin (g/cm3)
    real(kind=f)                       :: arat(NBIN)            ! projected area ratio
    integer                            :: i
    integer                            :: j
    real(kind=f)                       :: wave(NWAVE)               ! CAM band wavelength centers (cm)
    integer                            :: fid
    integer                            :: wave_did
    integer                            :: wave_vid
    integer                            :: real_vid
    integer                            :: imag_vid
    character(len=256)                 :: efile                     ! refractive index file name
    real(kind=f)                       :: interp
    complex(kind=f)                    :: refidx_ice(NWAVE)         ! the refractive index at each CAM wavelength    
    integer                            :: LUNOPRT
    logical                            :: do_print
    complex(kind=f)                    :: refidx_liq(NWAVE)         ! the refractive index for liquid water at each CAM wavelength [added by LW (Mar 2022)]    
    real(kind=f)                       :: liq_real(NWAVE)=(/2.200349_f,1.529309_f,1.398169_f,1.24938_f,1.135773_f,1.116059_f,1.206729_f,1.242789_f,&
                                                            1.271952_f,1.295193_f,1.330121_f,1.29755_f,1.311588_f,1.319948_f,1.333929_f,1.39326_f,&
                                                            1.39326_f,1.127959_f,1.275295_f,1.294919_f,1.304155_f,1.312888_f,1.316645_f,1.321521_f,&
                                                            1.328769_f,1.336943_f,1.350438_f,1.369839_f,1.410702_f,1.268163_f/)         ! real part of the liq water refractive index [added by LW (Mar 2022)] 
    real(kind=f)                       :: liq_imag(NWAVE)=(/7.13e-01_f,3.48e-01_f,4.26e-01_f,4.04e-01_f,3.26e-01_f,1.13e-01_f,4.72e-02_f,3.88e-02_f,&
                                                            3.37e-02_f,3.20e-02_f,1.17e-01_f,1.01e-02_f,1.55e-02_f,8.45e-03_f,4.60e-03_f,1.32e-02_f,&
                                                            1.32e-02_f,1.02e-01_f,5.43e-04_f,7.39e-04_f,1.22e-04_f,3.60e-04_f,1.09e-05_f,2.69e-06_f,&
                                                            3.35e-08_f,1.89e-09_f,1.66e-09_f,3.83e-09_f,1.10e-08_f,3.43e-02_f/)         ! imaginary part part of the liq water refractive index [added by LW (Mar 2022)] 
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! LW (Sep 2022): added for homo. & hetero. nuc.
!    type(physics_buffer_desc), pointer    :: pbuf(:)      !! physics buffer

    character(len=*), parameter :: routine = 'CARMA_DefineModel'
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Default return code.

    rc = RC_OK
    
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT, wave=wave)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_Get failed.')
    
    ! Report model specific configuration parameters.

    if (masterproc) then
      if (do_print) then
        write(LUNOPRT,*) ''
        write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
        write(LUNOPRT,*) '  carma_mice_file       = ', trim(carma_mice_file)
        write(LUNOPRT,*) '  carma_sulfate_method  = ', trim(carma_sulfate_method)
      end if
    end if

    ! Get the refractive index for ice as a function of wavelength for particle heating
    ! calculations.
    !
    ! NOTE: These values probably should be a band average, but for now just do band centers.
    
    ! Read the values in from Warren et al. 2008.

    !if (carma_do_pheat) then
      if (masterproc) then 
      
        ! Open the netcdf file (read only)

        call getfil(carma_mice_file, efile, fid)
        if (do_print) write(LUNOPRT,*) 'carma_init(): Reading ice refractive indexes from ', efile
  
        call wrap_open(efile, 0, fid)
  
        ! Alocate the table arrays

        call wrap_inq_dimid(fid, "wavelength", wave_did)
        call wrap_inq_dimlen(fid, wave_did, warren_nwave)
      endif
      
#if ( defined SPMD )
        call mpibcast(warren_nwave, 1, mpiint, 0, mpicom)
#endif
  
        allocate(warren_wave(warren_nwave))
        allocate(warren_real(warren_nwave))
        allocate(warren_imag(warren_nwave))
  
        if (masterproc) then
          
          ! Read in the tables.

          call wrap_inq_varid(fid, 'wavelength', wave_vid)
          call wrap_get_var_realx(fid, wave_vid, warren_wave)
          warren_wave = warren_wave * 1e-4          ! um -> cm
  
          call wrap_inq_varid(fid, 'm_real', real_vid)
          call wrap_get_var_realx(fid, real_vid, warren_real)
  
          call wrap_inq_varid(fid, 'm_imag', imag_vid)
          call wrap_get_var_realx(fid, imag_vid, warren_imag)
  
          ! Close the file.

          call wrap_close(fid)
        end if
  
#if ( defined SPMD )
        call mpibcast(warren_wave,  warren_nwave, mpir8, 0, mpicom)
        call mpibcast(warren_real,  warren_nwave, mpir8, 0, mpicom)
        call mpibcast(warren_imag,  warren_nwave, mpir8, 0, mpicom)
#endif
      
      ! Interpolate the values.

      do i = 1, NWAVE
        do j = 1, warren_nwave
          if (wave(i) <= warren_wave(j)) then
            if ((j > 1) .and. (wave(i) /= warren_wave(j))) then
              interp = (wave(i) - warren_wave(j-1)) / (warren_wave(j) - warren_wave(j-1))
              refidx_ice(i) = cmplx(warren_real(j-1) + interp*(warren_real(j) - warren_real(j-1)), &
                   warren_imag(j-1) + interp*(warren_imag(j) - warren_imag(j-1)))
            else
              refidx_ice(i) = cmplx(warren_real(j), warren_imag(j))
            endif
            
            exit
          end if
        end do
      end do
    !end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Mar 2022)
    ! Get the refractive index for water
    refidx_liq = cmplx(liq_real, liq_imag)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !end if

    ! Define the Groups
    !
    ! NOTE: If NWAVE > 0 then the group should have refractive indices defined.
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.

    rmassmin = 4._f / 3._f * PI * rmin_cn ** 3 * RHO_CN

!    vmrat = 4.0_f     ! For 16 bins
!    vmrat = 2.8_f     ! For 21 bins
!    vmrat = 2.16_f    ! For 28 bins
!    vmrat = 2.0_f     ! For 32 bins
!    vmrat = 1.7_f     ! For 48 bins JAS 2020 Nov 20

    !  Be careful.  This is the definition for aerosol.  A second definition for hydrometeors
    ! appears below.

    !vmrat = 1.85_f     ! For 48 bins JAS 2021 Jun 29
    vmrat = 1.555_f     ! For 48 bins sulfate aersol, size range 0.1~1.01 µm (LW 07/2023)
    
    ! Since these sulfates are prescribed, don't sediment them. This will save some
    ! processing time.

    call CARMAGROUP_Create(carma, I_GRP_CRCN, "Sulfate CN", rmin_cn, vmrat, I_SPHERE, 1._f, .false., &
                           rc, shortname="CRCN", rmassmin=rmassmin, do_mie=.false., &
                           cnsttype=I_CNSTTYPE_DIAGNOSTIC, do_vtran=.false.)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')

    ! NOTE: For freezing and melting, the ice and water bins need to have the same mass.

    rmassmin = 4._f / 3._f * PI * rmin_ice ** 3 * RHO_I                 ! JAS 2020 Nov 20

    !vmrat = 2.055_f     ! For 28 bins, Heysmfield Ice Density, cold
    !vmrat = 1.7_f       ! For 48 bins, Heysmfield Ice Density, cold     ! JAS 2020 Nov 20
    vmrat = 1.85_f       ! For 48 bins, Rmax ~ 4.4 mm, JAS 2021 Jun 29

    ! Make the aged detrained ice have a variable density to represent the complex set of
    ! possible shapes that we can't represent. This is based upon Heymsfield and 
    ! Westfield [2010] and Heysfield and Schmitt [2010].

    call CARMAGROUP_Create(carma, I_GRP_CRDICE, "Detrained Ice, Aged", rmin_ice, vmrat, I_SPHERE, 1._f, .true., &
                           rc, shortname="CRDICE", rmassmin=rmassmin, do_mie=carma_do_pheat, refidx=refidx_ice, &
                           ifallrtn=I_FALLRTN_HEYMSFIELD2010, imiertn=I_MIERTN_BOHREN1983, is_cloud=.true.)

    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')

    is_convtran1 ( I_Grp_CrDIce ) = .true.

    ! Make the in-situ ice a plate, AR=6. This is based upon observations from Lawson
    ! et al. [2008]. AR=6 is for larger particles, so AR=3 is a compromise that is
    ! part way between that and more spheroidal particles that are likely at smaller sizes.
    !
    ! NOTE: All cloud particles should be convectively transported in the first phase of
    ! convection.
    !
    ! NOTE: All ice particles have the last bin as the one that gets autoconverted to
    ! snow at the end of the timestep and thus it does not need to be a prognostic bin.

!    call CARMAGROUP_Create(carma, I_GRP_CRSICE, "In-situ Ice", rmin_ice, vmrat, I_SPHERE, 1._f, .true., &
!    call CARMAGROUP_Create(carma, I_GRP_CRSICE, "In-situ Ice", rmin_ice, vmrat, I_HEXAGON, 1._f / 6._f, .true., &

    call CARMAGROUP_Create &
    (carma, I_GRP_CRSICE, "In-situ Ice", rmin_ice, vmrat, I_HEXAGON, 1._f / 3._f, .true., &
    rc, shortname="CRSICE", rmassmin=rmassmin, do_mie=carma_do_pheat, refidx=refidx_ice, &
    ifallrtn=I_FALLRTN_HEYMSFIELD2010, imiertn=I_MIERTN_BOHREN1983, &
    is_cloud=.true.)

    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
    is_convtran1 ( I_Grp_CrSIce ) = .true.

    ! Water drops are spherical.

    call CARMAGROUP_Create &
    (carma, I_GRP_CRLIQ, "Water Drop", rmin_ice, vmrat, I_Sphere, 1._f, .false. &
    , rc, shortname="CRLIQ", rmassmin = rmassmin, do_mie = .true., refidx=refidx_liq  &        ! JAS 2020 Oct
    , cnsttype=I_CnstType_Prognostic, is_cloud = .true., do_vtran = .true., imiertn=I_MIERTN_BOHREN1983)        ! JAS 2020 Oct

    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
    is_convtran1 ( I_Grp_CrLiq ) = .true.

    ! Graupel   
    call CARMAGROUP_Create &
    (carma, I_GRP_CRGRP, "Graupel", rmin_ice, vmrat, I_Sphere, 1._f, .true. &        !+++ ! is_ice = true
    , rc, shortname="CRGRP", rmassmin = rmassmin, do_mie = carma_do_pheat &        ! JAS 2020 Oct
    , cnsttype=I_CnstType_Prognostic, is_cloud = .true., do_vtran = .true.)        ! JAS 2020 Oct

    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
    is_convtran1 ( I_GRP_CRGRP ) = .true.

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Sep 2022)
    ! Dust - aerosol
    rmassmin = (4._f / 3._f) * PI * (rmin_dust ** 3) * RHO_DUST
    call CARMAGROUP_Create(carma, I_GRP_CRDUST, "Dust", rmin_dust, vmrat_dust, I_SPHERE, 1._f, .false., &
                           rc, shortname="CRDUST", rmassmin=rmassmin, do_mie=.false., &
                           cnsttype=I_CNSTTYPE_DIAGNOSTIC, do_vtran=.false.)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
    is_convtran1(I_GRP_CRDUST) = .false.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.

    call CARMAELEMENT_Create(carma, I_ELEM_CRCN, I_GRP_CRCN, "Sulfate CN", RHO_CN, &
         I_INVOLATILE, I_H2SO4, rc, shortname="CRCN", isolute=I_SOL_CRH2SO4)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    ! The density of ice is changed based on the maximum dimensions of ice particles
    ! as a function of mass from Heymsfield and Schmitt [2010].

!    call rhoice_heymsfield2010(carma, RHO_I, I_GRP_CRDICE, "conv", rhoelem, arat, rc)
    call rhoice_heymsfield2010(carma, RHO_I, I_GRP_CRDICE, "warm", rhoelem, arat, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::rhoice_heymsfield2010 failed.')
    
    call CARMAELEMENT_Create(carma, I_ELEM_CRDICE, I_GRP_CRDICE, "Detrained Ice", RHO_I, &
         I_VOLATILE, I_ICE, rc, shortname="CRDICE", rhobin=rhoelem, arat=arat)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')
    
    call CARMAELEMENT_Create(carma, I_ELEM_CRSICE, I_GRP_CRSICE, "In-situ Ice", RHO_I, &
         I_VOLATILE, I_ICE, rc, shortname="CRSICE")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! LW (Sep 2022): modified for hetero. nuc.
!    call CARMAELEMENT_Create(carma, I_ELEM_CRCORE, I_GRP_CRSICE, "Core Mass", RHO_CN, &
!         I_COREMASS, I_H2SO4, rc, shortname="CRCORE", isolute=I_SOL_CRH2SO4)
!    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CRCORE, I_GRP_CRSICE, "Sulfate Core Mass", & 
         RHO_CN, I_COREMASS, I_H2SO4, rc, shortname="CRCORE", isolute=I_SOL_CRH2SO4)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CRDCOR, I_GRP_CRSICE, "Dust Core Mass", &
         RHO_DUST, I_COREMASS, I_DUST, rc, shortname="CRDCOR")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    call CARMAELEMENT_Create(carma, I_ELEM_CRLIQ, I_GRP_CRLIQ, "Water Drop", RHO_W, &
         I_VOLATILE, I_WATER, rc, shortname="CRLIQ")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CRGRP, I_GRP_CRGRP, "Graupel", RHO_I, &
         I_VOLATILE, I_ICE, rc, shortname="CRGRP")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! LW (Sep 2022): modified for hetero. nuc.

    call CARMAELEMENT_Create(carma, I_ELEM_CRDUST, I_GRP_CRDUST, "Dust", RHO_DUST, &
         I_INVOLATILE, I_DUST, rc, shortname="CRDUST")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Define the Solutes

    call CARMASOLUTE_Create(carma, I_SOL_CRH2SO4, "Sulfuric Acid", 2, &
         98._f, 1.38_f, rc, shortname="CRH2SO4")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMASOLUTE_Create failed.')
    
    ! Define the Gases

    call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, &
         I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, rc, shortname="Q", ds_threshold=-0.2_f)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')
    
    ! Define the Processes
    
    ! Detrained Ice, Aged
    call CARMA_AddGrowth(carma, I_ELEM_CRDICE, I_GAS_H2O, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddNucleation(carma, I_ELEM_CRDICE, I_ELEM_CRLIQ, I_ICEMELT, &
         -latice*1e4_f, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_CRDICE, I_GRP_CRDICE, I_GRP_CRDICE, &
         I_COLLEC_DATA, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    ! In-Situ Ice
    call CARMA_AddGrowth(carma, I_ELEM_CRSICE, I_GAS_H2O, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    ! NOTE: For now, assume the latent heat for nucleation is the latent of of fusion of
    ! water, using the CAM constant (scaled from J/kg to erg/g).
    !
    ! NOTE: Since the sulfates are not seen as part of the water/energy budget in CAM, don't
    ! include any latent heat from the freezing of the sulfate liquid. The latent heat of
    ! the gas associated with nucleation is accounted for.
    !
    ! NOTE: The MAM3 vs. MAM4 changes should be double checked and the code should be set
    ! up to handle both based upon the namelist setting. I don't think that is currently
    ! the case.
    !
    ! NOTE: In CAM the MAM version are usually referred to as MAM3 and MAM4. The number of
    ! modes can be queried from rad_constiuents, so just need to say it is model.
    if (carma_sulfate_method .eq. "fixed") then
      call CARMA_AddNucleation(carma, I_ELEM_CRCN, I_ELEM_CRCORE, I_AERFREEZE + I_AF_KOOP_2000, &
        0._f, rc, igas=I_GAS_H2O, ievp2elem=I_ELEM_CRCN)
      if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')
    else if (carma_sulfate_method .eq. "modal") then
      call CARMA_AddNucleation(carma, I_ELEM_CRCN, I_ELEM_CRCORE, I_AERFREEZE + I_AF_KOOP_2000, &
        0._f, rc, igas=I_GAS_H2O, ievp2elem=I_ELEM_CRCN)
      if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

      call CARMA_AddNucleation(carma, I_ELEM_CRDUST, I_ELEM_CRDCOR, I_HETNUC, 0._f, rc, &
        igas=I_GAS_H2O, ievp2elem=I_ELEM_CRDUST)
      if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')   
    end if
    
    call CARMA_AddNucleation(carma, I_ELEM_CRSICE, I_ELEM_CRLIQ, I_ICEMELT, -latice*1e4_f, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_CRSICE, I_GRP_CRSICE, I_GRP_CRSICE, I_COLLEC_DATA, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

!    call CARMA_AddNucleation(carma, I_ELEM_CRLIQ, I_ELEM_CRSICE, I_IMMNUCL, 0._f, rc)
!    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    ! Water Drop
    call CARMA_AddGrowth(carma, I_ELEM_CRLIQ, I_GAS_H2O, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddNucleation(carma, I_ELEM_CRLIQ, I_ELEM_CRDICE, I_DROPFREEZE, latice*1e4_f, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_CRLIQ, I_GRP_CRLIQ, I_GRP_CRLIQ, I_COLLEC_DATA, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    
    !+++ for Graupel

    call CARMA_AddGrowth(carma, I_ELEM_CRGRP, I_GAS_H2O, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth for Graupel failed.')
    
    call CARMA_AddCoagulation(carma, I_GRP_CRLIQ, I_GRP_CRDICE, I_GRP_CRGRP, I_COLLEC_DATA, rc) 
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation for Graupel failed.')
    
    call CARMA_AddCoagulation(carma, I_GRP_CRLIQ, I_GRP_CRSICE, I_GRP_CRGRP, I_COLLEC_DATA, rc) 
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation for Graupel failed.')

    !+++ #


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Sep 2022): added for hetero. nuc.

    ! Define the indexing variables for MAM 4 variables

      ! Get dgnum for each mode
      call rad_cnst_get_info(0, nmodes=nmodes)
!      dgnum_idx   = pbuf_get_index('DGNUM')
!      call pbuf_get_field(pbuf, dgnum_idx,    dgnum,    start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )

      !Set up indices for modes/species
      ! Mode indices
      do m = 1, nmodes
        call rad_cnst_get_info(0, m, mode_type=str32)

        select case (trim(str32))
        case ('accum')
            mode_accum_idx = m
        case ('aitken')
            mode_aitken_idx = m
        case ('coarse')
            mode_coarse_idx = m
        end select
      end do

     ! Species indice (only get coarse mode until everything is working correctly)
      ! Sulfate
      call rad_cnst_get_info(0, mode_accum_idx, nspec=nspec)
      do nn = 1, nspec
         call rad_cnst_get_info(0, mode_accum_idx, nn, spec_type=str32)
         select case (trim(str32))
         case ('sulfate')
            so4_accum_idx = nn
         end select
      end do


      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do nn = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, nn, spec_type=str32)
         select case (trim(str32))
         case ('sulfate')
            so4_coarse_idx = nn
         end select
      end do

      call rad_cnst_get_info(0, mode_aitken_idx, nspec=nspec)
      do nn = 1, nspec
         call rad_cnst_get_info(0, mode_aitken_idx, nn, spec_type=str32)
         select case (trim(str32))
         case ('sulfate')
            so4_aitken_idx = nn
         end select
      end do

      ! Sea Salt
      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do nn = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, nn, spec_type=str32)
         select case (trim(str32))
         case ('seasalt')
            sslt_coarse_idx = nn
         end select
      end do

      call rad_cnst_get_info(0, mode_accum_idx, nspec=nspec)
      do nn = 1, nspec
         call rad_cnst_get_info(0, mode_accum_idx, nn, spec_type=str32)
         select case (trim(str32))
         case ('seasalt')
            sslt_accum_idx = nn
         end select
      end do

      ! Dust  
      call rad_cnst_get_info(0, mode_accum_idx, nspec=nspec)
      do nn = 1, nspec
         call rad_cnst_get_info(0, mode_accum_idx, nn, spec_type=str32)
         select case (trim(str32))
         case ('dust')
            dust_accum_idx = nn
         end select
      end do

      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do nn = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, nn, spec_type=str32)
         select case (trim(str32))
         case ('dust')
            dust_coarse_idx = nn
         end select
      end do

      ! Primary Organics
      call rad_cnst_get_info(0, mode_accum_idx, nspec=nspec)
      do nn = 1, nspec
         call rad_cnst_get_info(0, mode_accum_idx, nn, spec_type=str32)
         select case (trim(str32))
         case ('p-organic')
            pom_accum_idx = nn
         end select
      end do

      ! Secondary Organics
      call rad_cnst_get_info(0, mode_accum_idx, nspec=nspec)
      do nn = 1, nspec
         call rad_cnst_get_info(0, mode_accum_idx, nn, spec_type=str32)
         select case (trim(str32))
         case ('s-organic')
            soa_accum_idx = nn
         end select
      end do

      ! Black Carbon
      call rad_cnst_get_info(0, mode_accum_idx, nspec=nspec)
      do nn = 1, nspec
         call rad_cnst_get_info(0, mode_accum_idx, nn, spec_type=str32)
         select case (trim(str32))
         case ('black-c')
            bc_accum_idx = nn
         end select
      end do

      ! Check that the proper mode species were found
      if ( (so4_coarse_idx == -1) .or. (so4_accum_idx == -1)) then
         write(LUNOPRT,*) routine//': ERROR required so4 mode-species type not found - indicies:', &
            so4_coarse_idx, so4_accum_idx, so4_aitken_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
      if ( (sslt_coarse_idx == -1) .or. (sslt_accum_idx == -1)) then
         write(LUNOPRT,*) routine//': ERROR required sea salt mode-species type not found - indicies:', &
            sslt_coarse_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
      if ( (dust_coarse_idx == -1) .or. (dust_accum_idx == -1)) then
         write(LUNOPRT,*) routine//': ERROR required dust mode-species type not found - indicies:', &
            dust_coarse_idx, dust_accum_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
      if ( pom_accum_idx == -1) then
         write(LUNOPRT,*) routine//': ERROR required primary organics mode-species type not found - indicies:', &
            pom_accum_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
      if ( soa_accum_idx == -1) then
         write(LUNOPRT,*) routine//': ERROR required secondary organics mode-species type not found - indicies:', &
            soa_accum_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
      if ( bc_accum_idx == -1) then
         write(LUNOPRT,*) routine//': ERROR required black carbon mode-species type not found - indicies:', &
            bc_accum_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

       ! Get mode distribution sigma and diameter bounds
!write(LUNOPRT,*) ' get coarse mode props'
!       call rad_cnst_get_mode_props(0, mode_coarse_idx, sigmag=sigmag_coarse, dgnumlo=dgnumlo_co, dgnumhi=dgnumhi_co)
!write(LUNOPRT,*) ' get accumulation mode props'
!       call rad_cnst_get_mode_props(0, mode_accum_idx, sigmag=sigmag_accum) !, dgnumlo=dgnumlo_ac, dgnumhi=dgnumhi_ac)
!write(LUNOPRT,*) ' get aitken mode props'
!       call rad_cnst_get_mode_props(0, mode_aitken_idx, sigmag=sigmag_aitken)!,  dgnumlo=dgnumlo_at, dgnumhi=dgnumhi_at)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    return
  end subroutine CARMA_DefineModel
  

  !! Defines all the CARMA components (groups, elements, solutes and gases) and
  !! process (coagulation, growth, nucleation) that will be part of the
  !! microphysical model.
  !!
  !! NOTE: This currently won't work with the miacmic loop in physpkg, so it
  !! either needs to be fixed if we want to detrain ice (and possibly liquid).
  !! However, the best thing may just be to let CLUBB do the detrainment turn
  !! this off.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen (CB)
  !!
  !!  @version 2020 May
  !!  @author  Jamison A. Smith (JAS)
  subroutine CARMA_Detrain(carma, cstate, cam_in, dlf, det_cnd, det_ice, state, icol, dt, rc)
    use camsrfexch,         only: cam_in_t
    use physconst,          only: latice, latvap, cpair

    implicit none
    
    type(carma_type), intent(in)         :: carma !! the carma object
    type(carmastate_type), intent(inout) :: cstate !! the carma state object
    type(cam_in_t),  intent(in)          :: cam_in !! surface input
    real(r8), intent(in)                 :: dlf(pcols, pver) !! Detraining cld H20 from convection (kg/kg/s)
    real(r8), intent(out)                :: det_cnd(pcols) !! Detrained condensate (m/s)
    real(r8), intent(out)                :: det_ice(pcols) !! Detraining cld H20 from convection (m/s)
    type(physics_state), intent(in)      :: state !! physics state variables
    integer, intent(in)                  :: icol !! column index
    real(r8), intent(in)                 :: dt !! time step (s)
    integer, intent(out)                 :: rc !! return code, negative indicates failure
    
    real(kind=f) :: t(pver)             ! temperature (K)
    real(kind=f) :: mmr_ice(NBIN, pver) ! ice mass mixing ratio (kg/kg)
    real(kind=f) :: mmr_liq(NBIN, pver) ! liquid mass mixing ratio (kg/kg)
    real(kind=f) :: r_ice(NBIN)         ! ice radius bins (cm)
    real(kind=f) :: r_liq(NBIN)         ! liquid radius bins (cm)

    real(kind=f) :: ice_fraction  ! fraction of detrained condensate that is ice
    real(kind=f) :: mass_liq      ! detrainment rate of liquid (kg/kg/s)
    real(kind=f) :: mass_ice      ! detrainment rate of ice (kg/kg/s)
    real(kind=f) :: mass_dlf      ! detrained mass (m/s)

    integer      :: k             ! vertical index
    integer      :: ibin          ! bin index
    integer      :: itemp         ! temperature index

    real(r8) :: iceMass(pver)       ! ice mass mixing ratio (kg/kg)
    real(r8) :: iceNumber(pver)     ! ice number mixing ratio (#/kg)
    real(r8) :: snowSurface         ! snow on surface (kg/m2)
    real(r8) :: waterMass(pver)     ! ice mass mixing ratio (kg/kg)
    real(r8) :: waterNumber(pver)   ! ice number mixing ratio (#/kg)
    real(r8) :: rainSurface         ! rain on surface (kg/m2)
    
    logical  :: do_thermo           ! do thermodynamics?


    rc = rc_ok    ! Default return code.
    
    call CARMA_Get(carma, rc, do_thermo = do_thermo)
    if (rc < rc_OK) call endrun('carma_model_mod::Carma_Detrain - Carma_Get failed.')
    
    ! Put all of the detraining cloud water from convection into the large
    ! scale cloud. Put detraining cloud water into liq and ice based on
    ! temperature partition
    call CARMAGROUP_Get(carma, I_GRP_CRDICE, rc, r=r_ice(:))
    if (rc < rc_ok) call endrun( 'carma_model_mod::arma_Detrain - CARMAGROUP_Get failed, I_GRP_CRDICE')

    call CARMAGROUP_Get(carma, I_GRP_CRLIQ, rc, r=r_liq(:))
    if (rc < rc_ok) call endrun( 'carma_model_mod::Carma_Detrain - CarmaGroup_Get failed, I_GRP_CRLIQ')
    
    ! Account for the reserved ice that is being detrained in the precipitation.
    call CARMASTATE_GetState(cstate, rc, t=t)
    if (rc < rc_ok) call endrun ( 'carma_model_mod::Carma_Detrain - CarmaState_GetState failed.')
    
    !  Determine the amount of detrainment that could be used to saturate
    ! the atmosphere with respect to liquid. For GCM scales, assume that three
    ! things (could -JAS) happen to detrained condensate:
    !
    !   1) large particles will fallout as snow or rain
    !   2) will be converted to vapor
    !   3) will remain as ice
    !
    ! Because of the large scales of the GCM and because this is a stratiform
    ! parameterization, a lot of the condensate that hasn't fallen out will
    ! increase the humidity (i.e. detrained anvil evaporates or falls out
    ! entirely with 100 km of the convection).
    mmr_ice(:, :)  = 0._f
    mmr_liq(:, :)  = 0._f

    det_cnd(icol) = 0._f
    det_ice(icol) = 0._f
    
    do k = 1, pver
    
      ! Remove amount being detrained from rliq and prec_str.
      !
      ! NOTE: Since the clouds are now included in a loop, you can't
      ! change rliq or dlf as you want to have the same process rates
      ! for each substep. Also, CLUBB has already detrained the
      ! condensate, but it is not included in CMELIQ and CMEICE, so
      ! we will detrain it again into the bin structure.
      mass_dlf = dlf(icol, k) * (state%pdel(icol, k) / gravit)

      if (t(k) > 268.15_f) then
        ice_fraction = 0.0_f
      else if (t(k) < 238.15_f) then
        ice_fraction = 1.0_f
      else
        ice_fraction = (268.15_f - t(k)) / 30._f
      end if
      
      mass_liq  = dlf(icol, k) * (1._f - ice_fraction)
      mass_ice  = dlf(icol, k) * ice_fraction

      det_cnd(icol) = det_cnd(icol) + mass_dlf / 1000._f
      det_ice(icol) = det_ice(icol) + mass_dlf * ice_fraction / 1000._f

      ! Calculate the detrainment of ice and liquid into the appropriate
      ! CARMA bins.
      !      
      ! Scale the size based on whether the surface is land or ocean. This
      ! assumes that there are more aerosols over land, reducing the detrainment
      ! size. This is similar to the c0_lnd and c0_ocn parameter split done in
      ! the convective parameterization.
      !
      ! NOTE: This should really be tied to aerosol amount, not land fraction.
      do ibin = 1, NBIN
      
        ! Assume detrained cloud water is monodisperse.
        if (r_liq(ibin) >= r_dliq_ocn) then
          mmr_liq(ibin, k) = mmr_liq(ibin, k) + mass_liq * dt * (1._f - cam_in%landfrac(icol))
          exit
        end if
      end do
      
      do ibin = 1, NBIN

        ! Assume detrained cloud water is monodisperse.
        if (r_liq(ibin) >= r_dliq_lnd) then
          mmr_liq(ibin, k) = mmr_liq(ibin, k) + mass_liq * dt * cam_in%landfrac(icol)
          exit
        end if
      end do

      ! Detrain cloud ice into the bins according to the predefined distribution.
      itemp = max(-max(MIN_DTEMP, nint(t(k)-T0)), 0) + 1

      do ibin = 1, NBIN
   
        ! Detrain using a size distribution (log-normal in mass). The table has
        ! already been setup during initialization indicating the fraction of
        ! the mass that goes into each bin.
        mmr_ice(ibin, k) = mmr_ice(ibin, k) + dice_bin_fraction(ibin, itemp) * mass_ice * dt
      end do
           
      ! Account for latent heat release during freezing. By default the
      ! detrained condensate is assumed to be liquid for energy balance.
      !
      ! NOTE: Since CLUBB has already detrained, this may already be accounted for.
      ! Need to check energy conservation.
      !t(k) = t(k) + (mass_ice * latice * dt / cpair)
    end do  ! k


    do ibin = 1, NBIN
      call CARMASTATE_SetDetrain(cstate, I_ELEM_CRLIQ, ibin, mmr_liq(ibin, :), rc)
      if (rc < rc_ok) call endrun ( 'carma_model_mod::Carma_Detrain - CarmaState_SetDetrain failed, I_ELEM_CRLIQ')
  
      call CARMASTATE_SetDetrain(cstate, I_ELEM_CRDICE, ibin, mmr_ice(ibin, :), rc)
      if (rc < rc_ok) call endrun ( 'carma_model_mod::Carma_Detrain - CarmaState_SetDetrain failed, I_ELEM_CRDICE')
    end do
    
    ! Calculate the total column of detrain condensate and ice
    do ibin = 1, NBIN
    end do    
        
    if (do_thermo) then
      call CARMASTATE_SetState(cstate, rc, t(:))
    end if
    
    ! Check for total water conservation by CARMA.
    if (carma_do_mass_check2) then
      call CARMA_GetTotalWaterAndRain(carma, cstate, waterMass, waterNumber, rainSurface, rc)
      call CARMA_GetTotalIceAndSnow(carma, cstate, iceMass, iceNumber, snowSurface, rc)
      call CARMA_CheckMassAndEnergy(carma, cstate, "Carma_Detrain", state, icol, dt, dlf, waterMass, rainSurface, iceMass, snowSurface, rc)    
    end if
    
    return
  End Subroutine CARMA_Detrain


  !! For diagnostic groups, sets up up the CARMA bins based upon the CAM state.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen (CB) 
  !!  @version 2019 November
  !!  @author  Jamison A. Smith (JAS)
 
  !subroutine CARMA_DiagnoseBins (carma, cstate, state, pbuf, icol, dt &
  !                              , rc, rliq, prec_str, snow_str)
  subroutine CARMA_DiagnoseBins(carma, cstate, state, pbuf, icol, dt, rc, &
            dlf, prec_str, snow_str, so4_ac_frac, so4_co_frac, dst_ac_frac, dst_co_frac)
            ! (LW: modified for hetero. nuc., Sep 2022)

    use time_manager,   only: is_first_step
    use micro_mg_utils, only: size_dist_param_basic, size_dist_param_liq, &
         mg_ice_props, mg_liq_props, mg_rain_props, mg_snow_props
    use constituents,     only: cnst_get_ind

    implicit none
    
    type(carma_type), intent(in)          :: carma        !! the carma object
    type(carmastate_type), intent(inout)  :: cstate       !! the carma state object
    type(physics_state), intent(in)    :: state        !! physics state variables
    type(physics_buffer_desc), pointer    :: pbuf(:)      !! physics buffer
    integer, intent(in)                   :: icol         !! column index
    real(r8), intent(in)                  :: dt           !! time step
    integer, intent(out)                  :: rc           !! return code, negative indicates failure
    real(r8), intent(in), optional        :: dlf(pcols, pver) !! detrained condensation rate (kg/kg/s)

    ! JAS: I do not want the mass fix applied, so this subroutine will no longer modify prec_str

    real(r8), intent(inout), optional     :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional     :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(out), optional       :: so4_ac_frac(pcols,pver,NBIN) !!(LW: added for hetero. nuc.)
    real(r8), intent(out), optional       :: so4_co_frac(pcols,pver,NBIN) !!(LW: added for hetero. nuc.)
    real(r8), intent(out), optional       :: dst_ac_frac(pcols,pver,NBIN) !!(LW: added for hetero. nuc.)
    real(r8), intent(out), optional       :: dst_co_frac(pcols,pver,NBIN) !!(LW: added for hetero. nuc.)

    real(r8)                              :: mu(pver)       ! spectral width parameter of droplet size distr
    real(r8)                              :: lambda(pver)   ! slope of cloud liquid size distr

    ! JAS:  MG 2.0 has more types of prognostic condensate, so Liquid is now
    !      the Sum of CLDLIQ and RAINQM, and Ice is now the Sum of CLDICE and
    !      SNOWQM.  mmr has been replaced by mmr1 and mmr2, so that each
    !      Liquid mmr may be stored and the sum may be computed.  Same for Ice

    real(r8)                              :: mmr(NBIN,pver)  ! element's mass mixing ratio
    real(r8)                              :: mmr1(NBIN,pver) ! element's mass mixing ratio
    real(r8)                              :: mmr2(NBIN,pver) ! element's mass mixing ratio

    real(kind=f)                          :: r(NBIN)      ! bin mean radius
    real(kind=f)                          :: dr(NBIN)     ! bin radius width
    real(kind=f)                          :: rmass(NBIN)  ! bin mass
    
    integer                               :: igroup       ! group index
    integer                               :: ielem        ! element index
    integer                               :: ibin         ! bin index
    integer                               :: k            ! vertical index

    ! This buffer exists purely to work around the fact that "state" is
    ! intent(in), but the size_dist_param function will try to change the
    ! input number concentrations.

    real(r8)                              :: limNumber(pver)

    real(r8)                              :: iceMass(pver)      ! ice mass mixing ratio (kg/kg)
    real(r8)                              :: iceNumber(pver)    ! ice number mixing ratio (#/kg)
    real(r8)                              :: snowMass(pver)     ! snow mass mixing ratio (kg/kg)
    real(r8)                              :: snowNumber(pver)   ! snow number (#/kg)
    real(r8)                              :: snowSurface        ! snow on surface (kg/m2)
    real(r8)                              :: carma_ice          ! total cldice from CARMA bins (kg/kg)
    real(r8)                              :: waterMass(pver)    ! ice mass mixing ratio (kg/kg)
    real(r8)                              :: waterNumber(pver)  ! ice number mixing ratio (#/kg)
    real(r8)                              :: rainSurface        ! rain on surface (kg/m2)
    real(r8)                              :: carma_water        ! total cldliq from CARMA bins (kg/kg)
    real(r8)                              :: diff

    ! Aerosol size distribution

    real(r8), parameter                   :: n    = 100._r8     ! concentration (cm-3) 
    real(r8), parameter                   :: r0   = 2.5e-6_r8   ! mean radius (cm)
    real(r8), parameter                   :: rsig = 1.5_r8      ! distribution width
    
    real(r8)                              :: arg1(NBIN)
    real(r8)                              :: arg2(NBIN)
    real(r8)                              :: rhop(NBIN)         ! particle mass density (kg/m3)
    real(r8)                              :: totalrhop          ! total particle mass density (kg/m3)
    real(kind=f)                          :: rhoa_wet(pver)     ! air density (g/cm3)

    ! vertical integral of liquid not yet in q(ixcldliq)

    real(r8)                              :: rliq_new(pcols)    

    integer                               :: LUNOPRT
    logical                               :: do_print
    real(r8)                              :: lat
    real(r8)                              :: lon

    real(r8), pointer, dimension(:, :)    :: sulf               ! last saturation wrt ice
    integer                               :: lchnk              ! chunk identifier (LW: added for hetero. nuc.)
    integer                               :: itim_old

    integer                               :: icnst_q 
  
    character(len=8)                      :: c_name             ! constituent name
    
    real(r8), pointer, dimension(:,:)     :: activationrate     ! #/kg/s
    integer                               :: ienconc                        

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Sep 2022): added for hetero. nuc.
!    real(r8), pointer                     :: so4_mmr_co(:,:)    ! interstitial sulfate mmr, coarse mose (kg/kg)
!    real(r8), pointer                     :: so4_mmr_ac(:,:)    ! interstitial sulfate mmr, accumulation mode (kg/kg)
!    real(r8), pointer                     :: sslt_mmr_co(:,:)    ! interstitial sea salt mmr, coarse mode (kg/kg)
!    real(r8), pointer                     :: dust_mmr_co(:,:)   ! interstitial dust mmr, coarse mode (kg/kg)
!    real(r8), pointer                     :: num_coarse(:,:)    ! number mixing ratio of coarse mode (#/kg)

    ! Indexing variables for modes and species (3 species so4, dust, sea salt)
!    integer :: m, nn, ii, l

    ! Aerosol constituents info
!    character(len=32)     :: str32
!    character(len=32)   :: tmpname
!    character(len=32)   :: tmpname_cw
!    character(len=*), parameter :: routine = 'CARMA_DiagnoseBins'
!    integer               :: nspec                       ! number of chemical species
!    integer               :: nmodes                       ! number of aerosol modes
    integer               :: dgnum_idx   = -1
    real(r8), pointer     :: dgnum(:,:,:)                        ! aerosol mode dry diameter (m)
    real(r8)              :: dgnum_accum
    real(r8)              :: dgnum_coarse

    ! Aerosol distribution variables
    real(r8)              :: so4_arg1(NBIN)
    real(r8)              :: so4_arg2(NBIN)
    real(r8)              :: so4_arg3(NBIN)
    real(r8)              :: so4_arg4(NBIN)
    real(r8)              :: dust_arg1(NBIN)
    real(r8)              :: dust_arg2(NBIN)
    real(r8)              :: dust_arg3(NBIN)
    real(r8)              :: dust_arg4(NBIN)
    real(r8)              :: nmr_co_bin(NBIN,pver)
    real(r8)              :: nmr_ac_bin(NBIN,pver)
    real(r8)              :: mmr_frac(NBIN,pver)
    real(r8)              :: mmr_sum(pver)
    real(r8)              :: nmr_so4
    real(r8)              :: nmr_dust
    real(r8)              :: mmr_tot_ac
    real(r8)              :: mmr_tot_co
    real(r8)              :: bin_frac_ac(NBIN,pver)
    real(r8)              :: bin_frac_co(NBIN,pver)
    real(r8)              :: mmr_ac(NBIN,pver)
    real(r8)              :: mmr_co(NBIN,pver)
    real(r8)              :: so4_frac
    real(r8)              :: dust_frac
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  1 format(/,'CARMA_DiagnoseBins::ERROR - CAM ice mass conservation error, icol=',i4,', iz=',i4,',lat=',&
              f7.2,',lon=',f7.2,',cam=',e17.10,',carma=',e17.10,',rer=',e10.3)
  2 format(/,'CARMA_DiagnoseBins::ERROR - CAM liquid mass conservation error, icol=',i4,', iz=',i4,',lat=',&
              f7.2,',lon=',f7.2,',cam=',e17.10,',carma=',e17.10,',rer=',e10.3)

    ! Default return code.

    rc = RC_OK

    mmr = 0.0_f
    mmr1 = 0.0_f
    mmr2 = 0.0_f
    
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    
    ! Get the air density.

    call CARMASTATE_GetState(cstate, rc, rhoa_wet=rhoa_wet)
    if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_GetState failed.')

    ! Aerosols
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Sep 2022): added for hetero. nuc.
!    ! Before reading in individual species, set up indexing

      ! Get dgnum for each mode
!      call rad_cnst_get_info(0, nmodes=nmodes)
      dgnum_idx   = pbuf_get_index('DGNUM')
      call pbuf_get_field(pbuf, dgnum_idx,    dgnum,    start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
!       write(*,*) 'dgnum_idx:'
!       write(*,*) dgnum_idx

!       ! Get mode distribution sigma
       call rad_cnst_get_mode_props(0, mode_coarse_idx, sigmag=sigmag_coarse)
!       write(*,*) 'check mode_coarse_idx:'
!       write(*,*) mode_coarse_idx
       call rad_cnst_get_mode_props(0, mode_accum_idx, sigmag=sigmag_accum)
!       write(*,*) 'mode_accum_idx:'
!       write(*,*) mode_accum_idx
!
        ! Get mode number mixing ratio
        call rad_cnst_get_mode_num(0, mode_coarse_idx, 'a', state, pbuf, num_coarse)
        call rad_cnst_get_mode_num(0, mode_accum_idx, 'a', state, pbuf, num_accum)

        ! Get aerosol species mmr
        ! Sulfates
        call rad_cnst_get_aer_mmr(0, mode_coarse_idx, so4_coarse_idx, 'a', state, pbuf, so4_mmr_co)
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, so4_accum_idx, 'a', state, pbuf, so4_mmr_ac)
!        write(*,*) 'so4_coarse_idx:'
!        write(*,*) so4_coarse_idx
!        write(*,*) 'so4_accum_idx:'
!        write(*,*) so4_accum_idx
        
        ! Sea Salt
        call rad_cnst_get_aer_mmr(0, mode_coarse_idx, sslt_coarse_idx, 'a', state, pbuf, sslt_mmr_co)
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, sslt_accum_idx, 'a', state, pbuf, sslt_mmr_ac)
!        write(*,*) 'sslt_coarse_idx:'
!        write(*,*) sslt_coarse_idx
!        write(*,*) 'sslt_accum_idx:'
!        write(*,*) sslt_accum_idx
        ! Dust
        call rad_cnst_get_aer_mmr(0, mode_coarse_idx, dust_coarse_idx, 'a', state, pbuf, dust_mmr_co)     
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, dust_accum_idx, 'a', state, pbuf, dust_mmr_ac)
!        write(*,*) 'dust_coarse_idx:'
!        write(*,*) dust_coarse_idx
!        write(*,*) 'dust_accum_idx:'
!        write(*,*) dust_accum_idx
        ! Black Carbon
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, bc_accum_idx, 'a', state, pbuf, bc_mmr_ac)
!        write(*,*) 'bc_accum_idx:'
!        write(*,*) bc_accum_idx
        ! Primary Organics
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, pom_accum_idx, 'a', state, pbuf, pom_mmr_ac)
!        write(*,*) 'pom_accum_idx:'
!        write(*,*) pom_accum_idx
        ! Secondary Organics
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, soa_accum_idx, 'a', state, pbuf, soa_mmr_ac)
!        write(*,*) 'soa_accum_idx:'
!        write(*,*) soa_accum_idx


! Get aerosol (sulfate)  mmr
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    igroup = 1
    ielem  = 1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Sep 2022): added for hetero. nuc.
    ! Use the MAM 4 aerosols
    !
    ! NOTE: Should query rad_constituents for number of modes and set up the code
    ! so that it works for either 3 or 4 modes. If you don't want to support MAM3
    ! then can just check for MAM4 and fail the run if their aren't 4 modes.
    if (carma_sulfate_method == "modal") then 
      lchnk = state%lchnk
 
      ! Get the correct aerosol group
      call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')

     ! Assume external mixing for coarse mode aerosols
     do k = 1, pver

        ! An issue arose in the first time step where dgnum = 0 because the physic's buffer
        ! is yet to be created causing the model to break. Try to prevent that here
        dgnum_accum = dgnum(icol,k,mode_accum_idx)

        if (is_first_step()) then
          if (dgnum(icol,k,mode_accum_idx) .eq. 0.0_f) then
             dgnum_accum = 1.0e-7_f
          end if
        end if

        dgnum_coarse = dgnum(icol,k,mode_coarse_idx)

        if (is_first_step()) then
          if (dgnum(icol,k,mode_coarse_idx) .eq. 0.0_f) then
             dgnum_coarse = 1.0e-7_f
          end if
        end if

       ! Use the lognormal distribution (i.e. equation 7.34, page 421 form Seinfield and Pandis (1998)) to get bin dist.
       ! Note that MAM 4 also assumes lognormal distributions for aerosols
         nmr_ac_bin(:,k) = 0._f
  !       mmr_frac(:,k) = 0._f
         nmr_so4 = 0._f
  !       mmr_sum(k) = so4_mmr_co(icol,k) + so4_mmr_ac(icol,k)

         ! Accumulation Mode
         mmr_tot_ac = so4_mmr_ac(icol,k) + dust_mmr_ac(icol,k) + sslt_mmr_ac(icol,k) + pom_mmr_ac(icol,k) + soa_mmr_ac(icol,k) + bc_mmr_ac(icol,k)
         so4_frac = so4_mmr_ac(icol,k) / mmr_tot_ac
         nmr_so4 = nmr_so4 + num_accum(icol,k) * so4_frac

         so4_arg1(:) = nmr_so4 * dr(:) / (sqrt(2._f * PI) * r(:) * log(sigmag_accum))
         so4_arg2(:) = -((log(r(:)) - log(dgnum_accum * 0.5_r8 * 1.0e2_r8))**2) / (2._f*(log(sigmag_accum))**2)

         nmr_ac_bin(:,k) = nmr_ac_bin(:,k) + so4_arg1(:) * exp(so4_arg2(:))


         ! Coarse Mode
         nmr_co_bin(:,k) = 0._f
         nmr_so4 = 0._f

         mmr_tot_co = so4_mmr_co(icol,k) + dust_mmr_co(icol,k) + sslt_mmr_co(icol,k)
         so4_frac = so4_mmr_co(icol,k) / mmr_tot_co
         nmr_so4 = nmr_so4 + num_coarse(icol,k) * so4_frac

         so4_arg3(:) = nmr_so4 * dr(:) / (sqrt(2._r8 * PI) * r(:) * log(sigmag_coarse)) 
         so4_arg4(:) = -((log(r(:)) - log(dgnum_coarse * 0.5_r8 * 1.0e2_r8))**2) / (2._r8*(log(sigmag_coarse))**2)

         nmr_co_bin(:,k) = nmr_co_bin(:,k) +  so4_arg3(:) * exp(so4_arg4(:))

         ! Calculate the combined nmr fraction per bin which will be used to calculate
         ! how much mass to put into each bin
!         mmr_frac(:,k) = mmr_frac(:,k) + nmr_co_bin(:,k) + nmr_ac_bin(:,k)
   
         bin_frac_ac(:,k) =  nmr_ac_bin(:,k) / sum(nmr_ac_bin(:,k))
         bin_frac_co(:,k) =  nmr_co_bin(:,k) / sum(nmr_co_bin(:,k))

         mmr_ac(:,k) = so4_mmr_ac(icol,k) * bin_frac_ac(:,k)
         mmr_co(:,k) = so4_mmr_co(icol,k) * bin_frac_co(:,k)

  !        ! We need to keep track of how much aerosol came from each mode in order to properly return things to MAM later
  !        do ibin = 1, NBIN
  !          so4_ac_frac(icol,k,ibin) = nmr_ac_bin(ibin,k) / mmr_frac(ibin,k)
  !          so4_co_frac(icol,k,ibin) = nmr_co_bin(ibin,k) / mmr_frac(ibin,k)
  !        end do
!
!         mmr_frac(:,k) = mmr_frac(:,k) / sum(mmr_frac(:,k))

      end do  !do k = 1, pver
    end if  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    ! Use a fixed aerosols distribution.

    if ((carma_sulfate_method == "fixed") .or. (carma_sulfate_method == "bulk")) then
    
      call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')
    
      arg1(:) = n * dr(:) / (sqrt(2._f*PI) * r(:) * log(rsig))
      arg2(:) = -((log(r(:)) - log(r0))**2) / (2._f*(log(rsig))**2)

      ! kg/m3
      rhop(:)   = arg1(:) * exp(arg2(:)) * rmass(:) * 1e6_f / 1e3_f

    
      if (carma_sulfate_method == "bulk") then
        totalrhop = sum(rhop(:))

        ! Get the index for the prescribed sulfates. This gives the mmr that should be
        ! present at this location. Use this to scale the size distribution that CARMA
        ! will generate.
        lchnk = state%lchnk  !(LW: added for hetero. nuc.)
        itim_old  = pbuf_old_tim_idx()

        call pbuf_get_field(pbuf, pbuf_get_index('sulf'), sulf, (/1,1,itim_old/),(/pcols,pver,1/))
      end if
    end if
    
    do ibin = 1, NBIN
    
      ! Use a fixed mixing ratio.

      if (carma_sulfate_method == "fixed") then
        mmr(ibin, :) = rhop(ibin) / rhoa_wet(:)
      end if
      
      ! Since bulk aerosols don't have a size distribution, use the fixed
      ! distribution for the shape of the distribution, but scale the total
      ! mass to the prescribed value.

      if (carma_sulfate_method == "bulk") then
        mmr(ibin, :) = rhop(ibin) / totalrhop * sulf(icol, :)
      end if

      ! Use the CRCNxx fields from a special prescribed aerosol file that has
      ! results from a CARMA simulation of sulfates. This will set the magnitude
      ! and the size distribution.

      if (carma_sulfate_method == "carma") then

        ! Get the index for the prescribed sulfates.
        lchnk = state%lchnk  !(LW: added for hetero. nuc.)
        itim_old  = pbuf_old_tim_idx()
        write(c_name, '(A, I2.2)') "CRCN", ibin

        call pbuf_get_field(pbuf, pbuf_get_index(c_name), sulf, (/1,1,itim_old/),(/pcols,pver,1/))
        mmr(ibin, :) = sulf(icol, :)
      end if

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Sep 2022): added for hetero. nuc.
    ! Use MAM 4 sulfate size distribution
      if (carma_sulfate_method == "modal") then
        do k = 1, pver
         if ((dgnum(icol,k,mode_accum_idx) .eq. 0._r8) .and. (dgnum(icol,k,mode_coarse_idx) .eq. 0._r8)) then
           mmr(ibin,k) = 0.0_f
         else
!           mmr(ibin,k) = mmr_sum(k) * mmr_frac(ibin,k)
            mmr(ibin,k) = mmr_ac(ibin,k) + mmr_co(ibin,k)

            ! We need to keep track of how much aerosol came from each mode in order to properly return things to MAM later
             so4_ac_frac(icol,k,ibin) = mmr_ac(ibin,k) / mmr(ibin,k)
             so4_co_frac(icol,k,ibin) = mmr_co(ibin,k) / mmr(ibin,k)
         end if
        end do
      end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(ibin, :), rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_SetBin failed.')
    end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Sep 2022): added for hetero. nuc.

!   if (carma_sulfate_method == "modal") then
!      mmr_sulf_org(icol,:) = sum(mmr, DIM=1)
!   end if
 
!do k = 1, pver
!write(LUNOPRT,*) ' '
!write(LUNOPRT,*) ' orig so4 mmr =', mmr_sum(k)
!write(LUNOPRT,*) ' new so4 mmr =', sum(mmr(:,k))
!end do

   ! Aerosol (dust)
   if (carma_sulfate_method == "modal") then
    igroup = I_GRP_CRDUST
    ielem  = I_ELEM_CRDUST

    ! Use the MAM 4 aerosols
      ! Get the correct aerosol group
      call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')

     ! Assume external mixing for coarse mode aerosols
     do k = 1, pver

        ! An issue arose in the first time step where dgnum = 0 because the physic's buffer
        ! is yet to be created causing the model to break. Try to prevent that here
        dgnum_accum = dgnum(icol,k,mode_accum_idx)

        if (is_first_step()) then
          if (dgnum(icol,k,mode_accum_idx) .eq. 0.0_f) then
             dgnum_accum = 1.0e-7_f
          end if
        end if

        dgnum_coarse = dgnum(icol,k,mode_coarse_idx)

        if (is_first_step()) then
          if (dgnum(icol,k,mode_coarse_idx) .eq. 0.0_f) then
             dgnum_coarse = 1.0e-7_f
          end if
        end if


       ! Use the lognormal distribution (i.e. equation 7.34, page 421 form Seinfield and Pandis (1998)) to get bin dist.
       ! Note that MAM 4 also assumes lognormal distributions for aerosols
         nmr_ac_bin(:,k) = 0._f
!         mmr_frac(:,k) = 0._f
         nmr_dust = 0._f
!         mmr_sum(k) = dust_mmr_co(icol,k) + dust_mmr_ac(icol,k)

         ! Accumulation Mode
         mmr_tot_ac = so4_mmr_ac(icol,k) + dust_mmr_ac(icol,k) + sslt_mmr_ac(icol,k) + pom_mmr_ac(icol,k) + soa_mmr_ac(icol,k) + bc_mmr_ac(icol,k)
         dust_frac = dust_mmr_ac(icol,k) / mmr_tot_ac
         nmr_dust = nmr_dust + num_accum(icol,k) * dust_frac

         dust_arg1(:) = nmr_dust * dr(:) / (sqrt(2._f * PI) * r(:) * log(sigmag_accum))
         dust_arg2(:) = -((log(r(:)) - log(dgnum_accum * 0.5_r8 * 1.0e2_r8))**2) / (2._f*(log(sigmag_accum))**2)

         nmr_ac_bin(:,k) = nmr_ac_bin(:,k) +  dust_arg1(:) * exp(dust_arg2(:))

         ! Coarse Mode
         nmr_co_bin(:,k) = 0._f
         nmr_dust = 0._f

         mmr_tot_co = so4_mmr_co(icol,k) + dust_mmr_co(icol,k) + sslt_mmr_co(icol,k)
         dust_frac = dust_mmr_co(icol,k) / mmr_tot_co
         nmr_dust = nmr_dust + num_coarse(icol,k) * dust_frac

         dust_arg3(:) = nmr_dust * dr(:) / (sqrt(2._r8 * PI) * r(:) * log(sigmag_coarse))
         dust_arg4(:) = -((log(r(:)) - log(dgnum_coarse * 0.5_r8 * 1.0e2_r8))**2) / (2._r8*(log(sigmag_coarse))**2)

         nmr_co_bin(:,k) = nmr_co_bin(:,k) + dust_arg3(:) * exp(dust_arg4(:))

         ! Calculate the combined mass fraction per bin
         bin_frac_ac(:,k) =  nmr_ac_bin(:,k) / sum(nmr_ac_bin(:,k))
         bin_frac_co(:,k) =  nmr_co_bin(:,k) / sum(nmr_co_bin(:,k))

         mmr_ac(:,k) = dust_mmr_ac(icol,k) * bin_frac_ac(:,k)
         mmr_co(:,k) = dust_mmr_co(icol,k) * bin_frac_co(:,k)

!write(LUNOPRT,*) ' mmr_frac =', mmr_frac(:,k)
!write(LUNOPRT,*) ' nmmr_ac_bin =', nmr_ac_bin(:,k)
!write(LUNOPRT,*) ' nmr_co_bin =', nmr_co_bin(:,k)
!          ! We need to keep track of how much aerosol came from each mode in order to properly return things to MAM later
!          do ibin = 1, NBIN
!            dst_ac_frac(icol,k,ibin) = nmr_ac_bin(ibin,k) / mmr_frac(ibin,k)
!            dst_co_frac(icol,k,ibin) = nmr_co_bin(ibin,k) / mmr_frac(ibin,k)
!          end do
!
!         mmr_frac(:,k) = mmr_frac(:,k) / sum(mmr_frac(:,k))

!write(LUNOPRT,*) 'so4_ac_frac =', so4_ac_frac(icol,k,:)
!write(LUNOPRT,*) 'so4_co_frac =', so4_co_frac(icol,k,:)
!write(LUNOPRT,*) 'dst_co_frac =', dst_co_frac(icol,k,:)
!write(LUNOPRT,*) 'dst_ac_frac =', dst_ac_frac(icol,k,:)
     end do


     do ibin = 1, NBIN

        ! Use MAM 3 size distribution to set mmr
          do k = 1, pver
           if ((dgnum(icol,k,mode_accum_idx) .eq. 0._r8) .and. (dgnum(icol,k,mode_coarse_idx) .eq. 0._r8)) then
             mmr(ibin,k) = 0.0_f
           else
!             mmr(ibin,k) = mmr_sum(k) * mmr_frac(ibin,k)
            mmr(ibin,k) = mmr_ac(ibin,k) + mmr_co(ibin,k)

            ! We need to keep track of how much aerosol came from each mode in order to properly return things to MAM later
             dst_ac_frac(icol,k,ibin) = mmr_ac(ibin,k) / mmr(ibin,k)
             dst_co_frac(icol,k,ibin) = mmr_co(ibin,k) / mmr(ibin,k)
           end if
          end do

        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(ibin, :), rc)
        if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_SetBin failed.')
    
     end do

!     mmr_dust_org(icol,:) = sum(mmr, DIM=1)

!do k = 1, pver
!write(LUNOPRT,*) ' '
!write(LUNOPRT,*) ' orig dust mmr =', mmr_sum(k)
!write(LUNOPRT,*) ' new dust mmr =', sum(mmr(:,k))
!end do

   end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    ! Cloud Ice & Snow
    !
    ! Cloud ice is maintained in advected species in CARMA, and we are only
    ! concerned with snow production by CARMA.
    !
    ! NOTE: To allow this code to be tested when not doing the cloud ice, but
    ! either doing nothing or doing detrainment, use the ice properties to convert
    ! from the 2 moment values to a size distribution.
    !
    ! NOTE: To keep mass and energy conservation happy, on the first step when
    ! camra_do_cldice is true, we take the bulk values and convert them
    ! into bins; however, this might cause issues with the CARMA growth code.

    !  JAS: carma_do_cldice is not declared within this module. Somehow
    ! carma_do_cldice exists via the Use carma_flags_mod statement.
    !  physics/cam/carma_flags_mod.F90 Sub read_carma_nl reads
    ! carma_do_cldice from the Carma namelist

    ! JAS:  I want to simplify the logic
    !       If Carma ice is diagnostic ( .not. carma_do_cldice ) then Diagnose the Bins
    !       If it's the first timestep and Carma ice is prognostic then Diagnose the Bins
    !      from the MG ice variables

    ! JAS: Precedence of logical operators: 1) .not. 2) .and. 3) .or.

    !print*
    !print*, 'carma_model_mod.F90 carma_do_cldice #1', carma_do_cldice
    !print*

    if ((.not. carma_do_cldice) .or. (is_first_step() .and. carma_do_initice)) then

      mmr = 0.0_f
      mmr1 = 0.0_f
      mmr2 = 0.0_f

      igroup = I_GRP_CRDICE
      ielem  = I_ELEM_CRDICE
      
      call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')

      ! Have to copy this, because size_dist_param_basic has an intent(inout)
      ! argument, but state is intent(in).
      limNumber = state%q(icol, :, ixnumice)

      ! Subroutine from MG utilities.
      call size_dist_param_basic(mg_ice_props, state%q(icol, :, ixcldice), limNumber, lambda(:))

      ! For ice, assume mu is 0.
      mu = 0._r8
  
      call CARMA_GetMmrFromGamma(carma, r(:), dr(:), rmass(:), &
           state%q(icol, :, ixcldice), limNumber, mu(:), &
           lambda(:), mmr1(:, :), rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetMmrFromGamma failed.')
    
      ! JAS:  MG 2.0 has Cloud Ice and Snow for prognostic variables
      !      Add the SNOWQM to the CLDICE
      limNumber = state%q(icol, : , ixnumsno)

      call size_dist_param_basic(mg_Snow_props, state%q(icol, :,ixsnowqm), limNumber, lambda(:))
      mu = 0._r8

      call CARMA_GetMmrFromGamma(carma, r(:), dr(:), rmass(:), state%q(icol, :, ixsnowqm), &
         limNumber, mu(:), lambda(:), mmr2(:, :), rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetMmrFromGamma failed.')

      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr1( ibin, :) + mmr2(ibin, :), rc)
        if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_SetBin failed.')
      end do
      
    else
    
      ! If CARMA is keeping track of ice, then total up the detrained and in-situ ice to
      ! make sure that no one else has messed with the ice fields. If changes were made,
      ! then adjustments need to be made to the totals to prevent mass and energy conservation
      ! errors within CAM. The difference could be accounted for in snow_str and prec_str.
      !
      ! NOTE: Advection, diffusion, ... may have affected the tracer values since the
      ! previous time step, however, the tracer correlations need to remain intact for
      ! CARMA to work properly. Also, no special processing can occur on the cldice fields
      ! outside of CARMA, since they are merely diagnostic fields of the CARMA state.
      if (carma_do_mass_check3 .or. carma_do_mass_fix) then
      
        call CARMA_GetTotalIceAndSnow(carma, cstate, iceMass, iceNumber, snowSurface, rc)
        if (rc < RC_OK) call endrun( 'CARMA_DiagnoseBins :: CARMA_GetTotalIceAndSnow failed.')
      
        do k = 1, pver

          carma_ice = iceMass(k) + snowMass(k)

          if (carma_ice <= qsmall ) then
            carma_ice = 0._r8
          end if
          
          if (carma_ice /= state%q(icol, k, ixcldice)) then

            if (carma_do_mass_check3) then
              if (abs(carma_ice - state%q(icol, k, ixcldice)) / max(abs(carma_ice), &
                  abs(state%q(icol, k, ixcldice))) >= 1e-10_r8)  then
                if (do_print) then
                  call CARMASTATE_Get(cstate, rc, lat=lat, lon=lon)
                  if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_Get failed.')
                  
                  write(LUNOPRT,1) icol, k, lat, lon, state%q(icol, k, ixcldice), &
                    carma_ice, (carma_ice - state%q(icol, k, ixcldice)) / max(abs(carma_ice), &
                    abs(state%q(icol, k, ixcldice)))
        
                  write(LUNOPRT,*) "  CAM cldice   :  ", state%q(icol, k, ixcldice)
                  write(LUNOPRT,*) ""
                  write(LUNOPRT,*) "  CARMA cldice :  ", iceMass(k) + snowMass(k)
                  write(LUNOPRT,*) "  CARMA ice    :  ", iceMass(k)
                  write(LUNOPRT,*) "  CARMA snow   :  ", snowMass(k)
                end if
              end if
            end if
          
            if (carma_do_mass_fix) then
         
              diff = ( state%q (icol, k, ixcldice) &
                       - ( iceMass(k) + snowMass(k) ) ) &
                     * ( state % pdel(icol, k) / gravit) / dt / 1000._r8
  
              snow_str(icol) = snow_str(icol) + diff
              prec_str(icol) = prec_str(icol) + diff
              
              if (carma_do_print_fix) then
                if (do_print) then
                   write(LUNOPRT,*) "  CARMA_DiagnoseBins::&
                        &WARNING - Adjusting prec_str for ice mass difference", &
                        icol, k, (state%q(icol, k, ixcldice) - (iceMass(k) + snowMass(k)))
                end if
              end if
            end if
          end if
        end do
      end if
    end if
    

    !
    ! Use the CAM mass and number (CLDLIQ and NUMLIQ) to determine an initial
    ! size distribution.

    ! JAS:  I would like water droplets to be Diagnosed from CldLiq and NumLiq
    !      for the initialization, similarly to ice in Carma Cirrus
    !       The 'not carma_do_cldliq' will Diagnose liquid water bins for
    !      Carma every timestep if Carma is not prognosing the liquid water.
    !      This diagnosis is required every timestep for Carma Cirrus.
    if ((is_first_step() .and. carma_do_initliq) .or. (.not. carma_do_cldliq)) then
      mmr = 0.0_f
      mmr1 = 0.0_f
      mmr2 = 0.0_f

      igroup = I_GRP_CRLIQ
      ielem  = I_ELEM_CRLIQ
    
      call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')

      ! Prevent size_dist_param from trying to change state by passing it a
      ! copy of the number concentration.
      limNumber = state%q(icol, :, ixnumliq)

      ! Subroutine from MG utilities.
      call size_dist_param_liq(mg_liq_props, state%q(icol, :, ixcldliq), &
           limNumber, rhoa_wet(:), mu(:), lambda(:))

      call CARMA_GetMmrFromGamma(carma, r(:), dr(:), rmass(:), state%q(icol, :, ixcldliq), &
           limNumber, mu(:), lambda(:), mmr1(:, :), rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetMmrFromGamma failed.')

      ! JAS:  MG 2.0 has Cloud Liquid and Rain for prognostic variables
      !      Add the RainQm to the CldLiq
      limNumber = state%q(icol, : , ixnumrai)

      call size_dist_param_liq(mg_Rain_props, state%q(icol, :, ixrainqm), &
           limNumber, rhoa_wet(:), mu(:), lambda(:))

      call CARMA_GetMmrFromGamma(carma, r(:), dr(:), rmass(:), state%q(icol, :, ixrainqm), &
           limNumber, mu(:), lambda(:), mmr2(:, :), rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetMmrFromGamma failed.')

      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr1(ibin, :) + mmr2(ibin, :), rc)
        if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_SetBin failed.')
      end do
    end if
           
    if (carma_do_mass_check2 .or. carma_do_mass_check3) then
    
      ! Check to see of the mass that we get back adds up.
      call CARMA_GetTotalWaterAndRain(carma, cstate, waterMass, waterNumber, rainSurface, rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetTotalWaterAndRain failed.')
      
      if (carma_do_mass_check3) then
        do k = 1, pver
                
          carma_water = waterMass(k)
          if (carma_water <= qsmall ) then
            carma_water = 0._r8
          end if
    
          ! The routine that provides the modal properties for water has a miniumum of 1e-18.
          ! This causes problems in comparisons, since smaller qc values are seen in the data,
          ! but CARMA's bins won't have values that small.
          if (carma_water /= state%q(icol, k, ixcldliq)) then

            if (abs(carma_water - state%q(icol, k, ixcldliq)) &
                 / max(abs(carma_water), abs(state%q(icol, k, ixcldliq))) &
                 >= 1e-10_r8) then

              if (do_print) then
                call CARMASTATE_Get(cstate, rc, lat=lat, lon=lon)
                if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_Get failed.')

                write(LUNOPRT,2) icol, k, lat, lon, state%q(icol, k, ixcldliq), &
                  carma_water, (carma_water - state%q(icol, k, ixcldliq)) / max(abs(carma_water), &
                  abs(state%q(icol, k, ixcldliq)))
      
                write(LUNOPRT,*) "  CAM cldliq   :  ", state%q(icol, k, ixcldliq)
                write(LUNOPRT,*) ""
                write(LUNOPRT,*) "  CARMA cldliq :  ", waterMass(k)
              end if
            end if
          end if
        end do
      end if
      
  
      ! Check for total water conservation by CARMA.
      if (carma_do_mass_check2) then
        call CARMA_GetTotalIceAndSnow (carma, cstate, iceMass, iceNumber, snowSurface, rc)
        if (rc < RC_OK) call endrun('CARMA_DiagnoseBins :: CARMA_GetTotalIceAndSnow failed.')

        call CARMA_CheckMassAndEnergy(carma, cstate, "Carma_DiagnoseBins", state, icol, dt, dlf, waterMass, rainSurface, iceMass, snowSurface, rc)    
        if (rc < RC_OK) call endrun('CARMA_DiagnoseBins :: CARMA_CheckMassAndEnergy failed.' )
      end if
    end if
    
    
    ! CAM determines the activation rate for liquid clouds. The rate is specified in
    ! the pbuf variable NPCCN (#/kg/s).
    !
    ! We will set this into a specific bin that is determined by the namelist
    ! parameter carma_dropact_bin.
    !
    ! NOTE: CAM also changes droplet number based upon changes in cloud fraction and
    ! regeneration. Here we only want the positive side of this number, and probably
    ! want regeneration turned off. If it is desired to have cloud mass reduced when
    ! activation is negative, then this code would need to be modified to divide
    ! up the negative rates amongst the bins and the code to remove the mass would
    ! to be included. This is currently in newstate_calc.
    if (carma_do_cldliq) then
      call pbuf_get_field(pbuf, pbuf_get_index("NPCCN"), activationrate)
      
      call CARMAGROUP_Get(carma, I_GRP_CRLIQ, rc, ienconc=ienconc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')

      call CARMASTATE_GetBin(cstate, ienconc, carma_dropact_bin, mmr, rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_GetBin failed.')

      call CARMASTATE_SetBin(cstate, ienconc, carma_dropact_bin, mmr, rc, activationrate=activationrate(icol,:))
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_SetBin failed.')
    endif    
    
    return
  end subroutine CARMA_DiagnoseBins
  
  
  !! For diagnostic groups, determines the tendencies on the CAM state from the CARMA bins.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen (CB) 
  !!
  !!  @version 2019 December 
  !!  @author  Jamison A. Smith (JAS)
  !!  
  !!  JAS:  MG 2.0 now has prognostic Rain and Snow, so some changes are
  !!       necessary

  !subroutine Carma_DiagnoseBulk ( carma, cstate, cam_out, state, pbuf, ptend &
  !                              , icol, dt, rc &
  !                              , rliq &
  !                              , prec_str &
  !                              , snow_str &
  !                              , prec_sed, snow_sed, tnd_qsnow, tnd_nsnow &
  !                              , re_ice)
  subroutine CARMA_DiagnoseBulk(carma, cstate, cam_out, state, pbuf, ptend, icol, dt, rc, dlf, prec_str, snow_str, &
    prec_sed, snow_sed, re_ice, so4_ac_frac, so4_co_frac, dst_ac_frac, dst_co_frac)
    ! (LW: modified for homo. & hetero. nuc., Sep 2022)

    use camsrfexch,     only: cam_out_t
    use time_manager,   only: is_first_step

    implicit none
    

    type (carma_type)     , intent (in)    :: carma   !! the carma object
    type (carmastate_type), intent (inout) :: cstate  !! the carma state object
    type (cam_out_t)      , intent (inout) :: cam_out !! cam output to surface models
    type (physics_state)  , intent (in)    :: state   !! physics state variables
    type (physics_buffer_desc), pointer    :: pbuf(:)    !! physics buffer
    type(physics_ptend)  , intent(inout)   :: ptend     !! constituent tendencies
    integer              , intent(in)      :: icol      !! column index
    real(r8)             , intent(in)      :: dt        !! time step
    integer              , intent(out)     :: rc        !! return code, negative indicates failure
    real(r8), intent(in), optional         :: dlf(pcols,pver) !! detrained condensation rate (kg/kg/s)
    real(r8), Intent(inout), optional      :: prec_str (pcols)  !! Total sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional      :: snow_str(pcols)   !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(inout), optional      :: prec_sed (pcols)  !! Total precip from cloud sedimentation (m/s)
    real(r8), intent(inout), optional      :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(out), optional        :: re_ice(pcols,pver)    !! ice effective radius (m) 

    !!---------------- Zhu Sep 12 2024 -------------------
    real(r8), pointer,dimension(:,:)      :: degrp_ptr ! graupel effective diameter (m)
    real(r8), pointer,dimension(:,:)      :: cldgrp_ptr !incloud graupel mixng ratio (kg/kg)
    !!----------------------------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Sep 2022): added for homo. & hetero. nuc.
    real(r8), intent(in), optional       :: so4_ac_frac(pcols,pver,NBIN)
    real(r8), intent(in), optional       :: so4_co_frac(pcols,pver,NBIN) 
    real(r8), intent(in), optional       :: dst_ac_frac(pcols,pver,NBIN)
    real(r8), intent(in), optional       :: dst_co_frac(pcols,pver,NBIN) 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! These values are chosen to match up with how small cloud ice values are handled in
    ! micro_mg.

    real(r8), parameter                  :: omsm       = 0.99999_r8    ! Prevents roundoff errors


    integer                              :: igroup    ! group index
    integer                              :: ielem     ! element index
    integer                              :: ibin      ! bin index
    integer                              :: icore     ! core index
    integer                              :: icorelem(NELEM) ! core indexes for group
    integer                              :: ncore     ! number of core elements
    
    real(kind=f)                         :: iceMass(pver)      ! ice mass mixing ratio (kg/kg)
    real(kind=f)                         :: iceNumber(pver)    ! ice number mixing ratio (#/kg)
    real(kind=f)                         :: snowMass(pver)     ! snow mass mixing ratio (kg/kg)
    real(kind=f)                         :: snowNumber(pver)   ! snow number (#/kg)
    real(kind=f)                         :: snowSurface        ! snow on surface (kg/m2)
    real(kind=f)                         :: waterMass(pver)    ! water mass mixing ratio (kg/kg)
    real(kind=f)                         :: waterNumber(pver)  ! water number mixing ratio (#/kg)
    !real(kind=f)                         :: cloudwaterMass(pver)    ! cloud water mass mixing ratio (kg/kg)
    !real(kind=f)                         :: cloudwaterNumber(pver)  ! cloud water number mixing ratio (#/kg)
    real(kind=f)                         :: rainSurface        ! rain on surface (kg/m2)
    real(kind=f) :: iceRe(pver)        ! ice effective radius (m)
    real(kind=f) :: LiqRe(pver)        ! JAS 2020Feb25 liq eff rad (m)

    real(r8)     :: newRain            ! [Total] sfc flux of rain from stratiform (m/s)
    real(r8)     :: newSnow            ! [Total] sfc flux of snow from stratiform (m/s)

    real(kind=f)                         :: mmr(pver)          ! mass mixing ratio (#/kg)
    real(kind=f)                         :: mmrcore(pver)      ! core mass mixing ratio (#/kg)
    real(kind=f)                         :: nmr(pver)          ! number mixing ratio (#/kg)
    real(kind=f)                         :: r(NBIN)            ! radius (cm)
    real(kind=f)                         :: sfc                ! surface mass (kg/m2)
    real(kind=f)                         :: sfccore            ! core surface mass (kg/m2)
    
    !!---------------- Zhu Sep 12 2024 --------------
    real(kind=f)                         :: graupelicmass(pver), graupelnum(pver), graupelde(pver)
    !!-----------------------------------------------
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! LW (Sep 2022): added for homo & hetero. nuc.
! MAM 4 aerosols
!    integer :: m, nn, lptr_ac, lptr_at, lptr_co, nptr_ac, nptr_at, nptr_co, k
!    integer ::  lptr_ac, lptr_at, lptr_co, nptr_ac, nptr_at, nptr_co, k
    integer ::  lptr, nptr, k, d_loc_so4, d_loc_dst

!    character(len=32)     :: str32
    character(len=*), parameter :: routine = 'CARMA_DiagnoseBulk'
    integer               :: nspec                            ! number of chemical species
    integer               :: nmodes                           ! number of aerosol modes

    real(kind=f)                         :: so4_mmr_co_new(pcols,pver)  ! updated sulfate mmr, coarse mode (kg/kg)
    real(kind=f)                         :: so4_nmr_co_new(pcols,pver)  ! updated sulfate number mixing ratio, coarse mode (#/kg)

    real(kind=f)                         :: r_so4(NBIN)            ! radius (cm)
    real(kind=f)                         :: r_dst(NBIN)            ! radius (cm)
    real(kind=f)                         :: dr(NBIN)
    real(kind=f)                         :: rmass(NBIN)
    real(kind=f)                         :: mmr_sulf(NBIN, pver)
    real(kind=f)                         :: mmr_dust(NBIN, pver)
    real(kind=f)                         :: nmr_update_ac(pcols,pver)
    real(kind=f)                         :: nmr_update_co(pcols,pver)
    real(kind=f)                         :: max_d_so4
    real(kind=f)                         :: max_d_dst
    real(kind=f)                         :: flag
!    real(kind=f)                         :: so4_mmr_sum(pver)
!    real(kind=f)                         :: dust_mmr_sum(pver)
    real(kind=f)                         :: so4_mmr_sum_ac(pver)
    real(kind=f)                         :: so4_mmr_sum_co(pver)
    real(kind=f)                         :: dust_mmr_sum_ac(pver)
    real(kind=f)                         :: dust_mmr_sum_co(pver)
    real(kind=f)                         :: tmp_so4(NBIN)
    real(kind=f)                         :: tmp_dust(NBIN)


    real(kind=f)                         :: tmass_carma
    real(kind=f)                         :: tmass_state_ac
    real(kind=f)                         :: tmass_state_co
    real(kind=f)                         :: mmr_tot_old
    real(kind=f)                         :: mmr_tot_new
    real(kind=f)                         :: mmr_tot_old_co
    real(kind=f)                         :: mmr_tot_new_co
    real(kind=f)                         :: mmr_tot_old_ac
    real(kind=f)                         :: mmr_tot_new_ac
    real(kind=f)                         :: dst_so4_wght_old
    real(kind=f)                         :: dst_so4_wght_new
    real(kind=f)                         :: dst_so4_nmr_old
    real(kind=f)                         :: dst_so4_nmr_new
    real(kind=f)                         :: dst_wght_old
    real(kind=f)                         :: dst_wght_new
    real(kind=f)                         :: so4_wght_old
    real(kind=f)                         :: so4_wght_new
    real(kind=f)                         :: dst_nmr_old
    real(kind=f)                         :: dst_nmr_new
    real(kind=f)                         :: so4_nmr_old
    real(kind=f)                         :: so4_nmr_new
!    real(kind=f)                         :: state_mmr_sum_so4(pver)
!    real(kind=f)                         :: state_mmr_sum_dust(pver)

    integer                              :: LUNOPRT
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Default return code.

    rc = RC_OK
    
    ! Aerosols
    !
    ! Currently, we are just using a fixed aerosol size distribution, but in the
    ! future this could be linked to the modal aerosols.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Sep 2022): added for homo. & hetero. nuc.
    ! Linked to the modal aerosol model

    ! Determine the changes to sulfate mass and number due to nucleation
    ! of ice
    call CARMA_GetTotalIceAndSnow(carma, cstate, iceMass, iceNumber, snowSurface, rc, iceRe=iceRe)    

    ! Calculate the tendencies of sulfate mass and number
    if (carma_do_bulk_tend) then


      if (carma_sulfate_method == "modal") then 
        ! Get aerosol group and element information
        ! Get dgnum for determining which mode to put the mmr and nmr back into
        call rad_cnst_get_mode_props(0, mode_coarse_idx, dgnumlo=dgnumlo_co, dgnumhi=dgnumhi_co)
        call rad_cnst_get_mode_props(0, mode_aitken_idx, dgnumlo=dgnumlo_at, dgnumhi=dgnumhi_at)
        call rad_cnst_get_mode_props(0, mode_accum_idx, dgnumlo=dgnumlo_ac, dgnumhi=dgnumhi_ac)

        ! Sulfates
        call CARMAGROUP_Get(carma, I_GRP_CRCN, rc, r=r_so4)
        if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')

        ! NOTE: Should only really being doing this if the sulfate method is using
        ! MAM.
        do ibin = 1, NBIN

          call CARMASTATE_GetBin(cstate, I_ELEM_CRCN, ibin, mmr, rc)
          if (rc < RC_OK) call endrun('CARMA_DiagnoseBulk::CARMASTATE_GetBin failed.')

          mmr_sulf(ibin,:) = mmr(:)
        end do

        ! Dust
        call CARMAGROUP_Get(carma, I_GRP_CRDUST, rc, r=r_dst)
        if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')

 
        do ibin = 1, NBIN

          call CARMASTATE_GetBin(cstate, I_ELEM_CRDUST, ibin, mmr, rc)
          if (rc < RC_OK) call endrun('CARMA_DiagnoseBulk::CARMASTATE_GetBin failed.')

          mmr_dust(ibin,:) = mmr(:)
        end do

        ! Get the mmr and nmr for each mode
        ! Coarse mode
        call rad_cnst_get_mode_num(0, mode_coarse_idx, 'a', state, pbuf, mode_nmr_co)
        call rad_cnst_get_aer_mmr(0, mode_coarse_idx, so4_coarse_idx, 'a', state, pbuf, so4_mmr_co)
        call rad_cnst_get_aer_mmr(0, mode_coarse_idx, dust_coarse_idx, 'a', state, pbuf, dust_mmr_co)
        call rad_cnst_get_aer_mmr(0, mode_coarse_idx, sslt_coarse_idx, 'a', state, pbuf, sslt_mmr_co)

        ! Accumulation mode
        call rad_cnst_get_mode_num(0, mode_accum_idx, 'a', state, pbuf, mode_nmr_ac)
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, so4_accum_idx, 'a', state, pbuf, so4_mmr_ac)
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, dust_accum_idx, 'a', state, pbuf, dust_mmr_ac)
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, sslt_accum_idx, 'a', state, pbuf, sslt_mmr_ac)
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, pom_accum_idx, 'a', state, pbuf, pom_mmr_ac)
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, soa_accum_idx, 'a', state, pbuf, soa_mmr_ac)
        call rad_cnst_get_aer_mmr(0, mode_accum_idx, bc_accum_idx, 'a', state, pbuf, bc_mmr_ac)

        ! Update the mass and number mixing ratios...this is the tricky part
        do k = 1, pver

          ! Note that mmr coming out of CARMA was originally the sum of the accumulation and coarse modes for their 
          ! respective aerosol species. Therefore if no nucleation happens so4_mmr_sum will equal (so4_mmr_ac + so4_mmr_co).
          ! As as a result, extra mass will be sent back to MAM. Stop this from happening.
!write(LUNOPRT,*) ' So4 carma pre adjust =', so4_mmr_sum(k)
         

          ! Setup the mass so that it can be returned to MAM3
          if (is_first_step()) then
            so4_mmr_sum_ac(k) = so4_mmr_ac(icol,k)
            so4_mmr_sum_co(k) = so4_mmr_co(icol,k)
            dust_mmr_sum_ac(k) = dust_mmr_ac(icol,k)
            dust_mmr_sum_co(k) = dust_mmr_co(icol,k)
          else
            tmp_so4 = mmr_sulf(:,k) * so4_ac_frac(icol,k,:)
            so4_mmr_sum_ac(k) = sum(tmp_so4)
            tmp_so4 = mmr_sulf(:,k) * so4_co_frac(icol,k,:)
            so4_mmr_sum_co(k) = sum(tmp_so4)
            tmp_dust = mmr_dust(:,k) * dst_ac_frac(icol,k,:)
            dust_mmr_sum_ac(k) = sum(tmp_dust)
            tmp_dust = mmr_dust(:,k) * dst_co_frac(icol,k,:)
            dust_mmr_sum_co(k) = sum(tmp_dust)
          end if

          ! Determine the fraction of the mode nmr is due to a given species (assuming external mixing)
          mmr_tot_old = so4_mmr_ac(icol,k) + dust_mmr_ac(icol,k) + sslt_mmr_ac(icol,k) + pom_mmr_ac(icol,k) + soa_mmr_ac(icol,k) + bc_mmr_ac(icol,k)
  !        mmr_tot_new = so4_mmr_sum(k) + dust_mmr_sum(k) + sslt_mmr_ac(icol,k) + pom_mmr_ac(icol,k) + soa_mmr_ac(icol,k) + bc_mmr_ac(icol,k)
          mmr_tot_new = so4_mmr_sum_ac(k) + dust_mmr_sum_ac(k) + sslt_mmr_ac(icol,k) + pom_mmr_ac(icol,k) + soa_mmr_ac(icol,k) + bc_mmr_ac(icol,k)

          dst_so4_wght_old = so4_mmr_ac(icol,k)
          if (mmr_tot_old .gt. 0._f) then
            dst_so4_wght_old = (dst_so4_wght_old + dust_mmr_ac(icol,k)) / mmr_tot_old
          end if

!         dst_so4_wght_new = (so4_mmr_sum(k) + dust_mmr_sum(k)) / mmr_tot_new
          dst_so4_wght_new = so4_mmr_sum_ac(k)
          if (mmr_tot_new .gt. 0._f) then
            dst_so4_wght_new = (dst_so4_wght_new  + dust_mmr_sum_ac(k)) / mmr_tot_new
          end if

          dst_so4_nmr_old = mode_nmr_ac(icol,k) * dst_so4_wght_old
          dst_so4_nmr_new = mode_nmr_ac(icol,k) * dst_so4_wght_new

          nmr_update_ac(icol,k) = mode_nmr_ac(icol,k) - dst_so4_nmr_old + dst_so4_nmr_new
!         nmr_update_co(icol,k) = state%q(icol, k, ixnum_a3)

          ptend%q(icol, k, ixso4_a1) = (so4_mmr_sum_ac(k) - state%q(icol, k, ixso4_a1)) / dt
          ptend%q(icol, k, ixdust_a1) = (dust_mmr_sum_ac(k) - state%q(icol, k, ixdust_a1)) / dt

          ! Determine the fraction of the mode nmr is due to a given species (assuming external mixing)
          mmr_tot_old = so4_mmr_co(icol,k) + dust_mmr_co(icol,k) + sslt_mmr_co(icol,k)
          mmr_tot_new = so4_mmr_sum_co(k) + dust_mmr_sum_co(k) + sslt_mmr_co(icol,k)

          dst_so4_wght_old = (so4_mmr_co(icol,k) + dust_mmr_co(icol,k)) / mmr_tot_old
          dst_so4_wght_new = (so4_mmr_sum_co(k) + dust_mmr_sum_co(k)) / mmr_tot_new

          dst_so4_nmr_old = mode_nmr_co(icol,k) * dst_so4_wght_old
          dst_so4_nmr_new = mode_nmr_co(icol,k) * dst_so4_wght_new

          nmr_update_co(icol,k) = mode_nmr_co(icol,k) - dst_so4_nmr_old + dst_so4_nmr_new

          ptend%q(icol, k, ixso4_a3)  = (so4_mmr_sum_co(k) - state%q(icol, k, ixso4_a3)) / dt
          ptend%q(icol, k, ixdust_a3) = (dust_mmr_sum_co(k) - state%q(icol, k, ixdust_a3)) / dt
        end do   ! end levels loop

        ! Calculate the nmr tendency
        ptend%q(icol, :, ixnum_a1) = (nmr_update_ac(icol,:) - state%q(icol, :, ixnum_a1)) / dt
        ptend%q(icol, :, ixnum_a3) = (nmr_update_co(icol,:) - state%q(icol, :, ixnum_a3)) / dt
      end if    ! End updating the aerosols
    end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Cloud Ice & Snow
    !
    ! Determine the changes to cloud ice (mass and number) and snow (mass and number)
    ! by looking at the totals of the detrained and in-situ ice.
    
    ! Get the total ice.
    !+++ Noted by Cheng: I think now we missed iceMass and iceNumber from Graupel, due to new added codes
    !+++ We can fix it later by adding Graupel num/mass into CARMA_GetTotalIceAndSnow,
    !+++ or load and add from CARMA_GetTotalGraupel (but that is just for radiation codes) 
    call CARMA_GetTotalIceAndSnow(carma, cstate, iceMass, iceNumber, snowSurface, rc, iceRe=iceRe)
    
    ! Calculate the tendencies on CLDICE, NUMICE, QSNOW and NSNOW
    if (carma_do_bulk_tend) then
    
      ptend%q(icol, :, ixcldice) = (iceMass(:) - state%q(icol, :, ixcldice)) / dt
      ptend%q(icol, :, ixnumice) = (iceNumber(:) - state%q(icol, :, ixnumice)) / dt

      where (iceMass(:) < qsmall)
        ptend%q(icol, :, ixcldice) = - state%q(icol, :, ixcldice) / dt
        ptend%q(icol, :, ixnumice) = - state%q(icol, :, ixnumice) / dt
      end where

      ! JAS:  Snow is a constituent now
      !       1) Convert CldIce, NumIce, SnowQm, and NumSno into CrDIceXX prior
      !          to the first timestep ( XX = bins of Carma Ice, e.g. 01 - 28 )
      !       2) Let Carma process CrDIce
      !       3) Convert CrDIceXX into CldIce and NumIce
      !          a) CrDIceXX are the prognostic variables that are the actual
      !             state of the model ice
      !          b) CldIce may be used by Clubb in its generation of turbulence
      !          c) CrDIceXX will be transported by Clubb's turbulence
      !          d) CldIce will be used by Rrtmg for heating rates
      !          d) CrDiceXX will be transported by Finite Volume dynamics
      !       4) Zero the Snow Variables SnowQm and NumSno by applying
      !          tend = -state ( q ) / dt.  Where did the snow go?  See
      !          3a) CrDIceXX are the actual state of the model ice, so
      !          Snow is made of the large bins in this distribution.

      ! JAS:  If MG Snow Vars have been converted into Carma DIce Vars then
      !      zero the MG Snow Vars
      !       CARMA_DiagnoseBulk will have to decide how to repopulate the
      !      MG Snow Vars when that action is desired
      if (is_first_step() .and. carma_do_cldice) then
        ptend%q(icol, :, ixsnowqm) = -state%q(icol, :, ixsnowqm) / dt
        ptend%q(icol, :, ixnumsno) = -state%q(icol, :, ixnumsno) / dt
      endif

      ! Now we need to change the reserve liquid. This was indicating the amount of
      ! water than was not in the atmosphere because it was in convection (dlf). Now
      ! we have included that water, but we have removed water representing snow in
      ! the atmosphere. This needs to be communicated to the CAM microphysics which
      ! will take care of actually precipitating or evaporating the snow.
      !
      ! NOTE: Putting CARMA in a microphysics loop means that you can't modify
      ! rliq, since the rates need to reamin constant for the timestep.
!      if (present(rliq)) rliq(icol) = rliq(icol) + sum(snowMass(:) * (state%pdel(icol, :) / gravit)) / dt / 1000._r8
      
      ! The ice effective radius is used by the radiation code; however, it uses a mass
      ! weighted effective diameter in um.
      if (present(re_ice)) re_ice(icol, :) = iceRe

    end if

          
    ! Water Drops
    !
    ! Calculate the total mass and total number of the water drops, and then
    ! determine the appropriate tendencies.
    ! call CARMA_GetTotalWaterAndRain (carma, cstate, waterMass, waterNumber, rainSurface, rc, liqRe=liqRe )
    call CARMA_GetCloudWaterAndRain (carma, cstate, CloudWaterMass=waterMass, CloudWaterNumber=waterNumber, rainSurface=rainSurface, rc=rc, liqRe=liqRe )

    ! Calculate the tendencies on CLDLIQ and NUMLIQ
    if (carma_do_bulk_tend) then
    
      ! In CAM in cldwat2m, a couple of things are done:
      !
      !   1) If cldliq < qsmall, then the number desnity is set to 0.
      !   2) to keep from overshooting into negative values, they don't try to drive
      !      the value all the way to 0.
      ptend%q(icol, :, ixcldliq) = (waterMass(:) - state%q(icol, :, ixcldliq)) / dt
      ptend%q(icol, :, ixnumliq) = (waterNumber(:) - state%q(icol, :, ixnumliq)) / dt

      where(waterMass(:) < qsmall)
        ptend%q(icol, :, ixcldliq) = -state%q(icol, :, ixcldliq ) / dt
        ptend%q(icol, :, ixnumliq) = -state%q(icol, :, ixnumliq ) / dt
      end where

      ! JAS:  Rain is a constituent now
      !       1) Convert CldLiq, NumLiq, RainQm, and NumRai into CrLiqXX prior
      !          to the first timestep ( XX = bins of Carma Liq, e.g. 01 - 28 )
      !       2) Let Carma process CrLiq
      !       3) Convert CrLiqXX into CldLiq and NumLiq
      !          a) CrLiqXX are the prognostic variables that are the actual
      !             state of the model liquid
      !          b) CldLiq is used by Clubb in its generation of turbulence
      !          c) CrLiqXX will be transported by Clubb's turbulence
      !          d) CldLiq will be used by Rrtmg for heating rates
      !          d) CrLiqXX will be transported by Finite Volume dynamics
      !       4) Zero the Rain Variables RainQm and NumRai by applying
      !          tend = -state ( q ) / dt.  Where did the rain go?  See
      !          3a) CrLiqXX are the actual state of the model liquid, so
      !          Rain is made of the large bins in this distribution.

      ! JAS:  If MG Rain Vars have been converted into Carma Liq Vars then
      !      zero the MG Rain Vars

      if ( .not. carma_do_cldliq &
           .or. is_first_step() .and. carma_do_cldliq ) then   
        ptend%q(icol, :, ixrainqm) = -state%q(icol, :, ixrainqm) / dt
        ptend%q(icol, :, ixnumrai) = -state%q(icol, :, ixnumrai) / dt
      endif

      ! JAS:  Carma_DiagnoseBulk will have to decide how to repopulate the
      !      MG Rain Vars when that action is desired

    end if

    !! -------------------- Zhu Sep 12 2024 -------------------------
    ! get the total graupel for radiation code
    call CARMA_GetTotalGraupel (carma, cstate, graupelicmass, graupelnum, graupelde, rc)

    if (carma_do_cldice .and. carma_do_cldliq ) then                                                  
        call pbuf_get_field ( pbuf, pbuf_get_index("ICGRPCARMA"),cldgrp_ptr)                          
        call pbuf_get_field ( pbuf, pbuf_get_index("DEGCARMA"),degrp_ptr)

        cldgrp_ptr(icol, :) = graupelicmass
        degrp_ptr(icol, :) = graupelde
      
    endif
    !! ------------------------------------------------------------
        
    ! For mass balance, we also need to supply the total precipation and snow. Not
    ! all of the snow may make the ground, but that will be determined later in the
    ! MG microphysics. For now, we need to account for all condensate that is not
    ! in CLDICE or CLDLIQ.
    !
    ! Need the 1000. to convert from kg/m2/s to m/s
    newSnow = snowSurface
    newRain = rainSurface

    if (present(snow_sed)) snow_sed(icol) = snow_sed(icol) + newSnow / dt / 1000._r8
    if (present(prec_sed)) prec_sed(icol) = prec_sed(icol) + (newRain + newSnow) / dt / 1000._r8
    
    if (present(snow_str)) snow_str(icol) = snow_str(icol) + newSnow / dt / 1000._r8
    if (present(prec_str)) prec_str(icol) = prec_str(icol) + (newRain + newSnow) / dt / 1000._r8

    ! Check for total water conservation by CARMA.
    if (carma_do_mass_check) then
      
      ! The CAM state has not been updated yet, so compare the original CAM state
      ! with the new CARMA state.
      !
      ! NOTE: This needs to be fixed for rliq, prec_str, and snow_str not being
      ! passed in.
      call CARMA_CheckMassAndEnergy(carma, cstate, "CARMA_DiagnoseBulk", state, &
       icol, dt, dlf, waterMass, rainSurface, iceMass, snowSurface, rc)
      !+++ Check if we should use waterMass (total) or cloudwaterMass here for 0-50 µm.
    end if

    return
  end subroutine CARMA_DiagnoseBulk

  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in Carma_init.
  !!
  !! @author  Chuck Bardeen (CB)
  !! @version May-2009
  !! @author  Jamison A Smith (JAS)
  !! @version 2020 Jun

  subroutine CARMA_InitializeModel (carma, lq_carma, rc)
    use cam_history,      only: addfld,  horiz_only, add_default
    use constituents,     only: cnst_get_ind, pcnst

    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer                            :: ibin                  ! bin index
    integer                            :: i
    integer                            :: itemp                 ! temperature index
    integer                            :: LUNOPRT
    
    logical                            :: do_print_init
    logical                            :: do_grow
    logical                            :: do_detrain
    logical                            :: do_thermo

    real(kind=f)                       :: r(NBIN)               ! bin center radius (cm)
    real(kind=f)                       :: dr(NBIN)              ! bin width (cm)
    real(kind=f)                       :: rmass(NBIN)           ! bin mass (g)
    real(kind=f)                       :: sub_d                 ! integration substep diameter (um)
    real(kind=f)                       :: sub_dd                ! integration substep width (um)
    real(kind=f)                       :: eshape                ! particle aspect ratio (> 0 is prolate)
    real(kind=f)                       :: shapeFactor           ! shape factor for maximum radius
    real(kind=f)                       :: lambda                ! fit factor for H&S 2010 size distribution
    real(kind=f)                       :: temp                  ! temperature (C)
     
    ! Default return code.

    rc = 0
    
    call CARMA_Get(carma, rc, do_print_init=do_print_init, LUNOPRT=LUNOPRT, &
         do_grow=do_grow, do_detrain=do_detrain, do_thermo=do_thermo)
    if (rc < RC_OK) call endrun('CARMA_CheckMassAndEnergy::CARMA_Get failed.') 
    
    ! Lookup indices to other constituents that are needed.

    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('NUMICE', ixnumice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('NUMLIQ', ixnumliq)

    ! JAS: MG 2.0 has prognostic Rain and Snow

    call cnst_get_ind('RAINQM', ixrainqm)
    call cnst_get_ind('NUMRAI', ixnumrai)
    call cnst_get_ind('SNOWQM', ixsnowqm)
    call cnst_get_ind('NUMSNO', ixnumsno)
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! LW (Sep 2022): added for hetero. nuc.
    ! MAM 4 aerosols
    call cnst_get_ind('so4_a1', ixso4_a1)
    call cnst_get_ind('so4_a3', ixso4_a3)
    call cnst_get_ind('num_a1', ixnum_a1)
    call cnst_get_ind('num_a3', ixnum_a3)
    call cnst_get_ind('dst_a1', ixdust_a1)
    call cnst_get_ind('dst_a3', ixdust_a3)
    
    lq_carma(ixso4_a1) = .true.
    lq_carma(ixso4_a3) = .true.
    lq_carma(ixdust_a1) = .true.
    lq_carma(ixdust_a3) = .true.
    lq_carma(ixnum_a1) = .true.
    lq_carma(ixnum_a3) = .true.    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
    ! Add the CAM ice and liquid fields as some that could be modified by CARMA.
    lq_carma(ixcldice) = .true.
    lq_carma(ixnumice) = .true.
    lq_carma(ixcldliq) = .true.
    lq_carma(ixnumliq) = .true.      

    lq_carma(ixrainqm) = .true.
    lq_carma(ixnumrai) = .true.
    lq_carma(ixsnowqm) = .true.
    lq_carma(ixnumsno) = .true.      

    if (do_print_init) then
      write(LUNOPRT,*) ""
      write(LUNOPRT,*) "Initializing CARMA Detrainment"
      write(LUNOPRT,*) ""
      write(LUNOPRT,*) "  Using ice method = ", carma_dice_method
    end if
    
    ! For detrainment of ice, setup the fractions of ice that go into each bin and
    ! into snow. This can be done different ways:
    !
    !   - monodisperse
    !   - temperature dependent size distribution
    !
    ! In any of these, a fraction of the ice can go directly to snow, rather than
    ! going into bins first.
    !
    ! Puts all of the detraining cloud water from convection into the large scale cloud,
    ! and puts detraining cloud water into liquid and ice based on temperature partition
    call CARMAGROUP_Get ( carma, I_Grp_CrDIce, rc, r=r(:), dr=dr(:), eshape=eshape, rmass=rmass(:))
    if (rc < rc_ok) call endrun ('Carma_InitializeModel :: CarmaGroup_Get failed.')

    ! This size distribution in based upon the maximum diameter, so the ice particles have
    ! a shape then pass the largest dimension to the size distribution.
    !
    ! NOTE: This is assuming the shape is a spheroid. Should consider passing shape
    ! parameters out of setupvfall, so that f1 is available for this.
    if (eshape >= 1._f) then
      shapeFactor = eshape**(1._f / 3._f)
    else
      shapeFactor = eshape**(- 1._f / 3._f)
    end if
        
    dice_bin_fraction(:, :) = 0._f
    
    ! Heymsfield & Schmitt [2010] tmeperature dependent distribution
    if (carma_dice_method == "dist_hym2010") then
    
      ! Integrate over the defined temperaure range.
      do itemp = 1, NDTEMP
      
        temp   = 1._f - itemp       

        ! Determine the exponentianal factor of the number distribution from H&S eq. 7.\
        lambda = dist_hym2010_alpha * exp(temp * dist_hym2010_beta)

        ! Determine a mass distribution using from a size distribution using this
        ! lambda. The number distribution is N = N0 * exp(-lambda * D), with D in
        ! cm from H&S eq. 1. Since this is just to generate a PDF, just use an N0
        ! of 1.
        !
        ! NOTE: This mass distribution (dMdD) is based on the diameter in cm.
        do ibin = 1, NBIN
      
          ! Determine the fraction in each bin.
          !
          ! NOTE: The bins are wide realtive to this function, so sum over an interval

          sub_dd = 2._f * dr(ibin) * shapeFactor
          sub_d =  2._f * r(ibin)  * shapeFactor

          dice_bin_fraction(ibin, itemp) = dice_bin_fraction(ibin, itemp) + &
            rmass(ibin) / lambda * &
            (exp(-lambda * (sub_d - (sub_dd / 2._f))) - exp(-lambda * (sub_d + (sub_dd / 2._f))))
        end do
 
        ! The sum of the integral may not be exactly 1, so scale the total so as not to skew
        ! the amount going straight to snow.
        dice_bin_fraction(:, itemp) = dice_bin_fraction(:, itemp) / (sum(dice_bin_fraction(:, itemp)))
      end do
    
    ! Default to monodisperse
    else
    
      do ibin = 1, NBIN
        if (r(ibin) >= r_dice_mono) then
          dice_bin_fraction(ibin, :) = 1._f
          
          exit
        end if
      end do
    end if
    
        
    if (do_print_init) then
      do itemp = 1, NDTEMP, 10
      
        if ((itemp == 1) .or. (carma_dice_method == "dist_hym2010")) then
        
          if (carma_dice_method == "dist_hym2010") then
            write(LUNOPRT,*) ""
            write(LUNOPRT,*) "  Temperature = ", 1 - itemp, " C"
            write(LUNOPRT,*) ""
          end if
        
          write(LUNOPRT,*) ""
          write(LUNOPRT,*) "        ibin       r (um)                   fraction"

          do ibin = 1, NBIN
            write(LUNOPRT,*) ibin, r(ibin)*1e4_f, dice_bin_fraction(ibin, itemp)
          end do
        end if
      end do
    end if
    
    ! Log a warning message if doing growth or detrainment and not doing
    ! thermodynamics. This will cause an energy error to be reported by CAM.

    if ((do_grow .or. do_detrain) .and. .not. do_thermo) then
      if (do_print_init) then
         write(LUNOPRT,*) "CARMA_InitializeModel:&
              &WARNING - do_grow and/or do_detrain are selected without &
              &do_thermo which may result in energy conservation errors."
      end if
    end if

    
    ! Provide diagnostics for SO4 tendencies from other physics packages
    !
    ! NOTE: This can be useful for determining an SO4 budget and for debugging
    ! SO4 conservation.
    if (carma_do_budget_diags) then
    
      call addfld("QBD", horiz_only, 'A', 'kg/m2/s', 'Q burden')
      if (carma_diags_file > 0) call add_default("QBD", carma_diags_file, ' ')

      call addfld("QSF", horiz_only, 'A', 'kg/m2/s', 'Q surface flux')
      if (carma_diags_file > 0) call add_default("QSF", carma_diags_file, ' ')

      call addfld("CRDICEBD", horiz_only, 'A', 'kg/m2', 'Detrained ice burden')
      if (carma_diags_file > 0) call add_default("CRDICEBD", carma_diags_file, ' ')
      call addfld("CRSICEBD", horiz_only, 'A', 'kg/m2', 'In situ ice burden')
      if (carma_diags_file > 0) call add_default("CRSICEBD", carma_diags_file, ' ')
      call addfld("CRGRPBD", horiz_only, 'A', 'kg/m2', 'Graupel burden')
      if (carma_diags_file > 0) call add_default("CRGRPBD", carma_diags_file, ' ')
      call addfld("CRLIQBD", horiz_only, 'A', 'kg/m2', 'Liquid burden')
      if (carma_diags_file > 0) call add_default("CRLIQBD", carma_diags_file, ' ')

      call addfld("CRSCOREBD", horiz_only, 'A', 'kg/m2', 'Sulfate vore burden')
      if (carma_diags_file > 0) call add_default("CRSCOREBD", carma_diags_file, ' ')
      call addfld("CRDCOREBD", horiz_only, 'A', 'kg/m2', 'Dust core burden')
      if (carma_diags_file > 0) call add_default("CRDCOREBD", carma_diags_file, ' ')

      call addfld("CLDLIQBD", horiz_only, 'A', 'kg/m2', 'Cloud liquid burden')
      if (carma_diags_file > 0) call add_default("CLDLIQBD", carma_diags_file, ' ')
      call addfld("CLDICEBD", horiz_only, 'A', 'kg/m2', 'Cloud ice burden')
      if (carma_diags_file > 0) call add_default("CLDICEBD", carma_diags_file, ' ')
      call addfld("RAINQMBD", horiz_only, 'A', 'kg/m2', 'Rain burden')
      if (carma_diags_file > 0) call add_default("RAINQMBD", carma_diags_file, ' ')
      call addfld("SNOWQMBD", horiz_only, 'A', 'kg/m2', 'Snow burden')
      if (carma_diags_file > 0) call add_default("SNOWQMBD", carma_diags_file, ' ')

      call addfld("CRDICENBD", horiz_only, 'A', '#/m2', 'Detrained ice number burden')
      if (carma_diags_file > 0) call add_default("CRDICENBD", carma_diags_file, ' ')
      call addfld("CRSICENBD", horiz_only, 'A', '#/m2', 'In situ ice number burden')
      if (carma_diags_file > 0) call add_default("CRSICENBD", carma_diags_file, ' ')
      call addfld("CRGRPNBD", horiz_only, 'A', '#/m2', 'Graupel number burden')
      if (carma_diags_file > 0) call add_default("CRGRPNBD", carma_diags_file, ' ')
      call addfld("CRLIQNBD", horiz_only, 'A', '#/m2', 'Liquid number burden')
      if (carma_diags_file > 0) call add_default("CRLIQNBD", carma_diags_file, ' ')

      call addfld("NUMLIQBD", horiz_only, 'A', '#/m2', 'Cloud liquid number burden')
      if (carma_diags_file > 0) call add_default("NUMLIQBD", carma_diags_file, ' ')
      call addfld("NUMICEBD", horiz_only, 'A', '#/m2', 'Cloud ice number burden')
      if (carma_diags_file > 0) call add_default("NUMICEBD", carma_diags_file, ' ')
      call addfld("NUMRAIBD", horiz_only, 'A', '#/m2', 'Rain number burden')
      if (carma_diags_file > 0) call add_default("NUMRAIBD", carma_diags_file, ' ')
      call addfld("NUMSNOBD", horiz_only, 'A', '#/m2', 'Snow number burden')
      if (carma_diags_file > 0) call add_default("NUMSNOBD", carma_diags_file, ' ')

      call addfld("CRLIQRESBD", horiz_only, 'A', 'um', 'Liquid Effective Radius, small')
      if (carma_diags_file > 0) call add_default("CRLIQRESBD", carma_diags_file, ' ')
      call addfld("CRLIQREMBD", horiz_only, 'A', 'um', 'Liquid Effective Radius, medium')
      if (carma_diags_file > 0) call add_default("CRLIQREMBD", carma_diags_file, ' ')
      call addfld("CRLIQRELBD", horiz_only, 'A', 'um', 'Liquid Effective Radius, large')
      if (carma_diags_file > 0) call add_default("CRLIQRELBD", carma_diags_file, ' ')

      call addfld("CRICERESBD", horiz_only, 'A', 'um', 'Ice Effective Radius, small')
      if (carma_diags_file > 0) call add_default("CRICERESBD", carma_diags_file, ' ')
      call addfld("CRICEREMBD", horiz_only, 'A', 'um', 'Ice Effective Radius, medium')
      if (carma_diags_file > 0) call add_default("CRICEREMBD", carma_diags_file, ' ')
      call addfld("CRICERELBD", horiz_only, 'A', 'um', 'Ice Effective Radius, large')
      if (carma_diags_file > 0) call add_default("CRICERELBD", carma_diags_file, ' ')


    end if


    if (carma_do_package_diags) then
      
      ! Iterate of the packages that have be instrumented. These should match the calls
      ! in physpkg.f90.
      do i = 1, carma_ndiagpkgs

        call addfld("QTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Q total tendency')
        if (carma_diags_file > 0) call add_default("QTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("CRDICETC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Detrained ice tendency')
        if (carma_diags_file > 0) call add_default("CRDICETC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRSICETC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//',  In situ ice tendency')
        if (carma_diags_file > 0) call add_default("CRSICETC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRGRPTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Graupel tendency')
        if (carma_diags_file > 0) call add_default("CRGRPTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRLIQTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Liquid tendency')
        if (carma_diags_file > 0) call add_default("CRLIQTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("CLDLIQTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Cloud liquid tendency')
        if (carma_diags_file > 0) call add_default("CLDLIQTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CLDICETC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Cloud ice tendency')
        if (carma_diags_file > 0) call add_default("CLDICETC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("RAINQMTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Rain tendency')
        if (carma_diags_file > 0) call add_default("RAINQMTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("SNOWQMTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Snow tendency')
        if (carma_diags_file > 0) call add_default("SNOWQMTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("CRDICENTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2/s', trim(carma_diags_packages(i))//', Detrained ice number tendency')
        if (carma_diags_file > 0) call add_default("CRDICENTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRSICENTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2/s', trim(carma_diags_packages(i))//',  In situ ice number tendency')
        if (carma_diags_file > 0) call add_default("CRSICENTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRGRPNTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2/s', trim(carma_diags_packages(i))//', Graupel number tendency')
        if (carma_diags_file > 0) call add_default("CRGRPNTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRLIQNTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2/s', trim(carma_diags_packages(i))//', Liquid number tendency')
        if (carma_diags_file > 0) call add_default("CRLIQNTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("NUMLIQTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2/s', trim(carma_diags_packages(i))//', Cloud liquid tendency')
        if (carma_diags_file > 0) call add_default("NUMLIQTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("NUMICETC_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2/s', trim(carma_diags_packages(i))//', Cloud ice tendency')
        if (carma_diags_file > 0) call add_default("NUMICETC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("NUMRAITC_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2/s', trim(carma_diags_packages(i))//', Rain tendency')
        if (carma_diags_file > 0) call add_default("NUMRAITC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("NUMSNOTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2/s', trim(carma_diags_packages(i))//', Snow tendency')
        if (carma_diags_file > 0) call add_default("NUMSNOTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("QBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', Q burden')
        if (carma_diags_file > 0) call add_default("QBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("CRDICEBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', Detrained ice burden')
        if (carma_diags_file > 0) call add_default("CRDICEBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRSICEBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//',  In situ ice burden')
        if (carma_diags_file > 0) call add_default("CRSICEBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRGRPBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', Graupel burden')
        if (carma_diags_file > 0) call add_default("CRGRPBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRLIQBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', Liquid burden')
        if (carma_diags_file > 0) call add_default("CRLIQBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("CLDLIQBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', Cloud liquid burden')
        if (carma_diags_file > 0) call add_default("CLDLIQBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CLDICEBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', Cloud ice burden')
        if (carma_diags_file > 0) call add_default("CLDICEBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("RAINQMBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', Rain burden')
        if (carma_diags_file > 0) call add_default("RAINQMBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("SNOWQMBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', Snow burden')
        if (carma_diags_file > 0) call add_default("SNOWQMBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("CRDICENBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2', trim(carma_diags_packages(i))//', Detrained ice number burden')
        if (carma_diags_file > 0) call add_default("CRDICENBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRSICENBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2', trim(carma_diags_packages(i))//',  In situ ice number burden')
        if (carma_diags_file > 0) call add_default("CRSICENBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRGRPNBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2', trim(carma_diags_packages(i))//', Graupel number burden')
        if (carma_diags_file > 0) call add_default("CRGRPNBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("CRLIQNBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2', trim(carma_diags_packages(i))//', Liquid number burden')
        if (carma_diags_file > 0) call add_default("CRLIQNBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("NUMLIQBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2', trim(carma_diags_packages(i))//', Cloud liquid number burden')
        if (carma_diags_file > 0) call add_default("NUMLIQBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("NUMICEBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2', trim(carma_diags_packages(i))//', Cloud ice number burden')
        if (carma_diags_file > 0) call add_default("NUMICEBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("NUMRAIBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2', trim(carma_diags_packages(i))//', Rain number burden')
        if (carma_diags_file > 0) call add_default("NUMRAIBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("NUMSNOBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', '#/m2', trim(carma_diags_packages(i))//', Snow number burden')
        if (carma_diags_file > 0) call add_default("NUMSNOBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("FLXVAPTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Boundary flux vapor')
        if (carma_diags_file > 0) call add_default("FLXVAPTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("FLXCNDTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Boundary flux condensate')
        if (carma_diags_file > 0) call add_default("FLXCNDTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("FLXICETC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Boundary flux ice')
        if (carma_diags_file > 0) call add_default("FLXICETC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("FLXSENTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'W/m2', trim(carma_diags_packages(i))//', Boundary flux sensible heat')
        if (carma_diags_file > 0) call add_default("FLXSENTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("QSF_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Surface Flux, H2O')
        if (carma_diags_file > 0) call add_default("QSF_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

      end do
    end if


    return
  end subroutine CARMA_InitializeModel


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009

  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday, &
                             is_perpetual
    use camsrfexch,       only: cam_in_t

    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    real(r8), intent(out)              :: surfaceFlux(pcols)    !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer      :: ncol                    ! number of columns in chunk

    ! Default return code.

    rc = rc_ok

    ncol = state%ncol    
    
    ! Add any surface flux here.

    surfaceFlux(:ncol) = 0.0_r8
    
    ! For emissions into the atmosphere, put the emission here.

    tendency(:ncol, :pver) = 0.0_r8
    
    return
  end subroutine CARMA_EmitParticle


  !! Sets the initial condition for CARMA aerosol particles. By default, there are no
  !! particles, but this routine can be overridden for models that wish to have an
  !! initial value.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009

  subroutine CARMA_InitializeParticle(carma, ielem, ibin, latvals, lonvals, mask, q, rc)

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plat, plev, plon

    implicit none
    
    type(carma_type), intent(in)  :: carma      !! the carma object
    integer,          intent(in)  :: ielem      !! element index
    integer,          intent(in)  :: ibin       !! bin index
    real(r8),         intent(in)  :: latvals(:) !! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) !! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    !! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     !! mass mixing ratio (gcol, lev)
    integer,          intent(out) :: rc         !! return code, negative indicates failure
        

    ! Default return code.
    rc = RC_OK

    ! Add initial condition here.
    
    return
  end subroutine CARMA_InitializeParticle


  !! Called at the end of the timestep after all the columns have been processed to
  !! to allow additional diagnostics that have been stored in pbuf to be output.
  !!
  !! NOTE: For the cloud model, we want to have budgets for water vapor, hydrometeors,
  !! and possibly energy. For the hydrometeors, we may want to look at both the CAM
  !! aummary fields (e.g. CLDLIQ, CLDICE) as weel as the CARMA fields. They should be
  !! kept in synch, but might diverge.
  !!
  !! NOTE: Output occurs a chunk at a time.
  !!
  !!  @version January-2023
  !!  @author  Chuck Bardeen
  subroutine CARMA_OutputBudgetDiagnostics(carma, icnst4elem, icnst4gas, state, ptend, dt, pname, rc, old_cflux, cflux, flx_vap, flx_cnd, flx_ice, flx_sen)
    use cam_history,  only: outfld
    use constituents, only: pcnst, cnst_get_ind
    
    type(carma_type), intent(in)         :: carma        !! the carma object
    integer, intent(in)                  :: icnst4elem(NELEM, NBIN) !! constituent index for a carma element
    integer, intent(in)                  :: icnst4gas(NGAS)         !! constituent index for a carma gas
    type(physics_state), intent(in)      :: state        !! Physics state variables - before pname
    type(physics_ptend), intent(in)      :: ptend        !! indivdual parameterization tendencies
    real(r8), intent(in)                 :: dt           !! timestep (s)
    character(*), intent(in)             :: pname        !! short name of the physics package
    integer, intent(out)                 :: rc           !! return code, negative indicates failure
    real(r8), optional, intent(in)       :: old_cflux(pcols,pcnst)  !! cam_in%clfux from before the timestep_tend
    real(r8), optional, intent(in)       :: cflux(pcols,pcnst)  !! cam_in%clfux from after the timestep_tend

    ! These match the inputs to check_energy_chng(), so they keep track of changes in
    ! water and energy that aren't part of the state vector.
    real(r8), optional, intent(in)       :: flx_vap(pcols)  !! boundary flux of vapor         (kg/m2/s)
    real(r8), optional, intent(in)       :: flx_cnd(pcols)  !! boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), optional, intent(in)       :: flx_ice(pcols)  !! boundary flux of ice           (m/s) (snow?)
    real(r8), optional, intent(in)       :: flx_sen(pcols)  !! boundary flux of sensible heat (w/m2)
 
    integer                              :: icol         !! column index
    integer                              :: ibin         !! bin index
    integer                              :: i
    integer                              :: ncols
    integer                              :: icnst        !! constituent index
    integer                              :: ienconc      !! concentration element index
    integer                              :: ncore        !! number of cores
    integer                              :: icorelem(NELEM) !! core element index
    real(r8)                             :: mair(pver)   !! Mass of air column (kg/m2)
    real(r8)                             :: rmass(NBIN)  !! Mass of bin (g)
    real(r8)                             :: crdicetend(pcols) !! Tendency detrained ice (kg/m2/s)
    real(r8)                             :: crsicetend(pcols) !! Tendency in situ ice (kg/m2/s)
    real(r8)                             :: crgrptend(pcols) !! Tendency graupel (kg/m2/s)
    real(r8)                             :: crliqtend(pcols) !! Tendency liquid (kg/m2/s)
    real(r8)                             :: h2otend(pcols)  !! Tendency water vapor (kg/m2/s)
    real(r8)                             :: cldliqtend(pcols) !! Tendency cloud liquid (kg/m2/s)
    real(r8)                             :: cldicetend(pcols) !! Tendency cloud ice (kg/m2/s)
    real(r8)                             :: rainqmtend(pcols) !! Tendency rain (kg/m2/s)
    real(r8)                             :: snowqmtend(pcols) !! Tendency snow (kg/m2/s)
    real(r8)                             :: crdicentend(pcols) !! Tendency detrained ice number (#/m2/s)
    real(r8)                             :: crsicentend(pcols) !! Tendency in situ ice number (#/m2/s)
    real(r8)                             :: crgrpntend(pcols) !! Tendency graupel number (#/m2/s)
    real(r8)                             :: crliqntend(pcols) !! Tendency liquid number (#/m2/s)
    real(r8)                             :: numliqtend(pcols) !! Tendency cloud liquid number (#/m2/s)
    real(r8)                             :: numicetend(pcols) !! Tendency cloud ice number (#/m2/s)
    real(r8)                             :: numraitend(pcols) !! Tendency rain number (#/m2/s)
    real(r8)                             :: numsnotend(pcols) !! Tendency snow number (#/m2/s)
    real(r8)                             :: tottend(pver)   !! Total Tendency in situ ice (kg/m2/s)
    real(r8)                             :: crdicebd(pcols) !! Burden detrained ice (kg/m2)
    real(r8)                             :: crsicebd(pcols) !! Burden in situ ice (kg/m2)
    real(r8)                             :: crgrpbd(pcols) !! Burden graupel (kg/m2)
    real(r8)                             :: crliqbd(pcols) !! Burden liquid (kg/m2)
    real(r8)                             :: h2obd(pcols)  !! Burden water vapor (kg/m2)
    real(r8)                             :: cldliqbd(pcols) !! Burden cloud liquid (kg/m2)
    real(r8)                             :: cldicebd(pcols) !! Burden cloud ice (kg/m2)
    real(r8)                             :: rainqmbd(pcols) !! Burden cloud liquid (kg/m2)
    real(r8)                             :: snowqmbd(pcols) !! Burden cloud ice (kg/m2)
    real(r8)                             :: crdicenbd(pcols) !! Burden detrained ice number (#/m2)
    real(r8)                             :: crsicenbd(pcols) !! Burden in situ ice number (#/m2)
    real(r8)                             :: crgrpnbd(pcols) !! Burden graupel number (#/m2)
    real(r8)                             :: crliqnbd(pcols) !! Burden liquid number (#/m2)
    real(r8)                             :: numliqbd(pcols) !! Burden cloud liquid number (#/m2)
    real(r8)                             :: numicebd(pcols) !! Burden cloud ice number (#/m2)
    real(r8)                             :: numraibd(pcols) !! Burden rain number (#/m2)
    real(r8)                             :: numsnobd(pcols) !! Burden snow number (#/m2)
    real(r8)                             :: ch2oflux(pcols)  !! Burden surface flux (#/m2)

    ! Default return code.
    rc = RC_OK
    
    h2otend(:)  = 0._r8
    crdicetend(:) = 0._r8
    crsicetend(:) = 0._r8
    crgrptend(:) = 0._r8
    crliqtend(:) = 0._r8
    cldliqtend(:) = 0._r8
    cldicetend(:) = 0._r8
    rainqmtend(:) = 0._r8
    snowqmtend(:) = 0._r8
    crdicentend(:) = 0._r8
    crsicentend(:) = 0._r8
    crgrpntend(:) = 0._r8
    crliqntend(:) = 0._r8
    numliqtend(:) = 0._r8
    numicetend(:) = 0._r8
    numraitend(:) = 0._r8
    numsnotend(:) = 0._r8
    
    h2obd(:)  = 0._r8
    crdicebd(:) = 0._r8
    crsicebd(:) = 0._r8
    crgrpbd(:) = 0._r8
    crliqbd(:) = 0._r8
    cldliqbd(:) = 0._r8
    cldicebd(:) = 0._r8
    rainqmbd(:) = 0._r8
    rainqmbd(:) = 0._r8
    crdicenbd(:) = 0._r8
    crsicenbd(:) = 0._r8
    crgrpnbd(:) = 0._r8
    crliqnbd(:) = 0._r8
    numliqbd(:) = 0._r8
    numicebd(:) = 0._r8
    numraibd(:) = 0._r8
    numsnobd(:) = 0._r8

    ! Add up the sulfate tendencies.
    ncols = state%ncol
    do icol = 1, ncols
      
      ! Get the air mass in the column
      !
      ! NOTE convert GRAV from cm/s2 to m/s2.
      mair(:) = state%pdel(icol,:) / (GRAV / 100._r8)
      
      do ibin = 1, nbin
      
        ! For CRDICE, CRGRP, and CRLIQ it is just the tendency for the concentration element.
        call CARMAGROUP_Get(carma, I_GRP_CRDICE, rc, ienconc=ienconc, rmass=rmass)
        icnst = icnst4elem(ienconc, ibin)

        if (ptend%lq(icnst)) then
          crdicetend(icol) = crdicetend(icol) + sum(ptend%q(icol,:,icnst) * mair(:))
          crdicentend(icol) = crdicentend(icol) + sum(ptend%q(icol,:,icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)
        end if

        crdicebd(icol) = crdicebd(icol) + sum(state%q(icol,:,icnst) * mair(:))
        crdicenbd(icol) = crdicenbd(icol) + sum(state%q(icol,:,icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)


        call CARMAGROUP_Get(carma, I_GRP_CRGRP, rc, ienconc=ienconc, rmass=rmass)
        icnst = icnst4elem(ienconc, ibin)

        if (ptend%lq(icnst)) then
          crgrptend(icol) = crgrptend(icol) + sum(ptend%q(icol,:,icnst) * mair(:))
          crgrpntend(icol) = crgrpntend(icol) + sum(ptend%q(icol,:,icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)
        end if

        crgrpbd(icol) = crgrpbd(icol) + sum(state%q(icol,:,icnst) * mair(:))
        crgrpnbd(icol) = crgrpnbd(icol) + sum(state%q(icol,:,icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)


        call CARMAGROUP_Get(carma, I_GRP_CRLIQ, rc, ienconc=ienconc, rmass=rmass)
        icnst = icnst4elem(ienconc, ibin)

        if (ptend%lq(icnst)) then
          crliqtend(icol) = crliqtend(icol) + sum(ptend%q(icol,:,icnst) * mair(:))
          crliqntend(icol) = crliqntend(icol) + sum(ptend%q(icol,:,icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)
        end if
        
        crliqbd(icol) = crliqbd(icol) + sum(state%q(icol,:,icnst) * mair(:))
        crliqnbd(icol) = crliqnbd(icol) + sum(state%q(icol,:,icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)
        
        
        ! For CRSICE, it is the difference in mass between the concentration element
        ! and the sum of the core masses.
        call CARMAGROUP_Get(carma, I_GRP_CRSICE, rc, ienconc=ienconc, rmass=rmass, ncore=ncore, icorelem=icorelem)
        icnst = icnst4elem(ienconc, ibin)

        tottend(:) = 0._r8
        if (ptend%lq(icnst)) then
          tottend(:) = ptend%q(icol, :, icnst) * mair(:)
          crsicentend(icol) = crsicentend(icol) + sum(tottend(:)) / (rmass(ibin) / 1.e3_f)
        end if
        crsicebd(icol) = crsicebd(icol) + sum(state%q(icol,:,icnst) * mair(:))
        crsicenbd(icol) = crsicenbd(icol) + sum(state%q(icol,:,icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)
          
        do i = 1, ncore
          icnst = icnst4elem(icorelem(i), ibin)
          if (ptend%lq(icnst)) then
            tottend(:) = tottend(:) - ptend%q(icol,:,icnst) * mair(:)
          end if
          crsicebd(icol) = crsicebd(icol) - sum(state%q(icol,:,icnst) * mair(:))
        end do
        
        crsicetend(icol) = crsicetend(icol) + sum(tottend(:))
      end do
      
      ! Calculate the water vapor change.
      icnst = icnst4gas(I_GAS_H2O)
      if (ptend%lq(icnst)) then
        h2otend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if
      h2obd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

      call cnst_get_ind("CLDLIQ", icnst)
      if (ptend%lq(icnst)) then
        cldliqtend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if
      cldliqbd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

      call cnst_get_ind("CLDICE", icnst)
      if (ptend%lq(icnst)) then
        cldicetend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if
      cldicebd(icol) = sum(state%q(icol,:,icnst) * mair(:))

      call cnst_get_ind("NUMLIQ", icnst)
      if (ptend%lq(icnst)) then
        numliqtend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if
      numliqbd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

      call cnst_get_ind("NUMICE", icnst)
      if (ptend%lq(icnst)) then
        numicetend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if
      numicebd(icol) = sum(state%q(icol,:,icnst) * mair(:))

      call cnst_get_ind("RAINQM", icnst)
      if (ptend%lq(icnst)) then
        rainqmtend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if
      rainqmbd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

      call cnst_get_ind("SNOWQM", icnst)
      if (ptend%lq(icnst)) then
        snowqmtend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if
      snowqmbd(icol) = sum(state%q(icol,:,icnst) * mair(:))

      call cnst_get_ind("NUMRAI", icnst)
      if (ptend%lq(icnst)) then
        numraitend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if
      numraibd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

      call cnst_get_ind("NUMSNO", icnst)
      if (ptend%lq(icnst)) then
        numsnotend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if
      numsnobd(icol) = sum(state%q(icol,:,icnst) * mair(:))
      
      ! Calculate the surface flux of water
      if (present(cflux)) then 
        call cnst_get_ind("Q", icnst)
        ch2oflux(icol) = cflux(icol,icnst)
      end if
    end do

    
    ! Output the cloud and water tendencies for this physics package.
    call outfld("QTC_"//trim(pname), h2otend(:), pcols, state%lchnk)

    call outfld("CRDICETC_"//trim(pname), crdicetend(:), pcols, state%lchnk)
    call outfld("CRSICETC_"//trim(pname), crsicetend(:), pcols, state%lchnk)
    call outfld("CRGRPTC_"//trim(pname), crgrptend(:), pcols, state%lchnk)
    call outfld("CRLIQTC_"//trim(pname), crliqtend(:), pcols, state%lchnk)

    call outfld("CLDLIQTC_"//trim(pname), cldliqtend(:), pcols, state%lchnk)
    call outfld("CLDICETC_"//trim(pname), cldicetend(:), pcols, state%lchnk)

    call outfld("RAINQMTC_"//trim(pname), rainqmtend(:), pcols, state%lchnk)
    call outfld("SNOWQMTC_"//trim(pname), snowqmtend(:), pcols, state%lchnk)

    call outfld("CRDICENTC_"//trim(pname), crdicentend(:), pcols, state%lchnk)
    call outfld("CRSICENTC_"//trim(pname), crsicentend(:), pcols, state%lchnk)
    call outfld("CRGRPNTC_"//trim(pname), crgrpntend(:), pcols, state%lchnk)
    call outfld("CRLIQNTC_"//trim(pname), crliqntend(:), pcols, state%lchnk)

    call outfld("NUMLIQTC_"//trim(pname), numliqtend(:), pcols, state%lchnk)
    call outfld("NUMICETC_"//trim(pname), numicetend(:), pcols, state%lchnk)

    call outfld("NUMRAITC_"//trim(pname), numraitend(:), pcols, state%lchnk)
    call outfld("NUMSNOTC_"//trim(pname), numsnotend(:), pcols, state%lchnk)


    call outfld("QBD_"//trim(pname), h2obd(:), pcols, state%lchnk)

    call outfld("CRDICEBD_"//trim(pname), crdicebd(:), pcols, state%lchnk)
    call outfld("CRSICEBD_"//trim(pname), crsicebd(:), pcols, state%lchnk)
    call outfld("CRGRPBD_"//trim(pname), crgrpbd(:), pcols, state%lchnk)
    call outfld("CRLIQBD_"//trim(pname), crliqbd(:), pcols, state%lchnk)

    call outfld("CLDLIQBD_"//trim(pname), cldliqbd(:), pcols, state%lchnk)
    call outfld("CLDICEBD_"//trim(pname), cldicebd(:), pcols, state%lchnk)

    call outfld("RAINQMBD_"//trim(pname), rainqmbd(:), pcols, state%lchnk)
    call outfld("SNOWQMBD_"//trim(pname), snowqmbd(:), pcols, state%lchnk)

    call outfld("CRDICENBD_"//trim(pname), crdicenbd(:), pcols, state%lchnk)
    call outfld("CRSICENBD_"//trim(pname), crsicenbd(:), pcols, state%lchnk)
    call outfld("CRGRPNBD_"//trim(pname), crgrpnbd(:), pcols, state%lchnk)
    call outfld("CRLIQNBD_"//trim(pname), crliqnbd(:), pcols, state%lchnk)

    call outfld("NUMLIQBD_"//trim(pname), numliqbd(:), pcols, state%lchnk)
    call outfld("NUMICEBD_"//trim(pname), numicebd(:), pcols, state%lchnk)

    call outfld("NUMRAIBD_"//trim(pname), numraibd(:), pcols, state%lchnk)
    call outfld("NUMSNOBD_"//trim(pname), numsnobd(:), pcols, state%lchnk)

    ! Convert from m/s to kg/m2/s?
    if (present(flx_vap)) call outfld("FLXVAPTC_"//trim(pname), flx_vap(:ncols), pcols, state%lchnk)
    if (present(flx_cnd)) call outfld("FLXCNDTC_"//trim(pname), flx_cnd(:ncols) * 1000._f, pcols, state%lchnk)
    if (present(flx_ice)) call outfld("FLXICETC_"//trim(pname), flx_ice(:ncols) * 1000._f, pcols, state%lchnk)
    if (present(flx_sen)) call outfld("FLXSENTC_"//trim(pname), flx_sen(:ncols), pcols, state%lchnk)    
    if (present(cflux)) call outfld("QSF_"//trim(pname), ch2oflux(:ncols) * 1000._f, pcols, state%lchnk) 
      
    return
  end subroutine CARMA_OutputBudgetDiagnostics
  
    !! Called at the end of the timestep after all the columns have been processed to
  !! to allow additional diagnostics that have been stored in pbuf to be output.
  !!
  !! NOTE: Output occurs a chunk at a time.
  !!
  !!  @version January-2023
  !!  @author  Chuck Bardeen
  subroutine CARMA_OutputDiagnostics(carma, icnst4elem, state, ptend, pbuf, cam_in, rc, gpdiags)
    use cam_history,   only: outfld
    use constituents,  only: pcnst, cnst_get_ind
    use camsrfexch,    only: cam_in_t

    ! This is a hack as these are currently defined in carma_intr.F90, but I don't want to
    ! include that file. They may need to be moved into another file, but for now they
    ! are duplicated here.
    integer, parameter :: NGPDIAGS    = 27  ! Number of particle diagnostics ...
    integer, parameter :: GPDIAGS_ND  =  1  ! Number density
    integer, parameter :: GPDIAGS_AD  =  2  ! Surface area density
    integer, parameter :: GPDIAGS_MD  =  3  ! Mass density
    integer, parameter :: GPDIAGS_RE  =  4  ! Effective Radius
    integer, parameter :: GPDIAGS_RM  =  5  ! Mitchell [2002] Effective Radius
    integer, parameter :: GPDIAGS_JN  =  6  ! Nucleation Rate
    integer, parameter :: GPDIAGS_MR  =  7  ! Mass Mixing Ratio
    integer, parameter :: GPDIAGS_EX  =  8  ! Extinction
    integer, parameter :: GPDIAGS_OD  =  9  ! Optical Depth
    integer, parameter :: GPDIAGS_VM  = 10  ! Mass Weighted Fall Velocity
    integer, parameter :: GPDIAGS_PA  = 11  ! Projected Area
    integer, parameter :: GPDIAGS_AR  = 12  ! Area Ratio
    integer, parameter :: GPDIAGS_ADS = 13  ! Surface area density, small
    integer, parameter :: GPDIAGS_VDS = 14  ! Volume density, small 
    integer, parameter :: GPDIAGS_RES = 15  ! Effective radius, small
    integer, parameter :: GPDIAGS_MRS = 16  ! Mixing Ratio, small
    integer, parameter :: GPDIAGS_NDS = 17  ! Number Density, small
    integer, parameter :: GPDIAGS_ADM = 18  ! Surface area density, medium
    integer, parameter :: GPDIAGS_VDM = 19  ! Volume density, medium 
    integer, parameter :: GPDIAGS_REM = 20  ! Effective radius, medium
    integer, parameter :: GPDIAGS_MRM = 21  ! Mixing Ratio, medium
    integer, parameter :: GPDIAGS_NDM = 22  ! Number Density, medium
    integer, parameter :: GPDIAGS_ADL = 23  ! Surface area density, large
    integer, parameter :: GPDIAGS_VDL = 24  ! Volume density, large
    integer, parameter :: GPDIAGS_REL = 25  ! Effective radius, large
    integer, parameter :: GPDIAGS_MRL = 26  ! Mixing Ratio, medium
    integer, parameter :: GPDIAGS_NDL = 27  ! Number Density, medium 

    
    type(carma_type), intent(in)         :: carma        !! the carma object
    integer, intent(in)                  :: icnst4elem(NELEM, NBIN) !! constituent index for a carma element
    type(physics_state), intent(in)      :: state        !! Physics state variables - before CARMA
    type(physics_ptend), intent(in)      :: ptend        !! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer, intent(in)   :: pbuf(:)  !! physics buffer
    type(cam_in_t), intent(in)           :: cam_in       !! surface inputs
    integer, intent(out)                 :: rc           !! return code, negative indicates failure
    real(r8), intent(in)                 :: gpdiags(pcols, pver, NGROUP, NGPDIAGS) ! Carma group diagnostic output

 
    integer                              :: icol         !! column index
    integer                              :: ibin         !! bin index
    real(r8), pointer, dimension(:,:)    :: soacm        !! aerosol tendency due to gas-aerosol exchange  kg/kg/s
    real(r8), pointer, dimension(:,:)    :: soapt        !! aerosol tendency due to no2 photolysis  kg/kg/s
    character(len=16)                    :: binname      !! names bins
    integer                              :: i
    integer                              :: icnst        !! constituent index
    integer                              :: ienconc      !! concentration element index
    integer                              :: ncore        !! number of cores
    integer                              :: icorelem(NELEM) !! core element index
    real(r8)                             :: mair(pver)   !! Mass of air column (kg/m2)
    real(r8)                             :: rmass(NBIN)  !! Mass of bin (g)
    real(r8)                             :: h2obd(pcols)    !! Burden H2O gas (kg/m2)
    real(r8)                             :: h2osf(pcols)    !! H2O surface flux (kg/m2/s)
    real(r8)                             :: crdicebd(pcols) !! Burden Detrained ice (kg/m2)
    real(r8)                             :: crsicebd(pcols) !! Burden In Situ ice (kg/m2)
    real(r8)                             :: crgrpbd(pcols)  !! Burden Graupel (kg/m2)
    real(r8)                             :: crliqbd(pcols)  !! Burden Liquid (kg/m2)
    real(r8)                             :: crscorebd(pcols)!! Burden Sulfate Core (kg/m2)
    real(r8)                             :: crdcorebd(pcols)!! Burden Dust Core (kg/m2)
    real(r8)                             :: cldliqbd(pcols) !! Burden Cloud Liquid (kg/m2)
    real(r8)                             :: cldicebd(pcols) !! Burden Cloud Ice (kg/m2)
    real(r8)                             :: rainqmbd(pcols) !! Burden Rain (kg/m2)
    real(r8)                             :: snowqmbd(pcols) !! Burden Snow (kg/m2)
    real(r8)                             :: crdicenbd(pcols)!! Number Burden Detrained ice (kg/m2)
    real(r8)                             :: crsicenbd(pcols)!! Number Burden In Situ ice (kg/m2)
    real(r8)                             :: crgrpnbd(pcols) !! Number Burden Graupel (kg/m2)
    real(r8)                             :: crliqnbd(pcols) !! Number Burden Liquid (kg/m2)
    real(r8)                             :: numliqbd(pcols) !! Number Burden Cloud Liquid (#/m2)
    real(r8)                             :: numicebd(pcols) !! Number Burden Cloud Ice (#/m2)
    real(r8)                             :: numraibd(pcols) !! Number Burden Rain (#/m2)
    real(r8)                             :: numsnobd(pcols) !! Number Burden Snow (#/m2)
    real(r8)                             :: lresbd(pcols)   !! Burden liquid, small effective radius (um)
    real(r8)                             :: lrembd(pcols)   !! Burden liquid, medium effective radius (um)
    real(r8)                             :: lrelbd(pcols)   !! Burden liquid, large effective radius (um)
    real(r8)                             :: iresbd(pcols)   !! Burden ice, small effective radius (um)
    real(r8)                             :: irembd(pcols)   !! Burden ice, medium effective radius (um)
    real(r8)                             :: irelbd(pcols)   !! Burden ice, large effective radius (um)
    real(r8)                             :: dz(pver)       !! gridbox height (cm)
    real(r8)                             :: adtot
    character(len=16)                    :: shortname

    ! Default return code.
    rc = RC_OK
    
    ! Only do these if budget diagnostics have been turned on.
    if (carma_do_budget_diags) then

      ! Output the cloud burdens.
      crdicebd(:)  = 0._r8
      crsicebd(:)  = 0._r8
      crgrpbd(:)   = 0._r8
      crliqbd(:)   = 0._r8
      crscorebd(:) = 0._r8  
      crdcorebd(:) = 0._r8
      cldliqbd(:)  = 0._r8  
      cldicebd(:)  = 0._r8  
      rainqmbd(:)  = 0._r8  
      snowqmbd(:)  = 0._r8  
      crdicenbd(:) = 0._r8
      crsicenbd(:) = 0._r8
      crgrpnbd(:)  = 0._r8
      crliqnbd(:)  = 0._r8
      numliqbd(:)  = 0._r8  
      numicebd(:)  = 0._r8  
      numraibd(:)  = 0._r8  
      numsnobd(:)  = 0._r8  
      h2obd(:)     = 0._r8
      lresbd(:)    = CAM_FILL
      lrembd(:)    = CAM_FILL
      lrelbd(:)    = CAM_FILL
      iresbd(:)    = CAM_FILL
      irembd(:)    = CAM_FILL
      irelbd(:)    = CAM_FILL
    
      ! Add up the sulfate tendencies.
      do icol = 1, state%ncol
      
        ! Get the air mass in the column
        !
        ! NOTE convert GRAV from cm/s2 to m/s2.
        mair(:) = state%pdel(icol,:) / (GRAV / 100._r8)
      
        do ibin = 1, nbin
      
          ! For CRDICE, CRGRP, and CRLIQ it is just the tendency for the concentration element.
          call CARMAGROUP_Get(carma, I_GRP_CRDICE, rc, ienconc=ienconc, rmass=rmass)
          icnst = icnst4elem(ienconc, ibin)

          crdicebd(icol) = crdicebd(icol) + sum(state%q(icol,:,icnst) * mair(:))
          crdicenbd(icol) = crdicenbd(icol) + sum(state%q(icol,:,icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)


          call CARMAGROUP_Get(carma, I_GRP_CRGRP, rc, ienconc=ienconc, rmass=rmass)
          icnst = icnst4elem(ienconc, ibin)

          crgrpbd(icol) = crgrpbd(icol) + sum(state%q(icol,:,icnst) * mair(:))
          crgrpnbd(icol) = crgrpnbd(icol) + sum(state%q(icol,:,icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)


          call CARMAGROUP_Get(carma, I_GRP_CRLIQ, rc, ienconc=ienconc, rmass=rmass)
          icnst = icnst4elem(ienconc, ibin)

          crliqbd(icol) = crliqbd(icol) + sum(state%q(icol,:,icnst) * mair(:))
          crliqnbd(icol) = crliqnbd(icol) + sum(state%q(icol,:,icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)

         
          ! For CRSICE, it is the difference in mass between the concentration element
          ! and the sum of the core masses.
          call CARMAGROUP_Get(carma, I_GRP_CRSICE, rc, ienconc=ienconc, rmass=rmass, ncore=ncore, icorelem=icorelem)
          icnst = icnst4elem(ienconc, ibin)

          crsicebd(icol) = crsicebd(icol) + sum(state%q(icol, :, icnst) * mair(:))
          crsicenbd(icol) = crsicenbd(icol) + sum(state%q(icol, :, icnst) * mair(:)) / (rmass(ibin) / 1.e3_f)
          
          do i = 1, ncore
            icnst = icnst4elem(icorelem(i), ibin)
            crsicebd(icol) = crsicebd(icol) - sum(state%q(icol,:,icnst) * mair(:))

            call CARMAELEMENT_Get(carma, icorelem(i), rc, shortname=shortname)
            if (shortname .eq. "CRCORE") then
              crscorebd(icol) = crscorebd(icol) + sum(state%q(icol,:,icnst) * mair(:))
            else if (shortname .eq. "CRDCOR") then
              crdcorebd(icol) = crdcorebd(icol) + sum(state%q(icol,:,icnst) * mair(:))
            end if 

          end do
        end do
        
        ! Total the effective radii in the different size ranges in the vertical
        ! for ice and liquid.
        dz(:) = (state%zi(icol,1:pver) - state%zi(icol,2:pver+1)) * 100._r8

        adtot = sum(gpdiags(icol,:, I_GRP_CRLIQ, GPDIAGS_ADS) * dz(:))
        if (adtot .gt. 0._r8) then
          lresbd(icol) = sum(gpdiags(icol,:, I_GRP_CRLIQ, GPDIAGS_VDS) * dz(:)) / adtot * 3._r8
        end if

        adtot = sum(gpdiags(icol,:, I_GRP_CRLIQ, GPDIAGS_ADM) * dz(:))
        if (adtot .gt. 0._r8) then
          lrembd(icol) = sum(gpdiags(icol,:, I_GRP_CRLIQ, GPDIAGS_VDM) * dz(:)) / adtot * 3._r8
        end if

        adtot = sum(gpdiags(icol,:, I_GRP_CRLIQ, GPDIAGS_ADL) * dz(:))
        if (adtot .gt. 0._r8) then
          lrelbd(icol) = sum(gpdiags(icol,:, I_GRP_CRLIQ, GPDIAGS_VDL) * dz(:)) / adtot * 3._r8
        end if

        adtot = sum((gpdiags(icol,:, I_GRP_CRSICE, GPDIAGS_ADS) + gpdiags(icol,:, I_GRP_CRDICE, GPDIAGS_ADS)) * dz(:))
        if (adtot .gt. 0._r8) then
          iresbd(icol) = sum((gpdiags(icol,:, I_GRP_CRSICE, GPDIAGS_VDS) + gpdiags(icol,:, I_GRP_CRDICE, GPDIAGS_VDS)) * dz(:)) / adtot * 3._r8
        end if

        adtot = sum((gpdiags(icol,:, I_GRP_CRSICE, GPDIAGS_ADM) + gpdiags(icol,:, I_GRP_CRDICE, GPDIAGS_ADM)) * dz(:))
        if (adtot .gt. 0._r8) then
          irembd(icol) = sum((gpdiags(icol,:, I_GRP_CRSICE, GPDIAGS_VDM) + gpdiags(icol,:, I_GRP_CRDICE, GPDIAGS_VDM)) * dz(:)) / adtot * 3._r8
        end if

        adtot = sum((gpdiags(icol,:, I_GRP_CRSICE, GPDIAGS_ADL) + gpdiags(icol,:, I_GRP_CRDICE, GPDIAGS_ADL)) * dz(:))
        if (adtot .gt. 0._r8) then
          irelbd(icol) = sum((gpdiags(icol,:, I_GRP_CRSICE, GPDIAGS_VDL) + gpdiags(icol,:, I_GRP_CRDICE, GPDIAGS_VDL)) * dz(:)) / adtot * 3._r8
        end if

        
        ! Calculate the Water burden.
        call cnst_get_ind("Q", icnst)
        h2obd(icol) = sum(state%q(icol,:,icnst) * mair(:))      
        h2osf(icol) = cam_in%cflx(icol,icnst)      

        call cnst_get_ind("CLDLIQ", icnst)
        cldliqbd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

        call cnst_get_ind("CLDICE", icnst)
        cldicebd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

        call cnst_get_ind("NUMLIQ", icnst)
        numliqbd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

        call cnst_get_ind("NUMICE", icnst)
        numicebd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

        call cnst_get_ind("RAINQM", icnst)
        rainqmbd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

        call cnst_get_ind("SNOWQM", icnst)
        snowqmbd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

        call cnst_get_ind("NUMRAI", icnst)
        numraibd(icol) = sum(state%q(icol,:,icnst) * mair(:))      

        call cnst_get_ind("NUMSNO", icnst)
        numsnobd(icol) = sum(state%q(icol,:,icnst) * mair(:))      
      end do
    
      ! Output the total aerosol and gas burdens and the aerosol fluxes.
      call outfld("QBD", h2obd(:), pcols, state%lchnk)

      call outfld("CRDICEBD", crdicebd(:), pcols, state%lchnk)
      call outfld("CRSICEBD", crsicebd(:), pcols, state%lchnk)
      call outfld("CRGRPBD", crgrpbd(:), pcols, state%lchnk)
      call outfld("CRLIQBD", crliqbd(:), pcols, state%lchnk)
      
      call outfld("CRSCOREBD", crscorebd(:), pcols, state%lchnk)
      call outfld("CRDCOREBD", crdcorebd(:), pcols, state%lchnk)

      call outfld("RAINQMBD", rainqmbd(:), pcols, state%lchnk)
      call outfld("SNOWQMBD", snowqmbd(:), pcols, state%lchnk)

      call outfld("CLDLIQBD", cldliqbd(:), pcols, state%lchnk)
      call outfld("CLDICEBD", cldicebd(:), pcols, state%lchnk)

      call outfld("CRDICENBD", crdicenbd(:), pcols, state%lchnk)
      call outfld("CRSICENBD", crsicenbd(:), pcols, state%lchnk)
      call outfld("CRGRPNBD", crgrpnbd(:), pcols, state%lchnk)
      call outfld("CRLIQNBD", crliqnbd(:), pcols, state%lchnk)

      call outfld("NUMLIQBD", numliqbd(:), pcols, state%lchnk)
      call outfld("NUMICEBD", numicebd(:), pcols, state%lchnk)

      call outfld("NUMRAIBD", numraibd(:), pcols, state%lchnk)
      call outfld("NUMSNOBD", numsnobd(:), pcols, state%lchnk)

      call outfld("CRLIQRESBD", lresbd(:), pcols, state%lchnk)
      call outfld("CRLIQREMBD", lrembd(:), pcols, state%lchnk)
      call outfld("CRLIQRELBD", lrelbd(:), pcols, state%lchnk)

      call outfld("CRICERESBD", iresbd(:), pcols, state%lchnk)
      call outfld("CRICEREMBD", irembd(:), pcols, state%lchnk)
      call outfld("CRICERELBD", irelbd(:), pcols, state%lchnk)

      ! Output the total aerosol and gas burdens and the aerosol fluxes.
      call outfld("QSF", h2osf(:), pcols, state%lchnk)

    end if
    
    return
  end subroutine CARMA_OutputDiagnostics  
  
    
  !!  Called after wet deposition has been performed. Allows the specific model to add
  !!  wet deposition of CARMA aerosols to the aerosols being communicated to the surface.
  !!
  !!  @version July-2011 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_WetDeposition(carma, ielem, ibin, sflx, cam_out, state, rc)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: ielem       !! element index
    integer, intent(in)                  :: ibin        !! bin index
    real(r8), intent(in)                 :: sflx(pcols) !! surface flux (kg/m2/s)
    type(cam_out_t), intent(inout)       :: cam_out     !! cam output to surface models
    type(physics_state), intent(in)      :: state       !! physics state variables
    integer, intent(out)                 :: rc          !! return code, negative indicates failure
    
    integer    :: icol
 
    ! Default return code.
    rc = RC_OK
    
    return
  end subroutine CARMA_WetDeposition
  
  
  ! Using the specified parameters for the gamma distribution, determine the mass mixing ratio of particles

  subroutine CARMA_GetMmrFromGamma(carma, r, dr, rmass, qic, nic, mu, lambda, mmr, rc)
    use shr_spfn_mod, only           : gamma => shr_spfn_gamma

    implicit none
    
    type(carma_type), intent(in)       :: carma           !! the carma object
    real(kind=f), intent(in)           :: r(NBIN)         !! bin mean radius
    real(kind=f), intent(in)           :: dr(NBIN)        !! bin radius width
    real(kind=f), intent(in)           :: rmass(NBIN)     !! bin mass
    real(r8), intent(in)               :: qic(pver)       !! in-cloud cloud liquid mixing ratio
    real(r8), intent(in)               :: nic(pver)       !! in-cloud droplet number conc
    real(r8), intent(in)               :: mu(pver)        !! spectral width parameter of droplet size distr
    real(r8), intent(in)               :: lambda(pver)    !! slope of cloud liquid size distr
    real(r8), intent(out)              :: mmr(NBIN,pver)  !! elements mass mixing ratio
    integer, intent(out)               :: rc              !! return code, negative indicates failure
    
    integer                            :: k               ! z index
    integer                            :: ibin            ! bin index
    real(kind=f)                       :: totalMass       ! mmr of all particles (kg/kg)
    real(kind=f)                       :: n               ! number of particles (#/kg)
    real(kind=f)                       :: n0              ! number parameter for gamma distribution
    real(kind=f)                       :: d(NBIN)         ! bin diameter (m)
    real(kind=f)                       :: dd(NBIN)        ! diameter width of bin (m)


  
    ! Default return code.
    rc = RC_OK
    
    ! Their equations are in terms of diameter (in m)
    d(:)  = 2._r8 * r(:) * 1e-2_r8
    dd(:) = 2._r8 * dr(:) * 1e-2_r8

    do k = 1, pver
    
      ! From Morisson & Gettelman [2008] and cldwat2m
      !
      ! If there is a small mass, then there are no particles.

      if (qic(k) <= qsmall) then
        mmr(:, k) = 0._r8
      else

        !print*, 'JAS', nic(k), lambda(k), mu(k)

        n0 = (nic(k) * (lambda(k) ** (mu(k) + 1._r8)) / (gamma(mu(k) + 1._r8)))    
  
      
        ! Iterate over the bins.
        !
        ! NOTE: Just the functional fit can go negative for some bins with larger diameter, but this is not physical.
        do ibin = 1, NBIN
          n = n0 * (d(ibin)**mu(k)) * exp(-lambda(k) * d(ibin)) * dd(ibin)
          mmr(ibin, k) = n * rmass(ibin) * 1e-3_r8
        end do
      
        ! Adjust the number density so that we don't create mass. This will adjust for
        ! problems fitting the size distribution and for differences in the assumptions
        ! of the bulk density of the particles.
        totalMass = sum(mmr(:, k))
        if (totalMass /= 0._r8) then
          mmr(:, k) = mmr(:, k) * (qic(k) / totalMass)
        else
          mmr(:, k) = 0._r8
        end if
      end if
    end do

    return
  end subroutine CARMA_GetMmrFromGamma


  !! Determine the total cloud ice concentration and number stored in the bins that represent
  !! water within the Carma model.
  !!
  !!  @version Nov-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_GetTotalIceAndSnow(carma, cstate, iceMass, iceNumber, snowSurface, rc, iceRe)

    implicit none
    
    type(carma_type), intent(in)         :: carma     !! the carma object
    type(carmastate_type), intent(inout) :: cstate    !! the carma state object
    real(kind=f), intent(out)            :: iceMass(pver)      !! ice mass mixing ratio (kg/kg)
    real(kind=f), intent(out)            :: iceNumber(pver)    !! ice number mixing ratio (#/kg)
    real(kind=f), intent(out)            :: snowSurface        !! snow on surface (kg/m2)
    integer, intent(out)                 :: rc        !! return code, negative indicates failure

    Real (kind=f), Intent (out), Optional :: iceRe (pver)  !! ice effective radius (m)
 
    integer                              :: LUNOPRT              ! logical unit number for output
    logical                              :: do_print             ! do print output?

    integer                              :: igroup    ! group index
    integer                              :: ielem     ! element index
    integer                              :: ibin      ! bin index
    integer                              :: iz        ! vertical index
    integer                              :: icore     ! core index
    integer                              :: icorelem(NELEM) ! core indexes for group
    integer                              :: ncore     ! number of core elements
    
    real(kind=f)                         :: coreMass(pver)     ! core mass mixing ratio (kg/kg)
    real(kind=f)                         :: coreSurface        ! core on surface (kg/kg)

    real(kind=f)                         :: mmr(pver)          ! mass mixing ratio (#/kg)
    real(kind=f)                         :: mmrcore(pver)      ! core mass mixing ratio (#/kg)
    real(kind=f)                         :: nmr(pver)          ! number mixing ratio (#/kg)
    real(kind=f)                         :: r(NBIN)            ! radius (cm)
    real(kind=f)                         :: rmass(NBIN)        ! mass (g)
    real(kind=f)                         :: rrat(NBIN)         ! particle maximum radius ratio ()
    real(kind=f)                         :: arat(NBIN)         ! particle area ratio ()

    real(kind=f)                         :: sfc  ! surface mass (kg/m2)

    real(kind=f)                         :: sfccore            ! core surface mass (kg/m2)
    real(kind=f)                         :: nd(pver)           ! number density (#/cm3)
    real(kind=f)                         :: pa(pver)           ! projected area (cm2)
    real(kind=f)                         :: md(pver)           ! mass density (g/cm3)

    !+++
    real(kind=f) :: t(pver)             ! temperature (K)
    integer      :: k             ! vertical index


    ! Default return code.

    rc = rc_ok
    
    call CARMA_Get(carma, rc, do_print=do_print, lunoprt=LUNOPRT)
    if (rc < rc_ok) call endrun ( 'CARMA_GetTotalIceAndSnow :: Carma_Get failed') 
    
    iceMass(:)     = 0._f
    iceNumber(:)   = 0._f
    snowSurface    = 0._f
    pa(:)          = 0._f
    md(:)          = 0._f

    if (present(iceRe)) iceRe(:) = 0._f
      
      
    ! Detrained Ice, Aged
    igroup = I_GRP_CRDICE
    ielem  = I_ELEM_CRDICE
    
    call CARMAGROUP_Get(carma, igroup, rc, r=r, rmass=rmass, arat=arat, rrat=rrat)
    if (rc < rc_ok) call endrun ( 'GetTotalIceAndSnow :: CarmaGroup_Get failed')
    
    do ibin = 1, NBIN
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, nmr=nmr, surface=sfc, NumberDensity=nd)
      if (rc < rc_ok) call endrun ( 'GetTotalIceAndSnow :: CarmaState_GetBin failed')

      iceMass(:)   = iceMass(:)   + mmr(:)
      iceNumber(:) = iceNumber(:) + nmr(:)
      where (nd(:) > SMALL_PC)
      
        ! NOTE: This is following the definition of Dave Mitchell for
        ! effective diameter, Mitchell [2002], which indicates it needs to be
        ! scaled based on the effective ice density
        pa(:) = pa(:) + nd(:) * PI * ((r(ibin) * rrat(ibin)) ** 2) * arat(ibin)
        md(:) = md(:) + nd(:) * rmass(ibin)
      end where
    
      ! The particles that sedimented out of the bottom layer need to be included
      ! in the mass of snow
      snowSurface = snowSurface + sfc
    end do

    ! Detrained Ice, Fresh
    !
    ! NOTE: Fresh means being detrained this timestep (pcd) that hasn't
    ! already been put in the ice mass (pc). 
    do ibin = 1, NBIN
      call CARMASTATE_GetDetrain(cstate, ielem, ibin, mmr, rc, nmr=nmr, numberDensity=nd)
      if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMASTATE_GetBin failed.')

      ! Only calculate snow if CARMA is responsible for the cloud ice.
      iceMass(:)   = iceMass(:)   + mmr(:)
      iceNumber(:) = iceNumber(:) + nmr(:)
      where (nd(:) > SMALL_PC)
        pa(:) = pa(:) + nd(:) * PI * ((r(ibin) * rrat(ibin))**2) * arat(ibin)
        md(:) = md(:) + nd(:) * rmass(ibin)
      end where
    end do
    
    ! In-situ ice
    igroup = I_GRP_CRSICE
    ielem  = I_ELEM_CRSICE
    
    call CARMAGROUP_Get(carma, igroup, rc, r=r, ncore=ncore, icorelem=icorelem, arat=arat, rrat=rrat)
    if (rc < rc_ok) call endrun ('GetTotalIceAndSnow :: CarmaGroup_Get failed.')
    
    do ibin = 1, NBIN
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, nmr=nmr, surface=sfc, numberDensity=nd)
      if (rc < rc_ok) call endrun ('Carma_GetTotalIceAndSnow :: CarmaState_GetBin failed.')
      
      ! Determine how much of the mmr is related to core mass. This needs to
      ! be subtracted to get the amount of water in the ice.
      coreMass(:) = 0.0_f
      coreSurface = 0.0_f
      
      do icore = 1, ncore
        call CARMASTATE_GetBin(cstate, icorelem(icore), ibin, mmrcore, rc, surface=sfccore)
        if (rc < rc_ok) call endrun ( 'GetTotalIceAndSnow :: CarmaState_GetBin failed.')
      
        coreMass(:) = coreMass(:) + mmrcore(:)
        coreSurface = coreSurface + sfccore
      end do
      
      ! The core mass can't be more than the particle mass. If so, this indicates
      ! that are problem happened, perhaps during advection and the particle masses
      ! should be ignored. This should never happen from CARMA itself.
      if (carma_do_mass_fix) then
        do iz = 1, pver
        
          if (coreMass(iz) > mmr(iz)) then
            if (carma_do_mass_fix) then
            
              if (carma_do_print_fix .and. do_print) write(LUNOPRT,*) &
                 "  GetTotalIceAndSnow::WARNING - Adjusting particle for core mass error", &
                 iz, ielem, ibin, mmr(iz), coreMass(iz)

              ! It is hard to know what the right fix should be. You could reset
              ! the particle mass to the coremass, but this will create lots of
              ! small particles. It may be safer just to zero out both the particle
              ! count and all of the core masses, assuming that this is a particle
              ! that was created by diffusion in the transport and shouldn't really exist.
              mmr(iz) = coreMass(iz)
            end if
          end if
        end do
            
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr, rc)
        if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMASTATE_SetBin failed.')
      end if
            
      iceMass(:)   = iceMass(:)   + mmr(:) - coreMass(:)
      iceNumber(:) = iceNumber(:) + nmr(:)

      where (nd(:) > SMALL_PC)
        pa(:) = pa(:) + nd(:) * PI * ((r(ibin) * rrat(ibin)) ** 2 ) * arat(ibin)
        md(:)  = md(:) + nd(:) * rmass(ibin)
      end where
    
    
      ! The particles that sedimented out of the bottom layer need to be included
      ! in the mass of snow.
      snowSurface = snowSurface + sfc - sfccore

      ! Calculate the effective radius (total volume / total area).
      ! NOTE: cm -> m.
      if (present(iceRe)) then
        where (pa(:) > 0.0_r8)
          iceRe(:) = (3._f / 4._f) * (md(:) / (0.917_f * pa(:))) * 1e-2_f
        end where
        
      end if
    end do

    !+++ Cheng: 2025/7/25. Repartion rain and snow from Graupel based on T.
    call CARMASTATE_GetState(cstate, rc, t=t)
    !do k = 1, pver
    !if (t(k) < 273.15_f) then
        !!if (do_print) 
    !    write(LUNOPRT,*) "for snowSurface from Graupel: T[K] =", t(k), "vertical_index_k =", k
    !endif
    !enddo
 
    !+++ Note_1: Add "if camra_do_graupel then" in the future
    !+++ Note_2: Noted by Cheng: adding surface grapel to snow is OK. But adding here may induce mass inbalance when using diagnosis.
    !+++ surface(Snow) > surface(Dice and Sice), while graupel mass is calculated in a seperate subroutine: CARMA_GetTotalGraupel. 
    !!--- Zhu Sep 12 2024 add graupel surface precipitation to ice precipitation -----
    igroup = I_GRP_CRGRP
    ielem  = I_ELEM_CRGRP
    
    if (t(pver) < 273.15_f + 1.0_f) then   !+++ surface-lev T
        do ibin = 1, NBIN
            call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, surface=sfc)
            !+++ Cheng: fixed the bug for graupel precipitation
            snowSurface = snowSurface + sfc
        end do
        ! !+++ just the last bin for graupel
        !snowSurface = snowSurface + sfc 
    endif
    !!--------------------------------------------------------------------------------

    return
  end subroutine CARMA_GetTotalIceAndSnow
    
    
  !! Subroutine Carma_GetTotalWaterAndRain 
  !!
  !! Determine the total cloud water concentration and number stored in the bins that represent
  !! water within the CARMA model.
  !!
  !!  @version Nov-2009 
  !!  @author  Chuck Bardeen 
  !! --------------------------------------------------------------------------
  subroutine CARMA_GetTotalWaterAndRain (carma, cstate, waterMass, waterNumber, rainSurface, rc, liqRe )    

    implicit none
    
    type(carma_type), intent (in)          :: carma   !! the carma object
    type(carmastate_type), intent (inout)  :: cstate  !! the carma state object
    real(kind=f), intent(out)              :: waterMass(pver)  !! water mass mixing ratio (kg/kg)
    real(kind=f), intent(out)              :: waterNumber(pver)  !! water number mixing ratio (#/kg)
    real(kind=f), intent(out)              :: rainSurface !! rain on surface (kg/m2)
    integer, intent(out)                   :: rc  !! return code, negative indicates failure
    real(kind = f), intent(out), optional  :: liqRe(pver)   !! Effective radius (m) of liquid size distribution
    
    integer             :: igroup            ! group index
    integer             :: ielem             ! element index
    integer             :: ibin              ! bin index
  
    real(kind=f)        :: mmr(pver)         ! mass mixing ratio (kg/kg)
    real(kind=f)        :: nmr(pver)         ! number mixing ratio (#/kg)
    real(kind=f)        :: sfc               ! surface mass (kg/m2)

    !! Summations that are intermediate to calculating the effective radius 
    real (kind = f) :: SumRadCub (pver)  
    real (kind = f) :: SumRadSqu (pver)  
 
    !! Need the radius grid to calculate the effective radius, cm
    real (kind = f) :: r(NBin)  
    
    !+++
    real(kind=f) :: t(pver)             ! temperature (K)
    !integer      :: k             ! vertical index
  
  
    rc = RC_OK

    WaterMass    = 0.0_f
    WaterNumber  = 0.0_f
    RainSurface  = 0.0_f
    
    !  JAS, What is a good default effective radius?

    if (present (LiqRe)) LiqRe = 14.0e-6_f  ! 14 micrometers
    SumRadCub = 0.0_f 
    SumRadSqu = 0.0_f 
    r = 1.0_f  ! 1.0 is unlikely to produce an FPE
      
    igroup = I_GRP_CRLIQ
    ielem  = I_ELEM_CRLIQ

    call CARMAGROUP_Get(carma, igroup, rc, r=r)
    if (rc < RC_OK) call endrun ( &
     'carma_model_mod.F90, Carma_GetTotalWaterAndRain :: ' &
     // 'CARMAGROUP_Get failed.' )

    do ibin = 1, NBIN

      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, nmr=nmr, surface=sfc)
      if (rc < RC_OK) call endrun ( &
       'carma_model_mod.F90, Carma_GetTotalWaterAndRain :: ' &
       // 'CARMASTATE_GetBin failed.' )

      waterMass(:)   = waterMass(:)   + mmr(:)
      waterNumber(:) = waterNumber(:) + nmr(:)
      
      rainSurface = rainSurface + sfc

      ! Include the detrained liquid that hasn't been added to the particle
      ! bins yet.
      call CARMASTATE_GetDetrain(cstate, ielem, ibin, mmr, rc, nmr=nmr)
      if (rc < RC_OK) call endrun('CARMA_GetTotalWaterAndRain::CARMASTATE_GetDetrain failed.')
      
      waterMass   = waterMass   + mmr
      waterNumber = waterNumber + nmr

      !  JAS, Calculate the effective radius, Sum over bins of [ r**3 N(r) ] / 
      ! Sum over bins of [ r**2 N(r) ]. Convert cm to m.

      !  JAS, Is r(ibin) known here? It is now that I add the call to 
      ! CarmaGroup_Get. Need to be mindful of units, r [=] cm, LiqRe [=] m.
      if (present(liqRe)) then
        SumRadCub = SumRadCub + (r(ibin)**3) * nmr
        SumRadSqu = SumRadSqu + (r(ibin)**2) * nmr
      end if
    end do  ! ibin

    if (present(liqRe)) then
      where (SumRadSqu .gt. tiny(SumRadSqu))
        liqRe = SumRadCub / SumRadSqu / 100.0_f  ! in meters
      end where
    end if
    
    !+++ Cheng: 2025/7/25. Repartion rain and snow from Graupel based on T.
    call CARMASTATE_GetState(cstate, rc, t=t)
    
    igroup = I_GRP_CRGRP
    ielem  = I_ELEM_CRGRP
    
    if (t(pver) >= 273.15_f + 1.0_f) then   !+++ surface-lev T
        do ibin = 1, NBIN
            call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, surface=sfc)
            rainSurface = rainSurface + sfc 
        end do
    endif
    !+++ --------------------------------------------------------------------------------
    
    return
  end subroutine CARMA_GetTotalWaterAndRain

  !! +++ 
  !! Subroutine CARMA_GetCloudlWater(AndRain) 
  !!
  !! Determine the cloud water (0-50 µm) concentration and number stored in the bins that represent
  !! water within the CARMA model.
  !! We want to pass the cloud water (0-50 µm) to CLUBB with new CLDLIQ, instead of total water.
  !!
  !!  @version Feb-2025 
  !!  @author  Cheng-Cheng Liu 
  !! --------------------------------------------------------------------------
  subroutine CARMA_GetCloudWaterAndRain (carma, cstate, CloudWaterMass, CloudWaterNumber, rainSurface, rc, liqRe)    

    implicit none
    
    type(carma_type), intent (in)          :: carma   !! the carma object
    type(carmastate_type), intent (inout)  :: cstate  !! the carma state object
    real(kind=f), intent(out)              :: CloudWaterMass(pver)  !! water mass mixing ratio (kg/kg)
    real(kind=f), intent(out)              :: CloudWaterNumber(pver)  !! water number mixing ratio (#/kg)
    real(kind=f), intent(out)              :: rainSurface !! rain on surface (kg/m2)
    integer, intent(out)                   :: rc  !! return code, negative indicates failure
    real(kind = f), intent(out), optional  :: liqRe(pver)   !! Effective radius (m) of liquid size distribution
    
    integer             :: igroup            ! group index
    integer             :: ielem             ! element index
    integer             :: ibin              ! bin index
  
    real(kind=f)        :: mmr(pver)         ! mass mixing ratio (kg/kg)
    real(kind=f)        :: nmr(pver)         ! number mixing ratio (#/kg)
    real(kind=f)        :: sfc               ! surface mass (kg/m2)

    !! Summations that are intermediate to calculating the effective radius 
    real (kind = f) :: SumRadCub (pver)  
    real (kind = f) :: SumRadSqu (pver)  
 
    !! Need the radius grid to calculate the effective radius, cm
    real (kind = f) :: r(NBin)  

    ! +++
    real(r8)      :: r_wet(pver)
    real(kind=f), parameter :: r_limit_forCloudWater  =  50e-4_f   !! wet radius for cloud water < 50 µm ([µm] to [cm])
    
    !+++
    real(kind=f) :: t(pver)             ! temperature (K)
    !integer      :: k             ! vertical index
  
    rc = RC_OK

    CloudWaterMass    = 0.0_f
    CloudWaterNumber  = 0.0_f
    RainSurface  = 0.0_f
    
    !  JAS, What is a good default effective radius?

    if (present (LiqRe)) LiqRe = 14.0e-6_f  ! 14 micrometers
    SumRadCub = 0.0_f 
    SumRadSqu = 0.0_f 
    r = 1.0_f  ! 1.0 is unlikely to produce an FPE
      
    igroup = I_GRP_CRLIQ
    ielem  = I_ELEM_CRLIQ

    call CARMAGROUP_Get(carma, igroup, rc, r=r)
    if (rc < RC_OK) call endrun ( &
     'carma_model_mod.F90, Carma_GetTotalWaterAndRain :: ' &
     // 'CARMAGROUP_Get failed.' )

    do ibin = 1, NBIN

      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, nmr=nmr, surface=sfc, r_wet=r_wet)
      if (rc < RC_OK) call endrun ( &
       'carma_model_mod.F90, Carma_GetCloudlWaterAndRain :: ' &
       // 'CARMASTATE_GetBin failed.' )

      ! if (r(ibin) .le. r_limit_forCloudWater) then
      where(r_wet(:) .le. (r_limit_forCloudWater))
        CloudWaterMass(:)   = CloudWaterMass(:)   + mmr(:)
        CloudWaterNumber(:) = CloudWaterNumber(:) + nmr(:)
      end where
      ! end if
      
      rainSurface = rainSurface + sfc

      ! Include the detrained liquid that hasn't been added to the particle
      ! bins yet.
      call CARMASTATE_GetDetrain(cstate, ielem, ibin, mmr, rc, nmr=nmr)
      if (rc < RC_OK) call endrun('CARMA_GetTotalWaterAndRain::CARMASTATE_GetDetrain failed.')
      where(r_wet(:) .le. (r_limit_forCloudWater))
        CloudWaterMass   = CloudWaterMass   + mmr
        CloudWaterNumber = CloudWaterNumber + nmr
      end where

      !  JAS, Calculate the effective radius, Sum over bins of [ r**3 N(r) ] / 
      ! Sum over bins of [ r**2 N(r) ]. Convert cm to m.

      !  JAS, Is r(ibin) known here? It is now that I add the call to 
      ! CarmaGroup_Get. Need to be mindful of units, r [=] cm, LiqRe [=] m.
      if (present(liqRe)) then
        SumRadCub = SumRadCub + (r(ibin)**3) * nmr
        SumRadSqu = SumRadSqu + (r(ibin)**2) * nmr
      end if
    end do  ! ibin

    if (present(liqRe)) then
      where (SumRadSqu .gt. tiny(SumRadSqu))
        liqRe = SumRadCub / SumRadSqu / 100.0_f  ! in meters
      end where
    end if

    !+++ Cheng: 2025/7/25. Repartion rain and snow from Graupel based on T.
    call CARMASTATE_GetState(cstate, rc, t=t)
    
    igroup = I_GRP_CRGRP
    ielem  = I_ELEM_CRGRP
    
    if (t(pver) >= 273.15_f + 1.0_f) then   !+++ surface-lev T
        do ibin = 1, NBIN
            call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, surface=sfc)
            rainSurface = rainSurface + sfc 
        end do
    endif
    !+++ --------------------------------------------------------------------------------
    
    return
  end subroutine CARMA_GetCloudWaterAndRain

  !!---------------------- zhu Sep 12 2024 ----------------------------------
  subroutine CARMA_GetTotalGraupel(carma, cstate, graupelicmass, graupelnum, graupelde, rc)
    implicit none
    
    type(carma_type), intent(in)         :: carma            !! the carma object                       
    type(carmastate_type), intent(inout) :: cstate           !! the carma state object                 
    real(kind=f), intent(out)              :: graupelicmass(pver)  !! graupel in-cloud mass mixing ratio (kg/kg)
    real(kind=f), intent(out)              :: graupelnum(pver)  !! graupel number mixing ratio (#/kg)
    real(kind=f), intent(out)              :: graupelde(pver) !! graupel effective diameter (m)
    integer, intent(out)                   :: rc  !! return code, negative indicates failure

    integer             :: igroup            ! group index
    integer             :: ielem             ! element index
    integer             :: ibin              ! bin index

    real(kind=f)        :: mmr(pver)         ! mass mixing ratio (kg/kg)
    real(kind=f)        :: nmr(pver)         ! number mixing ratio (#/kg)
    real(kind=f)        :: cldfrc(pver)      ! cloud fraction
    real(kind=f)        :: sfc               ! surface mass (kg/m2)

    !! Summations that are intermediate to calculating the effective diameter
    real (kind = f) :: SumRadCub (pver)  
    real (kind = f) :: SumRadSqu (pver)  
    
    !! Need the radius grid to calculate the effective diameter, cm
    real (kind = f) :: r(NBin)  

    rc = RC_OK

    graupelicmass = 0.0_f
    graupelnum = 0.0_f
    graupelde = 0.0_f

    igroup = I_GRP_CRGRP
    ielem  = I_ELEM_CRGRP

    call CARMAGROUP_Get(carma, igroup, rc, r=r)
    if (rc < RC_OK) call endrun('carma_model_mod.F90, CARMA_GetTotalGraupel: CARMAGROUP_Get failed.')   

    call CARMASTATE_Get(cstate, rc, cldfrc=cldfrc)
    do ibin = 1, NBIN
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, nmr=nmr, surface=sfc)
      if (rc < RC_OK) call endrun('carma_model_mod.F90, CARMA_GetTotalGraupel: CARMASTATE_GetBin failed.')
    
      graupelicmass(:) = graupelicmass(:) + mmr(:)/max(cldfrc(:),1e-4_f)
      graupelnum(:) = graupelnum(:) + nmr(:)

      SumRadCub = SumRadCub + (r(ibin)**3) * nmr
      SumRadSqu = SumRadSqu + (r(ibin)**2) * nmr
    end do

    where (SumRadSqu .gt. tiny(SumRadSqu))
      graupelde = SumRadCub / SumRadSqu / 100.0_f*2._f  ! in meters
    end where

    return
  end subroutine CARMA_GetTotalGraupel
  !!-------------------------------------------------------------

  subroutine CARMA_CheckMassAndEnergy(carma, cstate, name, state, &
       icol, dt, dlf, waterMass, rainSurface, iceMass, snowSurface, rc)
    implicit none
    
    type(carma_type), intent(in)         :: carma            !! the carma object
    type(carmastate_type), intent(inout) :: cstate           !! the carma state object
    character*(*),intent(in)             :: name             !! test name
    type(physics_state), intent(in)      :: state            !! physics state variables
    integer, intent(in)                  :: icol             !! column index
    real(kind=f), intent(in)             :: dt               !! time step
    real(kind=f), intent(in)             :: dlf(pcols, pver) !! detrainment rate (kg/kg/s)
    real(kind=f), intent(in)             :: waterMass(pver)  !! water mass mixing ratio (kg/kg)
    real(kind=f), intent(in)             :: rainSurface      !! rain mass at surface (kg/m2)
    real(kind=f), intent(in)             :: iceMass(pver)    !! ice mass mixing ratio (kg/kg)
    real(kind=f), intent(in)             :: snowSurface      !! snow mass at surface (kg/m2)
    integer, intent(out)                 :: rc               !! return code, negative indicates failure
 
 
    integer                              :: LUNOPRT              ! logical unit number for output
    logical                              :: do_print             ! do print output?
    logical                              :: do_detrain           ! do convective detrainment?

    real(kind=f)                         :: mmr(pver)          ! mass mixing ratio (#/kg)
    real(kind=f)                         :: totalMass
    real(kind=f)                         :: totalMass2

    real(r8)                             :: lat
    real(r8)                             :: lon
    

  1 format(/,'CARMA_CheckMassAndEnergy::ERROR - CARMA mass conservation error, ',a,',icol=',i4,',lat=',&
              f7.2,',lon=',f7.2,',cam=',e16.10,',carma=',e16.10,',diff=',e16.10,',rer=',e9.3)

    ! Default return code.
    rc = RC_OK
    
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT, do_detrain=do_detrain)
    if (rc < RC_OK) call endrun('CARMA_CheckMassAndEnergy::CARMA_Get failed.') 

    ! Get the total mass that came in from CAM     
    totalMass = sum(state%q(icol, :, ixcldliq) * (state%pdel(icol, :) / gravit))
    totalMass = totalMass + sum(state%q(icol, :, ixcldice) * (state%pdel(icol, :) / gravit))
    totalMass = totalMass + sum(state%q(icol, :, 1) * (state%pdel(icol, :) / gravit))
    
    ! Need to also add the rain now that it is prognostic. This is just the rain in
    ! the atmosphere, not the rain at the surface.
    totalMass = totalMass + sum(state%q(icol, :, ixrainqm) * (state%pdel(icol, :) / gravit))
    totalMass = totalMass + sum(state%q(icol, :, ixsnowqm) * (state%pdel(icol, :) / gravit))
    
    if (abs((totalMass - state%tw_cur(icol))) / state%tw_cur(icol) > 1e14_f) then
      if (do_print) then
         write(LUNOPRT,*) "CARMA_CheckMassAndEnergy::&
              &WARNING Total water not conserved, ", &
              totalMass, state%tw_cur, (totalMass - state%tw_cur(icol)), &
              (totalMass - state%tw_cur(icol)) / state%tw_cur(icol)
      end if        
    end if

    ! dlf contains mass that is detraining into the atmosphere, but has not
    ! yet been included in CLDLIQ or CLDICE.
    !
    ! NOTE: CLUBB has already detrained these into cldliq and cldice, but since
    ! we put the CMELIQ and CMEICE back into vapor for CARMA when need to compare
    ! to state without the detrained liquid. However, we need to add it back to
    ! compare with CARMA which 
    if (do_detrain) then
      totalMass = totalMass + sum(dlf(icol, :) * (state%pdel(icol, :) / gravit)) * dt
    end if
    
    
    ! Get the total water coming out of CARMA
    call CARMASTATE_GetGas(cstate, I_GAS_H2O, mmr(:), rc)
    if (rc < RC_OK) call endrun('CARMA_CheckMassAndEnergy::CARMASTATE_GetGas failed.')

    totalMass2 = sum(mmr(:) * (state%pdel(icol, :) / gravit))

    totalMass2 = totalMass2 + sum(waterMass(:) * (state%pdel(icol, :) / gravit))
    totalMass2 = totalMass2 + sum((iceMass(:)) * (state%pdel(icol, :) / gravit))
    
    totalMass2 = totalMass2 + rainSurface
    totalMass2 = totalMass2 + snowSurface


    if (totalMass /= totalMass2) then

      if (totalMass /= 0._f) then

        if (abs((totalMass - totalMass2) / totalMass) > 1e-10_f)  then
          if (do_print) then
            call CARMASTATE_Get(cstate, rc, lat=lat, lon=lon)
            if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_Get failed.')

            write(LUNOPRT,1) name, icol, lat, lon, totalMass, totalMass2, &
                 totalMass2-TotalMass, (totalMass - totalMass2) / totalMass

            write(LUNOPRT,*) "  state tw :  ", state%tw_cur(icol)
            write(LUNOPRT,*) ""
            write(LUNOPRT,*) "  old vap  :  ", sum(state%q(icol, :, 1) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  old liq  :  ", sum(state%q(icol, :, ixcldliq) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  old ice  :  ", sum(state%q(icol, :, ixcldice) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  old rain :  ", sum(state%q(icol, :, ixrainqm) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  old snow :  ", sum(state%q(icol, :, ixsnowqm) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  dlf      :  ", sum(dlf(icol, :) * dt * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) ""
            write(LUNOPRT,*) "  new vap  :  ", sum(mmr(:) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  new liq  :  ", sum(waterMass(:) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  new ice  :  ", sum(iceMass(:) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  new rain :  ", rainSurface
            write(LUNOPRT,*) "  new snow :  ", snowSurface
           end if
        end if
      end if
    end if

    return
  end subroutine CARMA_CheckMassAndEnergy

end module
