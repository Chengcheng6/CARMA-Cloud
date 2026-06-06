! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! @author Chuck Bardeen
!! @version Jan-2010
subroutine immersenucl(carma, cstate, iz, rc)

  ! types
  use physconst,   only: mwh2o,rhoh2o,tmelt
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

   type(carma_type), intent(in)         :: carma   !! the carma object
   type(carmastate_type), intent(inout) :: cstate  !! the carma state object
   integer, intent(in)                  :: iz      !! vertical index
   integer, intent(inout)               :: rc      !! return code, negative indicates failure

   ! local variables
   real(kind=f)     :: aw(3)                           ! water activity [ ]
   real(kind=f)     :: molal(3)                        ! molality [moles/kg]
   real(kind=f)     :: vwice
   real(kind=f)     :: rho                             ! air density (kg m-3)
   real(kind=f)     :: con1,r3lx,qcic,ncic
   real(kind=f)     :: rgimm,rgimm_bc,rgimm_dust_a1,rgimm_dust_a3
   real(kind=f)     :: sigma_iw,tc,supersatice,rhoice
   real(kind=f)     :: dg0imm_bc, dg0imm_dust_a1, dg0imm_dust_a3
   real(kind=f)     :: Aimm_bc, Aimm_dust_a1, Aimm_dust_a3
   real(kind=f)     :: Jimm_bc,Jimm_dust_a1,Jimm_dust_a3
   real(kind=f)     :: m
   real(kind=f)     :: f_imm_bc, f_imm_dust_a1, f_imm_dust_a3
   real(kind=f)     :: bcimm_pct,duimm_pct1,duimm_pct2
   real(kind=f)     :: pvap_liq,pvap_ice

   real(kind=f), parameter :: amu = 1.66053886e-27_f
   real(kind=f), parameter :: Mso4 = 96.06_f
   real(kind=f), parameter :: mincld = 0.0001_f
   real(kind=f), parameter :: kboltz = 1.38e-23_f  
   real(kind=f), parameter :: n1 = 1.e19_f     ! number of water molecules in contact with unit area of substrate [m-2]
   real(kind=f), parameter :: hplanck = 6.63e-34_f
   real(kind=f), parameter :: rhplanck = 1._f/hplanck

   logical :: do_bc, do_dst1, do_dst3

   !  Local declarations
  integer                              :: igroup  !! group index
  integer                              :: ibin    !! bin index
  integer                              :: iepart  !! element for condensing group index
  integer                              :: inuc    !! nucleating element index
  integer                              :: ienucto !! index of target nucleation element
  integer                              :: ignucto !! index of target nucleation group


   !********************************************************
   ! Wang et al., 2014 fitting parameters
   !********************************************************
   ! freezing parameters for immersion freezing
   real(kind=f),parameter :: theta_imm_bc = 48.0_f            ! contact angle [deg], converted to rad later !DeMott et al (1990)
   real(kind=f),parameter :: dga_imm_bc = 14.15E-20_f         ! activation energy [J]
   real(kind=f),parameter :: theta_imm_dust = 46.0_f          ! contact angle [deg], converted to rad later !DeMott et al (2011) SD
   real(kind=f),parameter :: dga_imm_dust = 14.75E-20_f       ! activation energy [J]
   ! freezing parameters for deposition nucleation
   real(kind=f),parameter :: theta_dep_dust = 20.0_f          ! contact angle [deg], converted to rad later !Koehler et al (2010) SD
   real(kind=f),parameter :: dga_dep_dust = -8.1E-21_f        ! activation energy [J]
   real(kind=f),parameter :: theta_dep_bc = 28._f             ! contact angle [deg], converted to rad later !Moehler et al (2005), soot
   real(kind=f),parameter :: dga_dep_bc = -2.E-19_f           ! activation energy [J]

   rho = rhoa_wet(iz)*1.e3_f !kg/m3

   tc = t(iz) - tmelt
   sigma_iw = (28.5_f+0.25_f*tc)*1E-3_f
   rhoice = 916.7_f-0.175_f*tc-5.e-4_f*tc**2
   vwice = mwh2o*amu/rhoice

   if (t(iz) > 235.15_f .and. t(iz) < tmelt) then
! NOTE: Should remove lcldm from the sate and just use cldfrc.
     !lcldm = max(ast(iz), mincld)          ! Stratiform cloud fraction

     ! NOTE: this is being done for the liquid water group. It should not be hard coded
     ! like this, but should be driven as an attribute of the group. Then you iterate
     ! over all groups and do it for the ones that are configured for it. Perhaps this
     ! al belongs inside the loop at the bottom. Are these values used anywhere else?
     !
     ! NOTE: For now, just patching the bug in rmass
!     qcic = min(sum(pc(iz,:,2))*rmass(iz,2)/rhoa_wet(iz)/lcldm(iz), 5.e-3_f) !kg/kg water drop
!     ncic = max(sum(pc(iz,:,2))/rhoa_wet(iz)/lcldm(iz), 0._f)          !#/kg water drop
     qcic = min(sum(pc(iz,:,ienconc(2))*rmass(:,2))/rhoa_wet(iz)/lcldm(iz), 5.e-3_f) !kg/kg water drop
     ncic = max(sum(pc(iz,:,ienconc(2)))/rhoa_wet(iz)/lcldm(iz), 0._f)          !#/kg water drop

     con1 = 1._f/(1.333_f*pi)**0.333_f
     r3lx = con1*(rho*qcic/(rhoh2o*max(ncic*rho, 1.0e6_f)))**0.333_f ! in m
     r3lx = max(4.e-6_f, r3lx)

   !*****************************************************************************
   !                take water activity into account 
   !*****************************************************************************
   !   solute effect
   aw(:) = 1._f
   molal(:) = 0._f

   ! The heterogeneous ice freezing temperatures of all IN generally decrease with
   ! increasing total solute mole fraction. Therefore, the large solution concentration
   ! will cause the freezing point depression and the ice freezing temperatures of all
   ! IN will get close to the homogeneous ice freezing temperatures. Since we take into
   ! account water activity for three heterogeneous freezing modes(immersion, deposition, 
   ! and contact), we utilize interstitial aerosols(not cloudborne aerosols) to calculate 
   ! water activity. 
   ! If the index of IN is 0, it means three freezing modes of this aerosol are depressed.

   !calculate molality
   if ( immersebcit(iz) > 0._f ) then
     molal(1) = (1.e-6_f*awcambc(iz)*(1._f-awfacmbc(iz))/(Mso4*immersebcit(iz)*1.e6_f))/ &
            (4*pi/3*rhoh2o*(MAX(r3lx,4.e-6_f))**3)
     aw(1) = 1._f/(1._f+2.9244948e-2_f*molal(1)+2.3141243e-3_f*molal(1)**2 &
            +7.8184854e-7_f*molal(1)**3)
   end if

   if ( immersed1it(iz) > 0._f ) then
     molal(2) = (1.e-6_f*awcamd1(iz)*(1._f-awfacmd1(iz))/(Mso4*immersed1it(iz)*1.e6_f))/ &
            (4*pi/3*rhoh2o*(MAX(r3lx,4.e-6_f))**3)
     aw(2) = 1._f/(1._f+2.9244948e-2_f*molal(2)+2.3141243e-3_f*molal(2)**2 &
            +7.8184854e-7_f*molal(2)**3)
   end if

   if ( immersed2it(iz) > 0._f ) then
    molal(3) = (1.e-6_f*awcamd2(iz)*(1._f-awfacmd2(iz))/(Mso4*immersed2it(iz)*1.e6_f))/ &
            (4*pi/3*rhoh2o*(MAX(r3lx,4.e-6_f))**3)
     aw(3) = 1._f/(1._f+2.9244948e-2_f*molal(3)+2.3141243e-3_f*molal(3)**2 &
            +7.8184854e-7_f*molal(3)**3)
   end if

   !*****************************************************************************
   !                immersion freezing begin 
   !*****************************************************************************    
   pvap_liq = 10.0_f * exp(54.842763_f - (6763.22_f / t(iz)) - &
             (4.210_f * log(t(iz))) + (0.000367_f * t(iz)) + &
             (tanh(0.0415_f * (t(iz) - 218.8_f)) * &
             (53.878_f - (1331.22_f / t(iz)) - &
             (9.44523_f * log(t(iz))) + 0.014025_f * t(iz))))
   pvap_ice = 10.0_f * exp(9.550426_f - (5723.265_f / t(iz)) + &
             (3.53068_f * log(t(iz))) - (0.00728332_f * t(iz)))
   supersatice = pvap_liq/pvap_ice

   ! critical germ size
   rgimm = 2.*vwice*sigma_iw/(kboltz*t(iz)*log(supersatice))
   ! take solute effect into account
   rgimm_bc = rgimm
   rgimm_dust_a1 = rgimm
   rgimm_dust_a3 = rgimm

   ! if aw*Si<=1, the freezing point depression is strong enough to prevent freezing
   if (aw(1)*supersatice > 1._f ) then
      do_bc   = .true.
      rgimm_bc = 2*vwice*sigma_iw/(kboltz*t(iz)*log(aw(1)*supersatice))
      !write(*,*) "aw(1)",aw(1),"supersatice",supersatice
   else
      do_bc = .false.
   end if

   if (aw(2)*supersatice > 1._f ) then
      do_dst1 = .true.
      rgimm_dust_a1 = 2*vwice*sigma_iw/(kboltz*t(iz)*log(aw(2)*supersatice))
      !write(*,*) "aw(2)",aw(2),"supersatice",supersatice
   else
      do_dst1 = .false.
   end if

   if (aw(3)*supersatice > 1._f ) then
      do_dst3 = .true.
      rgimm_dust_a3 = 2*vwice*sigma_iw/(kboltz*t(iz)*log(aw(3)*supersatice))
      !write(*,*) "aw(3)",aw(3),"supersatice",supersatice
   else
      do_dst3 = .false.
   end if

   ! form factor
   ! only consider flat surfaces due to uncertainty of curved surfaces

   m = cos(theta_imm_bc*pi/180._f)
   f_imm_bc = (2+m)*(1-m)**2/4._f
   !if (.not. pdf_imm_in) then
      m = cos(theta_imm_dust*pi/180._f)
      f_imm_dust_a1 = (2+m)*(1-m)**2/4._f

      m = cos(theta_imm_dust*pi/180._f)
      f_imm_dust_a3 = (2+m)*(1-m)**2/4._f
   !end if

   ! homogeneous energy of germ formation
   dg0imm_bc = 4*pi/3._f*sigma_iw*rgimm_bc**2
   dg0imm_dust_a1 = 4*pi/3._f*sigma_iw*rgimm_dust_a1**2
   dg0imm_dust_a3 = 4*pi/3._f*sigma_iw*rgimm_dust_a3**2

   ! prefactor
   Aimm_bc = n1*((vwice*rhplanck)/(rgimm_bc**3)*sqrt(3._f/pi*kboltz*t(iz)*dg0imm_bc))
   Aimm_dust_a1 = n1*((vwice*rhplanck)/(rgimm_dust_a1**3)*sqrt(3._f/pi*kboltz*t(iz)*dg0imm_dust_a1))
   Aimm_dust_a3 = n1*((vwice*rhplanck)/(rgimm_dust_a3**3)*sqrt(3._f/pi*kboltz*t(iz)*dg0imm_dust_a3))

   ! nucleation rate per particle
   Jimm_bc = 0._f
   Jimm_dust_a1 = 0._f
   Jimm_dust_a3 = 0._f

   if(do_bc) Jimm_bc = Aimm_bc*hetrbc(iz)**2/sqrt(f_imm_bc)*exp((-dga_imm_bc-f_imm_bc*dg0imm_bc)/(kboltz*t(iz)))
   !if (.not. pdf_imm_in) then
      ! 1/sqrt(f)
      ! the expression of Chen et al. (sqrt(f)) may however lead to unphysical
      ! behavior as it implies J->0 when f->0 (i.e. ice nucleation would be
      ! more difficult on easily wettable materials). 
      if(do_dst1) Jimm_dust_a1 = Aimm_dust_a1*hetrd1(iz)**2/sqrt(f_imm_dust_a1)*exp((-dga_imm_dust-f_imm_dust_a1*dg0imm_dust_a1)/(kboltz*t(iz)))
      if(do_dst3) Jimm_dust_a3 = Aimm_dust_a3*hetrd2(iz)**2/sqrt(f_imm_dust_a3)*exp((-dga_imm_dust-f_imm_dust_a3*dg0imm_dust_a3)/(kboltz*t(iz)))
   !end if 

   ! percent
   bcimm_pct = immersebccb(iz)/(immersebccb(iz) + immersed1cb(iz) + immersed2cb(iz))
   duimm_pct1 = immersed1cb(iz)/(immersebccb(iz) + immersed1cb(iz) + immersed2cb(iz))
   duimm_pct2 = immersed2cb(iz)/(immersebccb(iz) + immersed1cb(iz) + immersed2cb(iz))

   !imm_prob = min(.99_f,(1._f-exp(-Jimm_bc*dtime))*bcimm_pct(iz)+(1._f-exp(-Jimm_dust_a1*dtime))*duimm_pct1+(1._f-exp(-Jimm_dust_a3*dtime))*distimm_pct2)/dtime
  else
    Jimm_bc =0._f
    Jimm_dust_a1 = 0._f
    Jimm_dust_a3 = 0._f

    bcimm_pct = 0._f
    duimm_pct1 = 0._f
    duimm_pct2 = 0._f
  end if

     ! Loop over particle groups.
  do igroup = 1,NGROUP
  
    iepart = ienconc( igroup )            ! particle number density element
  
    ! Calculate nucleation loss rates.
    do inuc = 1,nnuc2elem(iepart)
    
      ienucto = inuc2elem(inuc,iepart)
      
      if( ienucto .ne. 0 )then
        ignucto = igelem( ienucto )
  
        ! Only compute nucleation rate for droplet freezing
        if( inucproc(iepart,ienucto) .eq. I_IMMNUCL ) then
    
          ! Loop over particle bins.  
          do ibin = 1,NBIN
      
            ! Bypass calculation if few particles are present 
            if( pc(iz,ibin,iepart) .gt. FEW_PC )then
              !rnuclg(ibin,igroup,ignucto) = min(.99_f,min(0.01_f,(1._f-exp(-Jimm_bc*dtime)))*bcimm_pct+(1._f-exp(-Jimm_dust_a1*dtime))*duimm_pct1+(1._f-exp(-Jimm_dust_a3*dtime))*duimm_pct2)/dtime 
              rnuclg(ibin,igroup,ignucto) = rnuclg(ibin,igroup,ignucto) + min(.99_f,(1._f-exp(-Jimm_bc*dtime))*bcimm_pct+(1._f-exp(-Jimm_dust_a1*dtime))*duimm_pct1+(1._f-exp(-Jimm_dust_a3*dtime))*duimm_pct2)/dtime
              !write(*,*) "ibin",ibin,"igroup",igroup,"ignucto",ignucto,"rnuclg",rnuclg(ibin,igroup,ignucto)
              !write(*,*) "t(iz)",t(iz),"supersatice",supersatice,"aw(1)",aw(1),"aw(2)",aw(2),"aw(3)",aw(3)
              !write(*,*) "molal(1)",molal(1),"molal(2)",molal(2),"molal(3)",molal(3)
              !write(*,*) "awcambc(iz)",awcambc(iz),"awfacmbc(iz)",awfacmbc(iz),"immersebcit(iz)",immersebcit(iz),"r3lx",r3lx
              !write(*,*) "Jimm_bc",Jimm_bc,"dtime",dtime,"bcimm_pct",bcimm_pct,"Jimm_dust_a1",Jimm_dust_a1,"duimm_pct1",duimm_pct1,"Jimm_dust_a3",Jimm_dust_a3,"duimm_pct2",duimm_pct2
              !write(*,*) "Aimm_dust_a3",Aimm_dust_a3,"hetrd2(iz)",hetrd2(iz),"f_imm_dust_a3",f_imm_dust_a3,"dga_imm_dust",dga_imm_dust,"dg0imm_dust_a3",dg0imm_dust_a3,"kboltz",kboltz,"t(iz)",t(iz)
              !rnuclg(ibin,igroup,ignucto) = rnuclg(ibin,igroup,ignucto)/100.
            endif     ! pc(source particles) .gt. FEW_PC
          enddo      ! ibin = 1,NBIN
        endif       ! inucproc(iepart,ienucto) .eq. I_IMMNUCL
      endif
    enddo        ! inuc = 1,nnuc2elem(iepart)
  enddo         ! igroup = 1,NGROUP

  return
end 
