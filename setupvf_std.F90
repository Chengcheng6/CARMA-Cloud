! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle fall velocities, vf(k) [cm s^-1]
!! and reynolds' numbers (re) based on fall velocities, re(j,i,k) [dimensionless].
!! indices correspond to vertical level <k>, bin index <i>, and aerosol
!! group <j>.
!! The following velocities are designed for liquid water droplets (and probably other liquid droplets) and are
!! based on Beard(1976, JAS,33,851-864) 
!! Table 1 These values seem compatible with more recent studies such as Bohm (1992 Atmospheric Res. 27, 253-274)
!! and Lamb and Verlinde Physics and Chemistry of Clouds (2011) 
!! Method: first use Stokes Cunningham flow (with Fuchs' size corrections, 
!! valid only for Stokes flow) to estimate fall velocity, for Reynolds number less than 0.01,which should correspond
!! to  for particle radius less than 15 µm (Beard switched at 
!! 9.5µµ, this choice is in Lamb and Verlinde).   These equations are analytic, and in the cloud size range v varies
!! as r^2.  Above 15µm (Re>0.01) the fall velocity comes from the Reynolds number definition.
!! For radius between 15µm and 0.5 mm (re from 0.01 to 300) correct drag coefficient (Cd) for turbulent boundary
!! layer through standard trick to solving the drag problem using the Davies (Best) number (NDa). The Davies number is
!! CdRe^2 it only depends on gravity, air density, particle density, particle radius and air viscosity.
!! Use data to fit y = log( Re ) as a function of x = log( Cd Re^2 ).  Vfall is proportional to r.
!! For r larger than 0.5 mm (re>300) the drops are no longer spheres.  The Bond number measures the ratio of gravity to 
!! surface tension so it measures the flattening.  The physical property number is a function of the Davies 
!! number and the Bond number.  It only depends on air density, surface tension, gravity, particle density and air
!! Viscosity.  The Reynolds number is found from an empirical function of the Bond number and physical property
!! number.
!! This routine requires that vertical profiles of temperature <t>,
!! air density <rhoa>, and viscosity <rmu> are defined (i.e., initatm.f
!! must be called before this).  The vertical profile with ix = iy = 1
!! is used.
!!
!! Added support for the particle radius being dependent on the relative
!! humidity according to the parameterizations of Gerber [1995] and
!! Fitzgerald [1975]. The fall velocity is then based upon the wet radius
!! rather than the dry radius. For particles that are not subject to
!! swelling, the wet and dry radii are the same. 
!!
!! authors  Yunqian Zhu, Brian Toon
!! @version May 2024 
subroutine setupvf_std(carma, cstate, j, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma    !! the carma object
  type(carmastate_type), intent(inout) :: cstate   !! the carma state object
  integer, intent(in)                  :: j        !! group index
  integer, intent(inout)               :: rc       !! return code, negative indicates failure

  ! Local declarations
  integer                 :: i, k
  real(kind=f)            :: x, y, cdrag
  real(kind=f)            :: rhoa_cgs, vg, rmfp, rkn, expon, NDa, NB, Npp,surfacetension
  !NDa is the Davies number, NB is the Bond number, Npp is the physical property number
                                   
  ! Define formats
  1 format(/,'Non-spherical particles specified for group ',i3, &
      ' (ishape=',i3,') but spheres assumed in I_FALLRTN_STD.', &
      ' Suggest using non-spherical code in I_FALLRTN_STD_SHAPE.')

  !  Warning message for non-spherical particles!
  if( ishape(j) .ne. 1 )then
    if (do_print) write(LUNOPRT,1) j, ishape(j)
  endif
  !write(*,*) "test new setupvf_std..."
  
  ! Loop over all atltitudes.
  do k = 1, NZ

    ! This is <rhoa> in cartesian coordinates (good old cgs units)
    rhoa_cgs = rhoa(k) / (xmet(k)*ymet(k)*zmet(k))

    ! <vg> is mean thermal velocity of air molecules [cm/s]
    vg = sqrt(8._f / PI * R_AIR * t(k))

    ! <rmfp> is mean free path of air molecules [cm]
    rmfp = 2._f * rmu(k) / (rhoa_cgs * vg)

    ! Loop over particle size bins.
    do i = 1,NBIN
    
      ! <rkn> is knudsen number
      rkn = rmfp / (r_wet(k,i,j) * rrat(i,j))

      ! <bpm> is the slip correction factor, the correction term for
      ! non-continuum effects.  Also used to calculate coagulation kernels
      ! and diffusion coefficients.
      expon = -.87_f / rkn
      expon = max(-POWMAX, expon)
      bpm(k,i,j) = 1._f + (1.246_f*rkn + 0.42_f*rkn*exp(expon))

      ! Stokes fall velocity and Reynolds' number
      vf(k,i,j) = (ONE * 2._f / 9._f) * rhop_wet(k,i,j) * r_wet(k,i,j)**2 * GRAV * bpm(k,i,j) / rmu(k) / rprat(i,j)
      re(k,i,j) = 2. * rhoa_cgs * r_wet(k,i,j) * rprat(i,j) * vf(k,i,j) / rmu(k)

      !write(*,*) "re",re(k,i,j),"r_wet(k,i,j)",r_wet(k,i,j)
      if (re(k,i,j) .ge. 0.01_f.and.re(k,i,j) .le. 300.0_f ) then
      
        ! Compute Davies number (Best number)
	    NDA=(32.0_f/3.0_f)*(GRAV*rhoa_cgs*rhop_wet(k,i,j)**3._f)* r_wet(k,i,j)**3/(rmu(k)**2._f)
        x = log(NDA)            !IS LN THE NATURAL LOG?  WHAT IS RPRAT(1,J)/. WHY IS IT PRESENT IN STOKES FALL
        y= -3.18657_f+x*0.992696+x**2*(-1.53193E-3)+x**3*(-9.87059E-4)+x**4*(-5.78878E-4)+x**5*(8.55176E-5) &                 
             +x**6*(-3.27815E-6)                    
        re(k,i,j) = exp(y) 
        vf(k,i,j) = re(k,i,j) * rmu(k) / (2._f * r_wet(k,i,j) * rprat(i,j) * rhoa_cgs)
	      
	  else if (re(k,i,j) .gt. 300.0_f ) then
      
        ! compute the Bond number
        surfacetension = 76.10_f - 0.155_f * (t(k) - 273.16_f)
        NB = (16._f/3._f)*rhop_wet(k,i,j) * r_wet(k,i,j)**2 * GRAV/surfacetension
        ! compute the particle property number
        Npp=rhoa_cgs**2*surfacetension**3/(GRAV*rhop_wet(k,i,j)*rmu(k)**4)
	    x=log(NB*Npp**(1._f/6._f))
        y= -5.00015+x*5.23778+x**2*(-2.04914)+x**3*(0.475294)+x**4*(-5.42819E-2) +x**5*(2.38449E-3)                         
		re(k,i,j) = Npp**(1._f/6._f)*exp(y)	

        vf(k,i,j) =  re(k,i,j) * rmu(k) / (2._f * r_wet(k,i,j) * rprat(i,j) * rhoa_cgs)
  
      endif
      !write(*,*) "t",t(k),"p",p(k),"k",k,"i",i,"j",j,"r_wet mm",r_wet(k,i,j)*10.,"vf m/s",vf(k,i,j)/100.
    enddo      ! <i=1,NBIN>
  enddo      ! <k=1,NZ>

  ! Return to caller with particle fall velocities evaluated.
  return
end
