             ! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine evaluates supersaturations <supsatl> and <supsati> for all gases.
!!
!! @author Andy Ackerman, Chuck Bardeen
!! @version Dec-1995, Aug-2010
subroutine supersat(carma, cstate, iz, igas, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(in)                  :: iz      !! z index
  integer, intent(in)                  :: igas    !! gas index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  real(kind=f)  :: rvap
  real(kind=f)  :: gc_cgs

  ! Calculate vapor pressures.
  call vaporp(carma, cstate, iz, igas, rc)

  ! Define gas constant for this gas
  rvap = RGAS / gwtmol(igas)

  gc_cgs = gc(iz,igas) / (zmet(iz)*xmet(iz)*ymet(iz))

  supsatl(iz,igas) = (gc_cgs * rvap * t(iz) - pvapl(iz,igas)) / pvapl(iz,igas)
  supsati(iz,igas) = (gc_cgs * rvap * t(iz) - pvapi(iz,igas)) / pvapi(iz,igas)

  ! For subgrid scale clouds, the supersaturation needs to be increased to represent
  ! the saturation ratio in the cloudy part of the gridbox. The amount
  ! that it needs to be increased is calculated by the parent model.
  !
  ! NOTE: Since the model now includes liquid cloud, the supsatl can become
  ! greater than 0. If it was just an ice model, then supsatl (scaled by qsatfac)
  ! should not be greater than 0; however, it is left to the parent model to control
  ! this
  ! behavior.
  !
  ! NOTE: Potentially qsatfac should be a function of the gas, but for now
  ! this will serve the purpose.
  if (do_incloud) then
    supsatl(iz,igas) = (gc_cgs * rvap * t(iz) / qsatfac(iz) - pvapl(iz,igas)) / pvapl(iz,igas)
    supsati(iz,igas) = (gc_cgs * rvap * t(iz) / qsatfac(iz) - pvapi(iz,igas)) / pvapi(iz,igas)
  end if

  return
end
