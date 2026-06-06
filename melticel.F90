! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle loss rates due to nucleation <rnuclg>:
!! Ice crystal melting only.
!! 
!! The loss rates for all particle elements in a particle group are equal.
!!
!! @author Eric Jensen, Chuck Bardeen
!! @version Jan-2000, Nov-2009
subroutine melticel(carma, cstate, iz, rc)

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
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  !  Local declarations
  integer                              :: igroup  !! group index
  integer                              :: igas    !! gas index
  integer                              :: iepart  !! element for condensing group index
  integer                              :: ienucto !! index of target nucleation element
  integer                              :: ignucto !! index of target nucleation group
  integer                              :: inuc    !! nucleating element index
  real(kind=f)                         :: bin_ice !! ice mass in the bin
  real(kind=f)                         :: max_ice !! amount of ice that can melt before raising T above T0
  real(kind=f)                         :: total_ice(NGAS)      !! total ice
  real(kind=f)                         :: total_liquid(NGAS)   !! total liquid
  real(kind=f)                         :: fraction !! fraction of liquid that can be removed

  ! Is it warm enough that anything could melt?
  if (t(iz) .gt. T0) then
  
    ! How much ice is present?
    !
    ! NOTE: This assumes all ice can be melted. You could make a different
    ! version of this routine that only includes ice that could be melted.
    call totalcondensate(carma, cstate, iz, total_ice, total_liquid, rc)

    ! Loop over particle groups.
    do igroup = 1, NGROUP

      iepart = ienconc(igroup)            ! particle number density element

      do inuc = 1, nnuc2elem(iepart)

        ienucto = inuc2elem(inuc,iepart)
      
        if (ienucto .ne. 0) then
          ignucto = igelem(ienucto)

          ! Only compute nucleation rate for ice crystal melting
          if (inucproc(iepart,ienucto) .eq. I_ICEMELT) then
    
            igas = igrowgas(iepart)      ! condensing gas

            ! Is there any ice to melt?
            if (total_ice(igas) .gt. 0._f) then
            
              ! How much ice could melt and still keep the temperature above
              ! freezing?
              max_ice = (t(iz) - T0) / rlhm(iz,igas) * CP * rhoa(iz)
              
              ! Bypass calculation if few particles are present 
              if (pconmax(iz,igroup) .gt. FEW_PC) then
              
                fraction = max_ice / total_ice(igas)
              
                ! If we melt all the ice, will it drive the temperature below T0? Then
                ! only melt a fraction of the ice.
                !
                ! NOTE: This assumes an implicit solver for particles and only being 
                ! able to melt a fraction of the mass.
                !  lg = ((1 / fraction) - 1) / dtime
                !
                ! NOTE: Choosing 99 reduces the values to 1% of the initial value.
                if (fraction < 1._f) then
                  rnuclg(:,igroup,ignucto) = (fraction / (1._f - fraction)) / dtime_orig
                else
                  rnuclg(:,igroup,ignucto) = 100._f / dtime_orig
                end if
              end if   ! pconmax(ixyz,igroup) .gt. FEW_PC
            end if
          end if       ! inucproc(iepart,ienucto) .eq. I_DROPFREEZE
        end if
      end do       ! inuc = 1,nnuc2elem(iepart)
    end do         ! igroup = 1,NGROUP
  end if

  ! Return to caller with particle loss rates due to nucleation evaluated.
  return
end
