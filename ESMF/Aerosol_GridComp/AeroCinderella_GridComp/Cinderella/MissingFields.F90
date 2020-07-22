#include "MAPL_Generic.h"

module MissingFields_mod
    use, intrinsic :: iso_fortran_env, only: REAL32, REAL64

    use ESMF
    use MAPL

    implicit none
    private

    public get_ZLE

contains
    subroutine get_AIRDENS()
        
    end subroutine get_AIRDENS
    
    subroutine get_AIRDENS_DRYP(PLE, T, AIRDENS_DRYP)
        real(kind=REAL32), intent(in   ) :: PLE(:,:,:)
        real(kind=REAL32), intent(in   ) :: T(:,:,:)
        real(kind=REAL32), intent(  out) :: AIRDENS_DRYP(:,:,:)

        AIRDENS_DRYP = PLE/(MAPL_RDRY*T)
    end subroutine get_AIRDENS_DRYP

    subroutine get_DELP(PRESS, DELP)
        real(kind=REAL32), intent(in   ) :: PRESS(:,:,:)
        real(kind=REAL32), intent(  out) :: DELP(:,:,:)

        integer :: i, n(3)

        n = size(PRESS)

        do i=1, n(3)-1
            DELP(:,:,i) =  PRESS(:,:,i+1) - PRESS(:,:,i)
        end do

    end subroutine get_DELP
    
    subroutine get_FRACI()
        
    end subroutine get_FRACI
    
    subroutine get_FRLAKE()
        
    end subroutine get_FRLAKE
    
    subroutine get_FROCEAN()
        
    end subroutine get_FROCEAN
    
    subroutine get_PS()
        
    end subroutine get_PS

    subroutine get_RH2(RH2)
        real(kind=REAL32), intent(out) :: RH2(:,:,:)

    end subroutine

    subroutine get_DZ()

    end subroutine
    
    subroutine get_ZLE(ZLE)
        real(kind=REAL32), intent(inout) :: ZLE(:,:,:)

        ZLE = ZLE/MAPL_GRAV
    end subroutine get_ZLE

    subroutine get_PFI_LSAN()

    end subroutine get_PFI_LSAN

    subroutine get_PFL_LSAN()

    end subroutine get_PFL_LSAN
end module MissingFields_mod