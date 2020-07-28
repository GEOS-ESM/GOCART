#include "MAPL_Generic.h"

module MissingFields_mod
    use, intrinsic :: iso_fortran_env, only: REAL32, REAL64

    use ESMF
    use MAPL

    use GEOS_UtilsMod

    implicit none
    private

    public get_ZLE

contains
    subroutine get_AIRDENS(PLE, T, AIRDENS)
        real(kind=REAL32), intent(in   ) :: PLE(:,:,:)
        real(kind=REAL32), intent(in   ) :: T(:,:,:)
        real(kind=REAL32), intent(  out) :: AIRDENS(:,:,:)

        AIRDENS = PLE/(MAPL_RDRY*T)
    end subroutine get_AIRDENS

    subroutine get_AIRDENS_DRYP()
        
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

    subroutine get_RH2(T, PLE, Q, RH2)
        real(kind=REAL32), intent(in   ) :: T(:,:,:)
        real(kind=REAL32), intent(in   ) :: PLE(:,:,:)
        real(kind=REAL32), intent(in   ) :: Q(:,:,:)
        real(kind=REAL32), intent(  out) :: RH2(:,:,:)

        integer :: i, j, k, n(3)

        n = size(T)

        do i=1, n(1)
            do j=1, n(2)
                do k=1, n(3)
                    RH2(i,j,k) = Q(i,j,k) / GEOS_Qsat(T(i,j,k), PLE(i,j,k))
                end do
            end do
        end do

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