module aerosol_driver
    use ESMF
    use NUOPC

    implicit none
    private

    public SetServices
contains
    subroutine SetServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        rc = ESMF_SUCCESS
    end subroutine SetServices
end module aerosol_driver