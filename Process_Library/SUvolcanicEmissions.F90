   subroutine SUvolcanicEmissions (nVolc, vStart, vEnd, vSO2, vElev, vCloud, iPoint, &
                                   jPoint, nhms, SO2EMVN, SO2EMVE, SO2, nSO2, SU_emis, km, cdt, grav,&
                                   hghte, delp, area, vLat, vLon, rc)
! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: nVolc     ! number of emissions
   integer, dimension(:), intent(in) :: vStart ! emission start time [sec]
   integer, dimension(:), intent(in) :: vEnd   ! emission end time [sec]
   real, dimension(:), intent(in)    :: vSO2   ! volcanic emission from file [kg]
   real, dimension(:), intent(in)    :: vCloud ! top elevation of emissions [m]
   integer, dimension(:), intent(in) :: iPoint, jPoint ! sub-domain locations of volcanos
   integer, intent(in) :: nhms ! current model time [sec]
   integer, intent(in) :: nSO2   ! index of SO2 relative to other sulfate tracers
   integer, intent(in) :: km   ! number of model levels
   real, intent(in)    :: cdt  ! model time step [sec]
   real, pointer, dimension(:,:,:) :: hghte     ! top of layer geopotential height [m]
   real, intent(in)    :: grav ! gravity [m sec-1]
!   real, dimension(:,:,:), intent(in) :: airdens ! layer air density [kg/m^3]
   real, dimension(:,:,:), intent(in) :: delp  ! pressure thickness [Pa]
   real, dimension(:,:), intent(in)   :: area  ! area of grid cell [m^2]
   real, dimension(:), intent(in)     :: vLat  ! latitude specified in file [degree]
   real, dimension(:), intent(in)     :: vLon  ! longitude specified in file [degree]
! !INOUT PARAMETERS:
  real, pointer, dimension(:,:), intent(inout) :: SO2EMVN ! non-explosive volcanic emissions [kg m-2 s-1]
  real, pointer, dimension(:,:), intent(inout) :: SO2EMVE ! explosive volcanic emissions [kg m-2 s-1]
  real, pointer, dimension(:,:,:), intent(inout) :: SO2 ! SO2 [kg kg-1]
  real, pointer, dimension(:,:,:), intent(inout) :: SU_emis      ! SU emissions, kg/m2/s
  real, dimension(:), intent(inout) ::  vElev ! bottom elevation of emissions [m]

! !OUTPUT PARAMETERS:
  integer, optional, intent(out)   :: rc    ! Error return code:
                                            !  0 - all is well

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! 22July2020 E.Sherman
!
! !Local Variables
   integer  ::  i, j, it
   real, dimension(:,:,:), allocatable  :: emissions_point
   real :: so2volcano

   real :: hup, hlow, dzvolc, dz, z1, k
   real :: deltaSO2v
   real, dimension(:,:), allocatable :: z0
   real, allocatable, dimension(:,:) :: srcSO2volc
   real, allocatable, dimension(:,:) :: srcSO2volce

!EOP
!-------------------------------------------------------------------------
!  Begin

   if (nVolc > 0) then

   allocate(srcSO2volc, mold=area)
   allocate(srcSO2volce, mold=area)
   srcSO2volc = 0.
   srcSO2volce = 0.

   if (associated(SU_emis)) SU_emis = 0.0
   if (associated(SO2EMVN)) SO2EMVN = 0.
   if (associated(SO2EMVE)) SO2EMVE = 0.

   allocate(z0, mold=area)
   z0 = hghte(:,:,km)

    do it = 1, nVolc
       so2volcano = 0.
       i = iPoint(it)
       j = jPoint(it)

!      Skip this volcano?
       if (i<1 .or. j<1) cycle ! volcano not in sub-domain

!      Check time against time range of eruption
       if(nhms < vStart(it) .or. nhms >= vEnd(it)) cycle

!      Emissions per volcano
       if(area(i,j) > 1.) then
          so2volcano = vSO2(it) / area(i,j)     ! to kg SO2/sec/m2
          so2volcano = max(so2volcano,tiny(so2volcano))
       endif

!        Distribute in the vertical
!        Database provides altitude of top of volcano cone (vElev) and altitude
!        of plume top (vCloud).  If vCloud != vElev then distribute emissions
!        in top 1/3 of column extending from vElev to vCloud (case of explosive
!        eruption), else put emissions in grid cell containing vElev (degassing)
!        --------------------------
         hup  = vCloud(it)
         hlow = vElev(it)
         if (hup .ne. hlow) then
            hlow = hup - (hup-hlow)/3.
         endif

!        Diagnostic - sum of volcanos
!        ----------------------------
         if (hup .eq. hlow) then
            srcSO2volc(i,j) = srcSO2volc(i,j) + so2volcano
         else
            srcSO2volce(i,j) = srcSO2volce(i,j) + so2volcano
         endif

         dzvolc = hup-hlow
         do k = km, 1, -1
            z1 = hghte(i,j,k-1) ! geopotential altitude at gridbox top
            dz = z1-z0(i,j)     ! thickness of gridbox
            deltaSO2v = 0.

!           Volcano is above this level
!           ---------------------------
            if(z1 .lt. hlow) then
               z0(i,j) = z1
               cycle
            end if

!           Volcano is below this level (except at surface)
!           -----------------------------------------------
            if(z0(i,j) .gt. hup .and. k .ne. km) then
               z0(i,j) = z1
               cycle
            end if

!           Volcano is in this level
!           ------------------------
            if( (k .eq. km .and. z0(i,j) .gt. hup) .or. &     ! below surface
                 (z0(i,j) .le. hlow .and. z1 .ge. hup) ) then ! in level
               deltaSO2v = so2volcano

!           Volcano only partly in level                       ! Cell:
!           ----------------------------
            else if (z0(i,j) .lt. hlow .and. z1 .lt. hup) then ! has bottom of cloud
               deltaSO2v = (z1-hlow)/dzvolc*so2volcano

            else if (z0(i,j) .gt. hlow .and. z1 .gt. hup) then ! has top of cloud
               deltaSO2v = (hup-z0(i,j))/dzvolc*so2volcano

            else                                               ! is filled with cloud
               deltaSO2v = dz/dzvolc*so2volcano
            end if

            z0(i,j) = z1
            so2(i,j,k) = so2(i,j,k) + deltaSO2v*cdt*grav/delp(i,j,k)

      end do ! k
   enddo     ! it
  end if ! nVolc > 0

  if (associated(SO2EMVN)) SO2EMVN = SO2EMVN + srcSO2volc
  if (associated(SO2EMVE)) SO2EMVE = SO2EMVE + srcSO2volce
  if (associated(SU_emis)) SU_emis(:,:,nSO2) = SU_emis(:,:,nSO2) + srcSO2volc + srcSO2volce

  __RETURN__(__SUCCESS__)
  end subroutine SUvolcanicEmissions
