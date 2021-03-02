! ut_GOCART - Simple ESMF/MAPL example demonstrating how to call GOCART
!
! It assumes 2 processors, so typically you will run it as
!
! % mprirun -np 2 ut_GOCART.x
!
! Arlindo da Silva <arlindo.dasilva@nasa.gov>, December 2009
!----------------------------------------------------------------------------

#  include "MAPL_Generic.h"

   Program ut_GOCART

   use ESMF
   use MAPL
   use GOCART_GridCompMod, only: SetServices

   implicit NONE

!  Basic ESMF objects being used in this example
!  ---------------------------------------------
   type(ESMF_Grid)         :: grid     ! Grid
   type(ESMF_VM)           :: vm       ! ESMF Virtual Machine
   type(ESMF_Time)         :: Time     ! Time objects
   type(ESMF_TimeInterval) :: TimeStep ! used to define a clock

!  Grid Component Objects
!  ----------------------
   type(ESMF_GridComp) :: GC 
   type(ESMF_State)    :: IMPORT
   type(ESMF_State)    :: EXPORT
   type(ESMF_Clock)    :: CLOCK

!  Basic information about the parallel environment
!         PET = Persistent Execution Threads
!  In the current implementation, a PET is equivalent 
!  to an MPI process
!  ------------------------------------------------
   integer :: myPET   ! The local PET number
   integer :: nPET    ! The total number of PETs you are running on

   integer :: status, rc
   integer :: i, j, n, im, jm

   integer :: Nx = 2, Ny=1                            ! Layout
   integer :: IM_World=72, JM_World=46, LM_WORLD=72   ! Grid dimensions

!  Coordinate variables
!  --------------------
   real(kind=8), pointer, dimension(:,:) :: lons, lats

   character(len=ESMF_MAXSTR)    :: name
   real, pointer, dimension(:,:) :: Array, newArray 

   character(len=*), parameter :: Iam = 'ut_GOCART'

!                             -----
    
    call Main()

CONTAINS

    subroutine Main()

!   Initialize the ESMF. For performance reasons, it is important
!    to turn OFF ESMF's automatic logging feature
!   -------------------------------------------------------------
    call ESMF_Initialize (LogKindFlag=ESMF_LOGKIND_NONE, vm=vm, __RC__)

!   Check the number of processors
!   ------------------------------
    call ESMF_VMGet(vm, localPET=myPET, PETcount=nPET)  
    if ( nPET /= 2 ) then
       if ( MAPL_am_I_root() ) then
          print *, 'Error: expecting 2 PETs but found ', nPET, 'PETs'
          print *, 'Try:   mpirun -np 2 ut_GOCART.x'
       end if
       _ASSERT(.FALSE.,'needs informative message')
    end if

    if ( MAPL_am_I_root() ) then
         print *
         print *, 'Starting ' // Iam // ' with ', nPET, ' PETs ...'
         print *
    end if

!   Create a global 2D Lat-Lon grid on a 2x1 layout
!   ------------------------------------------------
    Grid = MAPL_LatLonGridCreate (Name='myGrid',         &
                                  Nx = Nx, Ny = Ny,      &
                                  IM_World = IM_World,   &
                                  JM_World = JM_World,   &
                                  LM_World = LM_World,   &
                                  __RC__ )

!   Validate grid
!   -------------
    call ESMF_GridValidate(Grid,__RC__)

!   Create a clock starting at 1/1/2001 0Z with a 30 min time step
!   --------------------------------------------------------------
    call ESMF_CalendarSetDefault ( ESMF_CALKIND_GREGORIAN )
    call ESMF_TimeSet(Time, yy=2007, mm=7, dd=1,  h=0,  m=0, s=0)
    call ESMF_TimeIntervalSet( TimeStep, h=0, m=30, s=0, __RC__ )
    CLOCK = ESMF_ClockCreate ( "Clock", timeStep=TimeStep, startTime=Time, __RC__ )


!   Create states and the component
!   -------------------------------
    IMPORT = ESMF_StateCreate ( name='impGOCART', __RC__ )
    EXPORT = ESMF_StateCreate ( name='expGOCART', __RC__ )
        GC = ESMF_GridCompCreate ( name='GOCART',             &
                                   Grid=Grid,                 &
!                                   GridCompType = ESMF_ATM,   &
                                   ConfigFile='MAPL.rc',      &
                                   __RC__  )


!   Set component services
!   ----------------------
    call ESMF_GridCompSetServices ( GC, SetServices, __RC__ )

!   Initialize component
!   --------------------
    call ESMF_GridCompInitialize ( GC, importState=IMPORT, exportState=EXPORT, clock=CLOCK, __RC__ )

!   Fill in IMPORT state with reasonable values
!   -------------------------------------------
    call Fill_Import_State_ (IMPORT,__RC__)

!   Since we are not reading restarts, set the internal state with
!   reasonable profiles so that we can exercise the code
!   ---------------------------------------------------------------
    call Fill_Internal_State_ (GC,__RC__)

!   Look at states
!   --------------
    if ( MAPL_AM_I_ROOT() ) then
       call ESMF_StatePrint(IMPORT)
       call ESMF_StatePrint(EXPORT)
    end if

!   Run component
!   -------------
    call ESMF_GridCompRun ( GC, importState=IMPORT, exportState=EXPORT, clock=CLOCK, __RC__ )

!   Finalize component
!   ------------------
!!!    call ESMF_GridCompFinalize ( GC, IMPORT, EXPORT, CLOCK, __RC__ )

!   All done
!   --------
    call ESMF_Finalize(__RC__)

  end subroutine Main

!............................................................................................
  subroutine Fill_Import_State_ (IMPORT,rc)
  type(ESMF_State), intent(inout) :: IMPORT
  integer, optional, intent(out)     :: rc

!                              ----

  integer :: i1, i2, j1, j2, k1, k2

  real, pointer, dimension(:,:,:) :: ple, zle, airdens, fcld, dqdt, t, u, v, o3, rh2

  real, pointer, dimension(:,:)   :: tropp, lwi, zpbl, frlake, fraci, wet1, lai, grn, cn_prcp, ncn_prcp, &
                                     ps, sh, tsoil1, u10m, v10m, ustar, z0h

!     Get Pointers to IMPORT state
!     ----------------------------
      call MAPL_GetPointer ( IMPORT, PLE, 'PLE', __RC__ )
      call MAPL_GetPointer ( IMPORT, ZLE, 'ZLE', __RC__ )
      call MAPL_GetPointer ( IMPORT, AIRDENS, 'AIRDENS', __RC__ )
      call MAPL_GetPointer ( IMPORT, FCLD, 'FCLD', __RC__ )
      call MAPL_GetPointer ( IMPORT, DQDT, 'DQDT', __RC__ )
      call MAPL_GetPointer ( IMPORT, T, 'T', __RC__ )
      call MAPL_GetPointer ( IMPORT, U, 'U', __RC__ )
      call MAPL_GetPointer ( IMPORT, V, 'V', __RC__ )
!      call MAPL_GetPointer ( IMPORT, O3, 'O3', __RC__ )
      call MAPL_GetPointer ( IMPORT, RH2, 'RH2', __RC__ )
      call MAPL_GetPointer ( IMPORT, TROPP, 'TROPP', __RC__ )
      call MAPL_GetPointer ( IMPORT, LWI, 'LWI', __RC__ )
      call MAPL_GetPointer ( IMPORT, ZPBL, 'ZPBL', __RC__ )
      call MAPL_GetPointer ( IMPORT, FRLAKE, 'FRLAKE', __RC__ )
      call MAPL_GetPointer ( IMPORT, FRACI, 'FRACI', __RC__ )
      call MAPL_GetPointer ( IMPORT, WET1, 'WET1', __RC__ )
      call MAPL_GetPointer ( IMPORT, LAI, 'LAI', __RC__ )
      call MAPL_GetPointer ( IMPORT, GRN, 'GRN', __RC__ )
      call MAPL_GetPointer ( IMPORT, CN_PRCP, 'CN_PRCP', __RC__ )
      call MAPL_GetPointer ( IMPORT, NCN_PRCP, 'NCN_PRCP', __RC__ )
      call MAPL_GetPointer ( IMPORT, PS, 'PS', __RC__ )
      call MAPL_GetPointer ( IMPORT, SH, 'SH', __RC__ )
      call MAPL_GetPointer ( IMPORT, TSOIL1, 'TSOIL1', __RC__ )
      call MAPL_GetPointer ( IMPORT, U10M, 'U10M', __RC__ )
      call MAPL_GetPointer ( IMPORT, V10M, 'V10M', __RC__ )
      call MAPL_GetPointer ( IMPORT, USTAR, 'USTAR', __RC__ )
      call MAPL_GetPointer ( IMPORT, Z0H, 'Z0H', __RC__ )

     i1 = lbound(u,1)
     j1 = lbound(u,2)
     k1 = lbound(u,3)
     i2 = ubound(u,1)
     j2 = ubound(u,2)
     k2 = ubound(u,3)

     _ASSERT( (k2-k1+1) == 72,'needs informative message')


!     Fill typical values
!     -------------------
      do j = j1, j2
         do i = i1, i2

!           3D
!           --
            PLE(i,j,:) = (/ 1, 2, 3, 4, 6, 8, 11, 15, 21, 27, 36, 47, 61, 79, 101, 130,       &
                           165, 208, 262, 327, 407, 504, 621, 761, 929, 1127, 1364, 1645,  &
                           1979, 2373, 2836, 3381, 4017, 4764, 5638, 6660, 7851, 9236,     &
                           10866, 12783, 15039, 17693, 20792, 24398, 28606, 33388, 37003,  &
                           40612, 44214, 47816, 51405, 54997, 58584, 62170, 65769, 68147,  &
                           70540, 72931, 75313, 77711, 79623, 81046, 82485, 83906, 85344,  &
                           86765, 88201, 89636, 91071, 92516, 93921, 95376 /)       

            ZLE(i,j,:) = (/ 78676, 74222, 71032, 68578, 66390, 64345, 62371, 60419, 58455, &
                            56469, 54463, 52449, 50446, 48476, 46563, 44718, 42946, 41256, &
                            39651, 38123, 36656, 35234, 33847, 32499, 31199, 29940, 28704, &
                            27494, 26310, 25151, 24017, 22905, 21815, 20745, 19691, 18656, &
                            17629, 16609, 15589, 14559, 13514, 12470, 11475, 10487, 9469, &
                            8438, 7731, 7076, 6463, 5889, 5348, 4838, 4355, 3898, 3464, &
                            3187, 2918, 2656, 2403, 2155, 1963, 1821, 1682, 1546, 1412, &
                            1280, 1149, 1022, 896, 773, 654, 535, 417 /)

            AIRDENS(i,j,:) = (/ 2.27987766266e-05, 4.03523445129e-05, 6.19888305664e-05, 8.63075256348e-05, &
                                0.000117659568787, 0.000159025192261, 0.000209808349609, 0.000270366668701, &
                                0.000345230102539, 0.000439167022705, 0.00055980682373, 0.000717163085938, &
                                0.000923156738281, 0.00120162963867, 0.00156402587891, 0.00202178955078, &
                                0.00262451171875, 0.00339889526367, 0.00437164306641, 0.00555419921875, &
                                0.00694274902344, 0.00857543945312, 0.0105895996094, 0.0131225585938, &
                                0.0160827636719, 0.0195617675781, 0.0237731933594, 0.0287780761719, &
                                0.0347290039062, 0.0416870117188, 0.0499267578125, 0.0596313476562, &
                                0.0711669921875, 0.084716796875, 0.100830078125, 0.11865234375, 0.138671875, &
                                0.1630859375, 0.190185546875, 0.22021484375, 0.25927734375, 0.318359375, &
                                0.3720703125, 0.42138671875, 0.47265625, 0.521484375, 0.5615234375, &
                                0.6005859375, 0.638671875, 0.677734375, 0.71875, 0.759765625, 0.8017578125, &
                                0.8447265625, 0.8798828125, 0.90625, 0.9326171875, 0.958984375, 0.986328125, &
                                1.013671875, 1.03515625, 1.052734375, 1.072265625, 1.08984375, 1.10546875, &
                                1.123046875, 1.140625, 1.162109375, 1.1953125, 1.21875, 1.234375, 1.25 /)

            FCLD(i,j,:) = 1e-2 * &
                          (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                             0, 0, 0, 0, 0, 0, 16, 21, 26, 28, 0, 0, 0, 0, 0, 0, 0, &
                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0 /)

            DQDT(i,j,:) = 1e-12 * &
                          (/ 9, 11, -3, -3, -2, -18, -10, 2, 0, -3, -6, -5, -3, -1, 1, &
                             1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, &
                             1, 1, 0, 0, -5, -22, -33, 95, 474, 348, 177, 3377, 11045, &
                             11788, -5267, -7756, -17491, -19790, -10884, -6082, 8120, 4381, &
                             -10346, 8033, 69151, 77650, 61351, 46508, 33936, 23022, 15658, &
                             11598, 6469, 4861, -846, -7974, -30500, -20663, -14930 /)

            T(i,j,:) = (/ 219, 221, 223, 228, 230, 230, 232, 238, 245, 253, 259, 263, &
                          264, 262, 258, 253, 247, 239, 233, 229, 227, 227, 226, 223, &
                          222, 221, 220, 219, 218, 217, 216, 215, 214, 213, 212, 212, &
                          214, 214, 216, 219, 219, 210, 210, 218, 227, 234, 240, 245, &
                          250, 254, 257, 260, 262, 263, 265, 266, 267, 268, 269, 270, &
                          270, 270, 270, 270, 271, 271, 271, 270, 267, 265, 266, 266 /)


            U(i,j,:) = (/ -18, -13, 0, 10, 26, 36, 39, 40, 38, 37, 36, 35, 32, 28,    &
                           23, 16, 6, -2, -9, -13, -15, -16, -14, -14, -12, -12, -11, &
                          -10, -9, -5, -3, -2, 0, 1, 3, 5, 9, 13, 17, 22, 24, 26,     &
                           25, 26, 26, 22, 19, 17, 14, 12, 12, 11, 11, 11, 11, 10, 9, &
                            8, 6, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -6, -6 /)

            V(i,j,:) = (/ 20, 13, 9, 4, -1, -9, -20, -24, -25, -27, -28, -28, -26,    &
                         -25, -27, -28, -28, -28, -27, -27, -25, -23, -19, -15, -11,  &
                         -10, -9, -8, -7, -7, -8, -9, -10, -12, -14, -15, -16, -18,   &
                         -21, -22, -22, -25, -29, -25, -23, -23, -22, -20, -17, -13,  &
                         -9, -6, -4, -4, -4, -3, -2, -1, 0, 0, 0, 1, 1, 1, 2, 2,      &
                          3, 3, 3, 4, 4, 3 /)

            O3(i,j,:) = 1.E-9 * & 
                        (/ 16182, 9700, 7294, 5781, 4164, 3017, 2440, 2287, 2324, 2514,  &
                           2838, 3304, 4030, 4924, 5915, 7033, 8434, 9894, 11101, 11414, &
                           10475, 9745, 10058, 9119, 8538, 9238, 9164, 10028, 10132, 10237, &
                           9447, 7972, 7174, 5222, 4008, 3296, 2231, 1320, 768, 628, 685, &
                           676, 202, 122, 96, 88, 86, 83, 83, 84, 84, 83, 82, 81, 79, &
                           79, 77, 76, 77, 80, 84, 87, 89, 90, 89, 88, 83, 76, 69, 65, &
                           64, 64 /)
 
            RH2(i,j,:) = 1e-6 * &
                         (/ 1, 2, 2, 2, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 6, 18, 51,          &
                            129, 267, 394, 502, 682, 1135, 1603, 2076, 2820, 3792, 5120,     &
                            6806, 8912, 11597, 15397, 20386, 28168, 29755, 28748, 33875,     &
                            34058, 28657, 43458, 401856, 947266, 932618, 902344, 657227,     &
                            371583, 203370, 235108, 317872, 413086, 511719, 691407, 686524,  &
                            601563, 456055, 475098, 626954, 590821, 483399, 380860, 297852,  &
                            230958, 183594, 144288, 111084, 96558, 136963, 369629, 770508,   &
                            793946, 799805 /)
!           2D
!           --
            TROPP(i,j) = 20363.5
            LWI(i,j) = 1.
            ZPBL(i,j) = 59.
            FRLAKE(i,j) = 0.
            FRACI(i,j) = 0.
            WET1(i,j) = 0.0
            LAI(i,j) = 0.280273
            GRN(i,j) = 0.5
            CN_PRCP(i,j) = 0.0
            NCN_PRCP(i,j) = 3.18323e-10
            PS(i,j) = 96825.3 
            SH(i,j) = -28.548
            TSOIL1(i,j) = 260.014
            U10M(i,j) = -3.5
            V10M(i,j) = 2.8
            USTAR(i,j) = 0.29
            Z0H(i,j) = 0.02005

         end do
      end do

    end subroutine Fill_Import_State_

!...............................................................................................


  subroutine Fill_Internal_State_ (GC,rc)
  type(ESMF_GridComp), intent(inout) :: GC
  integer, optional, intent(out)     :: rc

!                              ----

  type(MAPL_MetaComp), pointer :: MAPL
  type(ESMF_State) :: INTERNAL

  real, pointer :: tracer(:,:,:)
  integer :: i1, i2, j1, j2, k1, k2

! Dust only for now...
  real, pointer, dimension(:,:,:) :: du001, du002, du003, du004, du005


!     Get my internal stateb out of GC
!     --------------------------------
      call MAPL_GetObjectFromGC ( GC, MAPL, __RC__ )
      call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

!     Get Pointers to IMPORT state
!     ----------------------------
      call MAPL_GetPointer ( INTERNAL, du001, 'GOCART::du001', __RC__ )
      call MAPL_GetPointer ( INTERNAL, du002, 'GOCART::du002', __RC__ )
      call MAPL_GetPointer ( INTERNAL, du003, 'GOCART::du003', __RC__ )
      call MAPL_GetPointer ( INTERNAL, du004, 'GOCART::du004', __RC__ )
      call MAPL_GetPointer ( INTERNAL, du005, 'GOCART::du005', __RC__ )

!     Local bounds
!     ------------
      _ASSERT( associated(du001), 'needs informative message' )  
      tracer => du001
      i1 = lbound(tracer,1)
      j1 = lbound(tracer,2)
      k1 = lbound(tracer,3)
      i2 = ubound(tracer,1)
      j2 = ubound(tracer,2)
      k2 = ubound(tracer,3)

     _ASSERT( (k2-k1+1) == 72,'needs informative message')

!     Fill typical values
!     -------------------
      do j = j1, j2
         do i = i1, i2

!           Dust (15Jul2008, 18W, 22N)
!           --------------------------
            if (associated(du001)) du001(i,j,:) = 1.0e-12 * &
               (/ 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 10, 15, 39, &
                  280, 1417, 1874, 1750, 1892, 2700, 4643, 5450, 7116, 9183, &
                  8906, 9838, 7757, 4577, 5341, 17317, 32480, 63214, 82073, 94646, &
                  102213, 102097, 100700, 93948, 91736, 91387, 90106, 88010, 85566, &
                  81840, 80211, 77067, 74390, 69617, 53261, 23545, 9299 /)
            if (associated(du002)) du002(i,j,:) = 1.0e-12 * &
               (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 44, &
                  568, 2231, 3024, 4148, 6563, 12137, 14465, 19064, 25990, 24448, &
                  28202, 22731, 13316, 15222, 47207, 88942, 172528, 225846, 262168, &
                  285451, 287779, 282657, 269153, 262168, 259839, 255649, 249595, &
                  243076, 232831, 228640, 220025, 212342, 199304, 153436, 68569, &
                  26863 /)
            if (associated(du003)) du003(i,j,:) = 1.0e-12 * &
               (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, &
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, &
                  30, 233, 686, 1638, 3147, 7349, 9707, 13563, 20286, 19151, &
                  23196, 19471, 11366, 12471, 37021, 71712, 141329, 189990, 224915, &
                  247732, 258443, 256114, 252389, 247732, 246801, 243076, 237255, &
                  230503, 221190, 217231, 210247, 205357, 197208, 159257, 76253, &
                  29861 /)
            if (associated(du004)) du004(i,j,:) = 1.0e-12 * &
               (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, &
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                  1, 2, 7, 33, 122, 736, 1361, 2263, 4999, 5574, 8091, 7902, &
                  5625, 5407, 12195, 23633, 47323, 69151, 87079, 105938, 120141, &
                  131084, 139000, 146451, 148081, 147149, 144588, 141096, 137836, &
                  137138, 138535, 142027, 146917, 138535, 86264, 44355 /)
            if (associated(du005)) du005(i,j,:) = 1.0e-12 * &
               (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, &
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                  1, 1, 1, 1, 1, 15, 25, 27, 116, 247, 456, 754, 852, 803, &
                  906, 1279, 2154, 3027, 3708, 5006, 6462, 9241, 12486, 15833, &
                  18860, 20635, 21887, 22614, 24302, 25437, 31433, 42085, 55531, &
                  67289, 67172, 54890 /)
            
!           To do: add other tracers...
!           ---------------------------

         end do
      end do

! --- D E B U G --- D E B U G --- D E B U G --- D E B U G --- D E B U G --- D E B U G --- 
      if (MAPL_AM_I_ROOT()) then
              print *, '----- Inside Fill_Internal_State -----'
	      print *, '     state  b o u n d s            = ',  &
              lbound(du001,1), ubound(du001,1), &
              lbound(du001,2), ubound(du001,2), &
              lbound(du001,3), ubound(du001,3)
              print *, '----- Inside Fill_Internal_State -----'
      end if
! --- D E B U G --- D E B U G --- D E B U G --- D E B U G --- D E B U G --- D E B U G --- 

    end subroutine Fill_Internal_State_

end Program ut_GOCART

