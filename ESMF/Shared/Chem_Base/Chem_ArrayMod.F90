Module Chem_ArrayMod

   private
   public Chem_Array

!  The Chem_Array, a light weight ESMF-like array
!  ----------------------------------------------
   type Chem_Array
        integer :: rank
        logical :: did_alloc = .false. ! useful to keep track of allocations
        real, pointer :: data2d(:,:)   => null()
        real, pointer :: data3d(:,:,:) => null()

!       A per-tracer scavenging efficiency in convective updrafts [km-1]
        real :: fscav = 0.0

!       A per-tracer large-scale wet removal efficiency [fraction]
        real :: fwet = 0.0

!       A kluge for doing RH affected fall velocities in CARMA
        integer :: irhFlag = 0

   end type Chem_Array

 end Module Chem_ArrayMod
