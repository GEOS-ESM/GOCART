

!
! Simple testing of Mie Tables
!

program ut_Mie

   use Chem_MieMod
   use m_die
   implicit NONE
   

   real :: rh
   real :: pmom(751,4), tau, g, ssa, idxChannel
   type(Chem_mie) :: mieTables

   integer :: iq, rc, idxTable, m, i
   integer, parameter :: nq = 15
   character(len=*), parameter :: vname(nq) = &
   (/ 'du001   ', 'du002   ', 'du003   ', 'du004   ', 'du005   ', &
      'ss001   ', 'ss002   ', 'ss003   ', 'ss004   ', 'ss005   ', &
      'BCphobic', 'BCphilic', 'OCphobic', 'OCphilic', 'SO4     ' /)

! Create 1 time
!----------------------
   print *, '1ere creation '
   print *, '------------------------------------------------- '
   mieTables = Chem_MieCreate('Aod_EOS.rc', rc )
   if ( rc /= 0 ) call die('ut_mie','oh, oh, cannot create the table')
         
   print *, 'nMom = ', mieTables%nMom
   print *, 'nPol = ', mieTables%nPol
   print *, 'nch = ', mieTables%nch
    do m = 1, mieTables%nch
        if ( abs(550.0 - (1.e9)*mieTables%channels(m)) < 1. ) then
           idxChannel = m
           print *, 'idxchannel', idxChannel
           exit
         end if
      end do

   do iq = 1, nq
 
      idxTable = Chem_MieQueryIdx(mieTables,vname(iq),rc)
      if ( rc /= 0 ) &
           call die('ut_mie','cannot do Chem_MieQueryIdx')
      print *
!      print *, '---> idx for ',vname(iq), ' is ', idxTable
      rh = 0.99
      call Chem_MieQuery(mieTables, idxTable, idxChannel, &
                         1.0, rh, ssa =ssa, gasym=g, tau=tau, pmom = pmom)

!      print*, 'tau', tau, 'g', g, 'ssa', ssa
!      print*,'pmom', pmom(11,1), pmom(20,2), pmom(20,3), pmom(20,4)

!      do i = 1, 751, 10
!if 0
!         print *, vname(iq)//' --> ', i, &
!              maxval(abs(mietables%vtableuse%pmom(:,:,:,i,1))), & 
!              maxval(abs(mietables%vtableuse%pmom(:,:,:,i,2))), & 
!              maxval(abs(mietables%vtableuse%pmom(:,:,:,i,3))), & 
!              maxval(abs(mietables%vtableuse%pmom(:,:,:,i,4)))
!else
!         print *, vname(iq)//' --> ', i, pmom(i,:)
!endif
!      end do

   end do

! destroy for 1 time
!--------------------------
     call Chem_MieDestroy(mieTables, rc)
     if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
!     return
     end if

! create 2nd time
!---------------------------  

   print *, '2eme creation '
   print *, '------------------------------------------------- '
   mieTables = Chem_MieCreate('Aod_EOS.rc', rc )
   if ( rc /= 0 ) call die('ut_mie','oh, oh, cannot create the table')
         
   print *, 'nMom = ', mieTables%nMom
   print *, 'nPol = ', mieTables%nPol
   print *, 'nch = ', mieTables%nch
    do m = 1, mieTables%nch
        if ( abs(550.0 - (1.e9)*mieTables%channels(m)) < 1. ) then
           idxChannel = m
           print *, 'idxchannel', idxChannel
           exit
         end if
      end do

   do iq = 1, nq
 
      idxTable = Chem_MieQueryIdx(mieTables,vname(iq),rc)
      if ( rc /= 0 ) &
           call die('ut_mie','cannot do Chem_MieQueryIdx')
      print *
!      print *, '---> idx for ',vname(iq), ' is ', idxTable
      rh = 0.99
      call Chem_MieQuery(mieTables, idxTable, idxChannel, &
                         1.0, rh, ssa =ssa, gasym=g, tau=tau, pmom = pmom)

!      print*, 'tau', tau, 'g', g, 'ssa', ssa
!      print*,'pmom', pmom(11,1), pmom(20,2), pmom(20,3), pmom(20,4)
   end do

! destroy for 2 time
!--------------------------
     call Chem_MieDestroy(mieTables, rc)
     if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
!     return
     end if

! create 3rd time
!--------------------------- 

   print *, '3eme creation '
   print *, '------------------------------------------------- '
   mieTables = Chem_MieCreate('Aod_EOS.rc', rc )
   if ( rc /= 0 ) call die('ut_mie','oh, oh, cannot create the table')
         
   print *, 'nMom = ', mieTables%nMom
   print *, 'nPol = ', mieTables%nPol
   print *, 'nch = ', mieTables%nch
    do m = 1, mieTables%nch
        if ( abs(550.0 - (1.e9)*mieTables%channels(m)) < 1. ) then
           idxChannel = m
           print *, 'idxchannel', idxChannel
           exit
         end if
      end do

   do iq = 1, nq
 
      idxTable = Chem_MieQueryIdx(mieTables,vname(iq),rc)
      if ( rc /= 0 ) &
           call die('ut_mie','cannot do Chem_MieQueryIdx')
      print *
!      print *, '---> idx for ',vname(iq), ' is ', idxTable
      rh = 0.99
      call Chem_MieQuery(mieTables, idxTable, idxChannel, &
                         1.0, rh, ssa =ssa, gasym=g, tau=tau, pmom = pmom)

!      print*, 'tau', tau, 'g', g, 'ssa', ssa
!      print*,'pmom', pmom(11,1), pmom(20,2), pmom(20,3), pmom(20,4)
   end do

! destroy for 3rd time
!--------------------------
     call Chem_MieDestroy(mieTables, rc)
     if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
!    return
     end if

! create 4 time
!--------------------------- 

   print *, '4eme creation '
   print *, '------------------------------------------------- '
   mieTables = Chem_MieCreate('Aod_EOS.rc', rc )
   if ( rc /= 0 ) call die('ut_mie','oh, oh, cannot create the table')
         
   print *, 'nMom = ', mieTables%nMom
   print *, 'nPol = ', mieTables%nPol
   print *, 'nch = ', mieTables%nch
    do m = 1, mieTables%nch
        if ( abs(550.0 - (1.e9)*mieTables%channels(m)) < 1. ) then
           idxChannel = m
           print *, 'idxchannel', idxChannel
           exit
         end if
      end do

   do iq = 1, nq
 
      idxTable = Chem_MieQueryIdx(mieTables,vname(iq),rc)
      if ( rc /= 0 ) &
           call die('ut_mie','cannot do Chem_MieQueryIdx')
      print *
!      print *, '---> idx for ',vname(iq), ' is ', idxTable
      rh = 0.99
      call Chem_MieQuery(mieTables, idxTable, idxChannel, &
                         1.0, rh, ssa =ssa, gasym=g, tau=tau, pmom = pmom)

!      print*, 'tau', tau, 'g', g, 'ssa', ssa
!      print*,'pmom', pmom(11,1), pmom(20,2), pmom(20,3), pmom(20,4)
   end do

! destroy for 4rd time
!--------------------------
     call Chem_MieDestroy(mieTables, rc)
     if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
!    return
     end if

! create 5 time
!--------------------------- 

   print *, '5eme creation '
   print *, '------------------------------------------------- '
   mieTables = Chem_MieCreate('Aod_EOS.rc', rc )
   if ( rc /= 0 ) call die('ut_mie','oh, oh, cannot create the table')
         
   print *, 'nMom = ', mieTables%nMom
   print *, 'nPol = ', mieTables%nPol
   print *, 'nch = ', mieTables%nch
    do m = 1, mieTables%nch
        if ( abs(550.0 - (1.e9)*mieTables%channels(m)) < 1. ) then
           idxChannel = m
           print *, 'idxchannel', idxChannel
           exit
         end if
      end do

   do iq = 1, nq
 
      idxTable = Chem_MieQueryIdx(mieTables,vname(iq),rc)
      if ( rc /= 0 ) &
           call die('ut_mie','cannot do Chem_MieQueryIdx')
      print *
!      print *, '---> idx for ',vname(iq), ' is ', idxTable
      rh = 0.99
      call Chem_MieQuery(mieTables, idxTable, idxChannel, &
                         1.0, rh, ssa =ssa, gasym=g, tau=tau, pmom = pmom)

!      print*, 'tau', tau, 'g', g, 'ssa', ssa
!      print*,'pmom', pmom(11,1), pmom(20,2), pmom(20,3), pmom(20,4)
   end do

! destroy for 5rd time
!--------------------------
     call Chem_MieDestroy(mieTables, rc)
     if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
!    return
     end if



 end program ut_Mie
