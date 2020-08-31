!
! Simple code to read a GEOS-5 Chem_Registry.rc file and create a GrADS control
! (ctl) file for the binary "gocart_internal_restart" file.
!

      Program Ctl

        use Chem_RegistryMod

        integer im, jm, km
        character(len=512) :: reg, dset

        type(Chem_Registry) :: r

        call init_()

        r = Chem_RegistryCreate ( ier, rcfile = reg )

        call print_ctl()

 Contains

!......................................................................

        subroutine print_ctl()

        if ( trim(dset) .ne. '*none*' ) then

          write(*,'(a)') 'DSET   '//trim(dset)
          write(*,'(a)') 'TITLE Chem Restart File'
          write(*,'(a)') 'OPTIONS sequential'
          write(*,'(a)') 'UNDEF 1.0E20'
          write(*,'(a,i5,a,f10.6)') 'XDEF ', im, ' LINEAR -180 ', 360./im
          write(*,'(a,i5,a,f10.6)') 'YDEF ', jm, ' LINEAR  -90 ', 180./(jm-1)
          write(*,'(a,i5,a)') 'ZDEF ', km, ' LINEAR  1 1 '
          write(*,'(a)') 'TDEF 1 LINEAR  05feb1960 6hr'
          write(*,'(a,i5)') 'VARS ', r%n_GOCART
          do i = r%i_GOCART, r%j_GOCART
             write(*,'(a,2i3,1x,a)') trim(r%vname(i)), km, 0, trim(r%vtitle(i)) 
          end do
          write(*,'(a)') 'ENDVARS '

        else

          do i = r%i_GOCART, r%j_GOCART
             call lower_(r%vname(i))
          end do
          write(*,'(a)') ' '
          write(*,'(128(a))') &
            'lats4d.sh -i old.ctl -o new -format sequential -v -vars', &
            (' '//trim(r%vname(i)), i=r%i_GOCART, r%j_GOCART)
          write(*,'(a)') ' '

        end if

        end subroutine print_ctl

!......................................................................

        subroutine init_()

!       Parses CLI

          integer iargc, argc
          character(len=512) :: res
          
          argc = iargc()
          if(argc .lt. 1) call usage_(' ')
          call GetArg(1,reg)
          if(argc .gt. 1) then
	     call GetArg(2,dset)
          else
             dset = '*none*'
          end if
          if(argc .gt. 2) then
             call GetArg(3,res)
          else
             res = 'd72'
          end if

          if ( res(2:3) .eq. '72' ) then
             km = 72
          else
             call usage_('cannot handle levels in resolution '//trim(res))
          endif

          if ( res(1:1) .eq. 'd' ) then
             im = 576
             jm = 361
          else if ( res(1:1) .eq. 'D' ) then
             im = 540
             jm = 361
          else if ( res(1:1) .eq. 'c' ) then
             im = 288
             jm = 181
          else if ( res(1:1) .eq. 'b' ) then
             im = 144
             jm = 91
          else if ( res(1:1) .eq. 'e' ) then
             im = 1152
             jm = 721
          else if ( res(1:1) .eq. 'E' ) then
             im = 1080
             jm = 721
          else
             call usage_('cannot handle horizontal resolution '//trim(res))
          end if

        end subroutine init_

!......................................................................

        subroutine lower_(str)
          character(len=*) str
          integer ich
          do i = 1, len(str)
             ich = ichar(str(i:i))
             if ((ich .GE. ichar('A')).AND.(ich .LE. ichar('Z'))) then
                   str(i:i) = char(ichar('a') + ich - ichar('A'))
             end if
          end do
        end subroutine lower_
!......................................................................

        subroutine usage_(msg)
          character(len=*) msg

         print *, "ctl_crst - GrADS CTL for Chem Restarts"
         print *
         print *, "Usage:"
         print *
         print *, "    ctl_crst.x Chem_Registry.rc [rs_filename [res] ]"
         print *
         PRINT *, "where *res* is the resolution/levels: b72, c72, ..."
         print *, "with the following resolution definitions:"
         print *, "              a:   72 x  46"
         print *, "              b:  144 x  91"
         print *, "              c:  288 x 181"
         print *, "              D:  540 x 361"
         print *, "              d:  576 x 361"
         print *, "              E: 1080 x 721"
         print *, "              e: 1152 x 721"
         print *
         print *, "Omit the 'rs_filename' name to get a sample lats4d line"
         print *, "that can be used to convert an old restart file"
         print *, "to a new one, e.g., "
         print *
         print *, "    ctl_crst.x Chem_Registry.rc"
         print *
         print *, msg
         print *

         call exit(1)

      end subroutine usage_

    end Program Ctl
