program gocart_restart_editor

    use Chem_RegistryMod

    implicit none

    integer, parameter  :: MAX_STR_LENGTH = 1024
    integer, parameter  :: STR_LENGTH     = 256

    integer, parameter  :: SUCCESS    = 0
    integer, parameter  :: ERROR_NONE = SUCCESS

    integer, parameter  :: ERROR_SOURCE_CHEM_REGISTRY              =  1
    integer, parameter  :: ERROR_TARGET_CHEM_REGISTRY              =  2
    integer, parameter  :: ERROR_MEMORRY_ALLOCATION                =  3
    integer, parameter  :: ERROR_DIM_SIZES                         =  4
    integer, parameter  :: ERROR_GRID                              =  5
    integer, parameter  :: ERROR_UNRECOGNIZED_ARGUMENTS            =  6
    integer, parameter  :: ERROR_UNSPECIFIED_SOURCE_CHEM_REGISTRY  =  7
    integer, parameter  :: ERROR_UNSPECIFIED_TARGET_CHEM_REGISTRY  =  8
    integer, parameter  :: ERROR_UNSPECIFIED_SOURCE_GOCART_RESTART =  9
    integer, parameter  :: ERROR_UNSPECIFIED_TARGET_GOCART_RESTART = 10
    integer, parameter  :: ERROR_UNSPECIFIED_RESOLUTION            = 11
    integer, parameter  :: ERROR_UNSPECIFIED_LEVELS                = 12
    integer, parameter  :: ERROR_DUPLICATED_OPTIONS                = 13
    integer, parameter  :: ERROR_MISSING_OPTION                    = 14
    integer, parameter  :: ERROR_VERTICAL_LEVELS                   = 15
    integer, parameter  :: ERROR_CONVERSION_TO_REAL                = 16

    integer, parameter  :: REQUEST_HELP                            = 100

    
    character(len=MAX_STR_LENGTH) :: restart_src
    character(len=MAX_STR_LENGTH) :: registry_src
    character(len=MAX_STR_LENGTH) :: restart_dst
    character(len=MAX_STR_LENGTH) :: registry_dst

    integer  :: im = 0
    integer  :: jm = 0
    integer  :: lm = 0

    logical  :: dry_run = .false.

    logical  :: correct_co2 = .false.
    real     :: delta_co2 = 0.0

    integer  :: rc = 0  


    ! parse command-line arguments
    call get_arguments(source_restart_file  = restart_src,  &
                       source_registry_file = registry_src, &
                       target_registry_file = registry_dst, &
                       target_restart_file  = restart_dst,  &
                       agcm_im              = im,           &
                       agcm_jm              = jm,           &
                       agcm_lm              = lm,           &
                       dry_run              = dry_run,      &
                       correct_co2          = correct_co2,  &
                       delta_co2            = delta_co2,    &
                       rc                   = rc)

    if (rc /= SUCCESS) then
       if (rc == REQUEST_HELP) then
           stop
       else
           stop 1
       end if
    end if

    ! create GOCART restart file
    call create_restart(source_restart_file  = restart_src,  &
                        source_registry_file = registry_src, &
                        dest_registry_file   = registry_dst, &
                        dest_restart_file    = restart_dst,  &
                        agcm_im              = im,           &
                        agcm_jm              = jm,           &
                        agcm_lm              = lm,           &
                        dry_run              = dry_run,      &
                        correct_co2          = correct_co2,  &
                        delta_co2            = delta_co2,    &
                        rc                   = rc)

    if (rc /= SUCCESS) then
        stop 2
    end if

contains

subroutine get_arguments(source_restart_file, source_registry_file, &
                         target_registry_file, target_restart_file, &
                         agcm_im, agcm_jm, agcm_lm, dry_run,        &
                         correct_co2, delta_co2, rc)

    implicit none

    character(len=*), intent(out) :: source_restart_file
    character(len=*), intent(out) :: source_registry_file
    character(len=*), intent(out) :: target_registry_file
    character(len=*), intent(out) :: target_restart_file

    integer,          intent(out) :: agcm_im
    integer,          intent(out) :: agcm_jm
    integer,          intent(out) :: agcm_lm
 
    logical,          intent(out) :: dry_run

    logical,          intent(out) :: correct_co2
    real,             intent(out) :: delta_co2

    integer,          intent(out) :: rc

    ! local
    character(len=STR_LENGTH)     :: resolution
    character(len=STR_LENGTH)     :: levels
    character(len=MAX_STR_LENGTH) :: arg
    character(len=MAX_STR_LENGTH) :: buffer
    character(len=MAX_STR_LENGTH) :: delta_co2_
    integer                       :: i
    integer                       :: j
    integer                       :: argc, iargc


    source_registry_file = ''
    target_registry_file = ''
    source_restart_file  = ''
    target_restart_file  = ''

    resolution           = ''
    levels               = ''
    
    dry_run              = .false.

    correct_co2          = .false.
    delta_co2_           = ''
    delta_co2            = 0.0

    argc = iargc()

    select case(argc)
        case (0)
            rc = REQUEST_HELP
            call help()
            return

        case (1)
            call getarg(1, arg)

            if ((arg == '-h') .or. (arg == '--help') .or. (arg == '-help')) then
               rc = REQUEST_HELP
               call help()
            else
               rc = ERROR_UNRECOGNIZED_ARGUMENTS
               call error_message(rc)
            end if

            return

        case default
            
            call get_arg_('-s', '--source-registry', source_registry_file, rc)
            if (rc /= SUCCESS) then
                call error_message(rc)
                return
            end if    

            call get_arg_('-t', '--target-registry', target_registry_file, rc) 
            if (rc /= SUCCESS) then
                call error_message(rc)
                return
            end if    

            call get_arg_('-i', '--source-restart',  source_restart_file,  rc)
            if (rc /= SUCCESS) then
                call error_message(rc)
                return
            end if    

            call get_arg_('-o', '--target-restart',  target_restart_file,  rc)
            if (rc /= SUCCESS) then
                call error_message(rc)
                return
            end if

            call get_arg_('-r', '--resolution',      resolution,           rc)
            if (rc /= SUCCESS) then
                call error_message(rc)
                return
            end if    

            call get_arg_('-l', '--levels',          levels,               rc)
            if (rc /= SUCCESS) then
                call error_message(rc)
                return
            end if

            ! optional
            call get_arg_('-c', '--correct-co2',     delta_co2_,           rc)
            if (rc /= SUCCESS) then
                correct_co2 = .false.
                delta_co2_  = ''
            else
                correct_co2 = .true.
            end if

            ! options
            call get_opt_('-n', '--dry-run',         dry_run,              rc)
    end select 


    call parse_resolution(resolution, agcm_im, agcm_jm, rc)
    if (rc /= SUCCESS) then
        call error_message(rc)
        return
    end if

    call parse_levels(levels, agcm_lm, rc)
    if (rc /= SUCCESS) then
        call error_message(rc)
        return
    end if

    if (correct_co2) then
        call parse_real(delta_co2_, delta_co2, rc)
        if (rc /= SUCCESS) then
            call error_message(rc)
            return
        end if
    endif


#if (DEBUG)
    print *, 'Dry run         : ', dry_run 
    print *, 'Source Registry : ', trim(source_registry_file)
    print *, 'Target Registry : ', trim(target_registry_file)
    print *, 'Source Restart  : ', trim(source_restart_file)
    print *, 'Target Restart  : ', trim(target_restart_file)
    print *, 'Resolution      : ', trim(resolution), agcm_im, agcm_jm
    print *, 'Levels          : ', trim(levels), agcm_lm
    print *, 'Correct CO2     : ', trim(delta_co2_), correct_co2, delta_co2
#endif

    rc = 0

    return

end subroutine get_arguments


subroutine get_arg_(opt_short, opt_long, arg, rc)

    implicit none

    character(len=*), intent(in)    :: opt_short
    character(len=*), intent(in)    :: opt_long
    character(len=*), intent(inout) :: arg
    integer,          intent(out)   :: rc


    ! local
    integer                         :: i, j
    character(len=MAX_STR_LENGTH)   :: buffer
    character(len=STR_LENGTH)       :: opt_long_
    integer                         :: argc, iargc

    arg = ''
    rc  = ERROR_MISSING_OPTION

    opt_long_ = trim(opt_long) // '='

    argc = iargc()

    SHORT_OPTION: do i = 1, argc
        call getarg(i, arg)

        if ((arg == opt_short) .and. (i < argc)) then
            call getarg(i+1, arg)
            rc = SUCCESS
            return
        end if
    end do SHORT_OPTION

    LONG_OPTION: do i = 1, argc
        call getarg(i, buffer)

        j = index(buffer, trim(opt_long_))
        if (j == 1) then
            arg = buffer(len_trim(opt_long_)+1:)
            rc = SUCCESS
            return
        end if
    end do LONG_OPTION

    return 
end subroutine get_arg_


subroutine get_opt_(opt_short, opt_long, arg, rc)

    implicit none

    character(len=*), intent(in)    :: opt_short
    character(len=*), intent(in)    :: opt_long
    logical,          intent(out)   :: arg
    integer,          intent(out)   :: rc


    ! local
    integer                         :: i
    character(len=MAX_STR_LENGTH)   :: buffer
    integer                         :: argc, iargc

    arg = .false.
    rc  = SUCCESS

    argc = iargc()

    do i = 1, argc
        call getarg(i, buffer)

        if ((buffer == opt_short) .or. (buffer == opt_long)) then
            arg = .true.
            return
        end if
    end do

    return 
end subroutine get_opt_


subroutine create_restart(source_restart_file, source_registry_file, &
                          dest_registry_file, dest_restart_file,     &
                          agcm_im, agcm_jm, agcm_lm, dry_run,        &
                          correct_co2, delta_co2, rc)

    implicit none

    character(len=*), intent(in) :: source_restart_file
    character(len=*), intent(in) :: source_registry_file
    character(len=*), intent(in) :: dest_registry_file
    character(len=*), intent(in) :: dest_restart_file

    integer,          intent(in) :: agcm_im
    integer,          intent(in) :: agcm_jm
    integer,          intent(in) :: agcm_lm

    logical,          intent(in) :: dry_run

    logical,          intent(in) :: correct_co2
    real,             intent(in) :: delta_co2

    integer,          intent(out) :: rc

    ! local

    real(kind=4), allocatable, dimension(:,:) :: q  ! buffer

    type(Chem_Registry) :: r_src         ! source chem registry 
    type(Chem_Registry) :: r_dst         ! destination chem registry

    integer :: status                    ! status code
    integer :: i_src, i_dst              ! species index
    integer :: l                         ! vertical level

    logical :: found_match               ! matching tracer names

    integer, parameter  :: src = 8       ! source restart file unit
    integer, parameter  :: dst = 9       ! destination restart file unit


    rc = SUCCESS

    r_src = Chem_RegistryCreate(status, rcfile=source_registry_file)
    if (status /= SUCCESS) then
        rc = ERROR_SOURCE_CHEM_REGISTRY
        return
    end if

    r_dst = Chem_RegistryCreate(status, rcfile=dest_registry_file)
    if (status /= SUCCESS) then
        rc = ERROR_TARGET_CHEM_REGISTRY
        return
    end if
         
    ! open the input and output files
    open(unit=src, file=trim(source_restart_file), form='unformatted', access='sequential', status='old', action='read')
    open(unit=dst, file=trim(dest_restart_file),   form='unformatted', access='sequential', status='new', action='write')

    ! allocate for a 2D slice
    allocate(q(im,jm), stat=status)
    if (status /= SUCCESS) then
        rc = ERROR_MEMORRY_ALLOCATION
        return
    end if

    
    DESTINATION_TRACERS: do i_dst = r_dst%i_GOCART, r_dst%j_GOCART

        found_match = .false.
        
        print *, 'adding to restart: ' // trim(r_dst%vname(i_dst))

        rewind(unit=src)

        SOURCE_TRACERS: do i_src = r_src%i_GOCART, r_src%j_GOCART

            if (r_dst%vname(i_dst) == r_src%vname(i_src)) then
                found_match = .true.
                exit SOURCE_TRACERS
            else
                SKIP_TRACER: do l = 1, lm
                    read (unit=src) q
                end do SKIP_TRACER
            end if
        end do SOURCE_TRACERS

        if (found_match) then
            print *, trim(r_dst%vname(i_dst)) // ' = ' // trim(r_src%vname(i_src))

            COPY_TRACER: do l = 1, lm
                    read  (unit=src) q
                    if (.not. dry_run) then
                        if (trim(r_src%vname(i_src)) == 'CO2' .and. correct_co2) then
                            write (unit=dst) (q + delta_co2)
                        else
                            write (unit=dst) q
                        end if   
                    end if    
            end do COPY_TRACER
        else
            print *, trim(r_dst%vname(i_dst)) // ' will be initialized to = 0.0'

            INITIALIZE_TRACER: do l = 1, lm
                    q = 0.0
                    if (.not. dry_run) then 
                        write (unit=dst) q
                    end if    
            end do INITIALIZE_TRACER
        end if

        print *
    end do DESTINATION_TRACERS

    deallocate(q)

    close(unit=src)
    close(unit=dst)

    return 

end subroutine create_restart


subroutine help()

    implicit none

    print *
    print *, 'Usage: gogo.x -s SOURCE_CHEM_REGISTRY  \'
    print *, '              -t TARGET_CHEM_REGISTRY  \'
    print *, '              -i SOURCE_GOCART_RESTART \'
    print *, '              -o TARGET_GOCART_RESTART \'
    print *, '              -r resolution            \'
    print *, '              -l levels                \'
    print *, '              [-c DELTA_CO2]           \'
    print *, '              [-n]'
    print *, ''
    print *, '       gogo.x -h | --help'
    print *
    print *, 'Mandatory arguments to long options are mandatory for short options too.'
    print *, '  -s, --source-registry=REGISTRY         Source Chem_Registry.rc'
    print *, '  -t, --target-registry=REGISTRY         Target Chem_Registry.rc'
    print *, '  -i, --source-restart=GOCART            Source gocart_internal_rst'
    print *, '  -o, --target-restart=GOCART            Target gocart_internal_rst'
    print *, '  -r, --resolution=b|...|c48|...|IM,JM   GEOS-5 nominal resolution or IM,JM'
    print *, '  -l, --levels=LM                        Vertical levels. Typically 72 or 137.'
    print *, '  -c, --correct-co2=DELTA_CO2            Modify CO2 by adding uniform value of DELTA_CO2 [mol/mol].'
    print *, '  -n, --dry-run                          Trial run that does not write a'
    print *, '                                         restart file'
    print *, '  -h, --help                             Show this screen'
    print *

    return
end subroutine help


subroutine parse_resolution(res, im, jm, rc)

    implicit none

    character(len=*), intent(in)  :: res  ! resolution string
    integer,          intent(out) :: im   ! IM
    integer,          intent(out) :: jm   ! JM
    integer,          intent(out) :: rc   ! return code

    ! local
    integer                       :: i
    character(len=1), parameter   :: delimiter = ','

    im = -1
    jm = -1

    select case(res)
        case ('a')
                    im =   72;    jm =  36
        case ('b')
                    im =  144;    jm =  91
        case ('c')
                    im =  288;    jm = 181
        case ('d')
                    im =  576;    jm = 361
        case ('D')
                    im =  540;    jm = 361
        case ('e')
                    im = 1152;    jm = 721
        case ('E')
                    im = 1080;    jm = 721

        case ('C48', 'C90', 'C180', 'C360', 'C500', 'C720', 'C1000', 'C1440', 'C2000', 'C2880', &
              'c48', 'c90', 'c180', 'c360', 'c500', 'c720', 'c1000', 'c1440', 'c2000', 'c2880')
            read(res(2:), '(i4)', iostat=rc) im
            if (rc == 0) then
                jm = 6*im
            else
                rc = ERROR_GRID
                return
            end if

        case default
            i = index(res, delimiter)

            if (i > 0) then
                read(res(1:i) , '(i4)', iostat=rc) im
                if (rc /= 0) then
                    rc = ERROR_DIM_SIZES
                    return
                end if
 
                read(res(i+1:), '(i4)', iostat=rc) jm
                if (rc /= 0) then
                    rc = ERROR_DIM_SIZES
                    return
                end if
            else
                rc = ERROR_DIM_SIZES
                return
            end if
    end select

    if ((im > 0) .and. (jm > 0)) then
        rc = SUCCESS 
    else
        rc = ERROR_GRID
    end if

    return
end subroutine parse_resolution 


subroutine parse_levels(lev, lm, rc)

    implicit none

    character(len=*), intent(in)  :: lev  ! resolution string
    integer,          intent(out) :: lm   ! JM
    integer,          intent(out) :: rc   ! return code

    read(lev, '(i3)', iostat=rc) lm

    if (rc == 0) then
        rc = SUCCESS
    else
        rc = ERROR_VERTICAL_LEVELS
    end if

end subroutine parse_levels


subroutine parse_real(str, v, rc)

    implicit none

    character(len=*), intent(in)  :: str  ! string
    real,             intent(out) :: v    ! real
    integer,          intent(out) :: rc   ! return code

    read(str, *, iostat=rc) v

    if (rc == 0) then
        rc = SUCCESS
    else
        rc = ERROR_CONVERSION_TO_REAL
    end if

end subroutine parse_real


subroutine error_message(rc)

    implicit none

    integer, intent (in) :: rc

    select case (rc)
        case (ERROR_UNRECOGNIZED_ARGUMENTS)
            print *, 'Missing source and target files, horizontal resolution.'
            print *, "Try 'gogo.x --help' for more information."

        case (ERROR_DIM_SIZES)
            print *, 'Please provide the restart file resolution as two numbers'
            print *, "separated by comma, e.g., '-r 576,361'"

        case (ERROR_GRID)
            print *, "Unrecognized grid. Please use one of the nominal GEOS-5 grids, e.g., 'c90', 'd', ..."

        case (ERROR_UNSPECIFIED_SOURCE_CHEM_REGISTRY)
            print *, 'Please provide source chem registry file.'

        case (ERROR_UNSPECIFIED_TARGET_CHEM_REGISTRY)
            print *, 'Please provide target chem registry file.'

        case (ERROR_UNSPECIFIED_SOURCE_GOCART_RESTART)
            print *, 'Please provide input GOCART restart file.'

        case (ERROR_UNSPECIFIED_TARGET_GOCART_RESTART)
            print *, 'Please provide output GOCART restart file.'

        case (ERROR_UNSPECIFIED_RESOLUTION)
            print *, 'Please provide horizontal resolution.'
            
        case (ERROR_DUPLICATED_OPTIONS)
            print *, "One or more duplicated options. Try 'gogo.x --help' for more information."

        case (ERROR_MISSING_OPTION)
            print *, "Missing mandatory option. Try 'gogo.x --help' for more information."

        case default  
         print *, "Try 'gogo.x --help' for more information."
    end select 

    print *
    print *

end subroutine

end program gocart_restart_editor

