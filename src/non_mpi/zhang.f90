!######################### program zhang_polynomial #################################
!####################################################################################
program zhang_polynomial
!
! This program calculates the Zhang-Zhang polynomial for benzenoid structures;
! the ZZ polynomial contains beside other quantities also the Clar number,
! Clar count, and Kekule number of a given structure
!   
! Reference:
!   I. Gutman, B. Furtula, and A. T. Balaban
!   Polycyclic Aromatic Compounds 26 pp.17-35, 2006
!
    use types_module
    use structure_module
    use options_m
    use getopt_m
    use output
    use temp_space

    implicit none
    integer(kint) :: i, level
    type(structure) :: pah


    integer(kint) :: argc
    character(len=80) :: input_fname
    character :: okey

    argc = command_argument_count()

    if (argc < 1) then
        stop 'abnormal termination: wrong arguments'
    end if

    call initialize_options()

    do
        okey = getopt('Pl:b:')
        if(okey == '>') exit
        if(okey == '!') then
            write(*,*) 'unknown option: ', trim(optarg)
            stop
        end if
        
        if(okey == 'P') then
            options%print_intermediate_structures = .true.
        end if
        if(okey == 'l') then
            read(optarg, *) options%print_order
            if (options%print_order < 0) then 
                options%print_order = 0
            end if
        end if
        if(okey == 'b') then
            read(optarg, *) options%bondlistfile
            options%has_bondlistfile = .true.
        end if
        if(okey == '.') then
            input_fname = optarg
        end if
    end do


    ! ############################################################
    ! # read initial geometry data and create topological matrix #
    ! ############################################################
    call read_input(input_fname, pah)

    call initialize_temp_space(pah%nat)
    
    level = 0
    call find_ZZ_polynomial(pah, level)
    
    ! ###########################
    ! # print the ZZ polynomial #
    ! ###########################
    call print_ZZ_polynomial(pah)
    call close_file()

    ! #############################################################
    ! # find recursively the ZZ polynomial of the given structure #
    ! #############################################################


    call finalize_temp_space()


end
!####################################################################################
!###################### end of program zhang_polynomial #############################
