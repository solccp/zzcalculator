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
    use mpi
    use mpi_global
    implicit none
    integer(kint) :: i,nhex,level=0
    type(structure) :: pah


    integer(kint) :: argc
    integer(kint) :: input_unit
    character(len=80) :: input_fname
    character :: okey

    integer :: error

    call MPI_Init ( error )
    call mpi_global_init()

    argc = command_argument_count()

    if (argc < 1) then
        stop 'abnormal termination: wrong arguments'
    end if

    call initialize_options()

    do
        okey = getopt('Pl:')
        if(okey == '>') exit
        if(okey == '!') then
            write(*,*) 'unknown option: ', trim(optarg)
            stop
        end if
        
        if(okey == 'P') then
            options%print_intermediate_structures = .true.
        end if
        if(okey == 'l') then
            read(optarg, '(i)') options%print_order
            if (options%print_order < 0) then 
                options%print_order = 0
            end if
        end if
        if(okey == '.') then
            input_fname = optarg
        end if
    end do

    if (options%print_intermediate_structures .and. image_count>1) then
        if (image_id == 0) then
            print*, 'running with MPI doesn''t support intermediate structure output'
        end if
        call MPI_Finalize ( error )
        stop
    end if


    ! ############################################################
    ! # read initial geometry data and create topological matrix #
    ! ############################################################
    call read_input(input_fname, pah)

    call initialize_temp_space(pah%nat)

    if (image_id == 0) then
        call find_ZZ_polynomial(pah, level)
    else
        call recv_structure(pah, level)
        call find_ZZ_polynomial(pah, level)
        call send_polynomial(pah)
    end if
    print *, image_id, 'finished jobs'

    ! #############################################################
    ! # find recursively the ZZ polynomial of the given structure #
    ! #############################################################

!    call MPI_BARRIER(MPI_COMM_WORLD, error)

    ! ###########################
    ! # print the ZZ polynomial #
    ! ###########################
    if ( image_id == 0 ) then
        call print_ZZ_polynomial(pah)
        call close_file()
    end if

    call finalize_temp_space()

    call mpi_global_finalize()
    call MPI_Finalize ( error )

end
!####################################################################################
!###################### end of program zhang_polynomial #############################
