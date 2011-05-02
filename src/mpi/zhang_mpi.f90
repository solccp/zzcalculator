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
    use tree
    use database_m

    implicit none
    integer(kint) :: i, level
    type(structure), target :: pah


    integer(kint) :: argc
    character(len=80) :: input_fname
    character :: okey

    integer :: error, free_node, final_size
    logical :: reach_limit
    integer :: j
    integer, parameter :: max_tree_size = 3000
!    integer, parameter :: max_tree_size = 2187
    type(tree_node_ptr), dimension(max_tree_size) :: pah_array
    type(tree_node) :: root

    integer :: local_pah_count, local_ext_count, local_rank, local_calc_count
    integer, dimension(max_tree_size) :: local_index

    integer, dimension(1000) :: vec

    call MPI_Init ( error )
    call mpi_global_init()

    level = 0


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
            read(optarg, * ) options%print_order
            if (options%print_order < 0) then 
                options%print_order = 0
            end if
        end if
        if(okey == 'b') then
            read(optarg, '(a)') options%bondlistfile
            if (image_id == 0) then
                print *, 'Load bondlist file:', trim(options%bondlistfile)
                call flush(6)
            end if
            options%has_bondlistfile = .true.
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

    root%pah => pah

    call set_max_size(maxval(pah%indexmapping))

    call build_tree(root, image_count, max_tree_size, pah_array, final_size, reach_limit)

    if ( image_id == 0 ) then
        if (reach_limit) then
            write(0, '(a)') 'warning: number of running cores are too many to be effient.'
        end if
        print *, final_size, 'jobs were produced.'
        do i=1, final_size
            vec = 0
            do j = 1, pah_array(i)%node%pah%nat
                vec(pah_array(i)%node%pah%indexmapping(j)) = 1
            end do
            write (*, '(i5, 3x, 99999(i1))' ) pah_array(i)%node%pah%nat, vec(:max_size)
        end do
        print *, 'finish printing'
    end if

    local_pah_count = final_size / image_count
    local_ext_count = mod(final_size, image_count)
    local_rank = image_count - image_id
    local_calc_count = local_pah_count
    if ( local_rank <= local_ext_count ) then   
        local_calc_count = local_calc_count + 1
    end if


    local_index = 0
    j = local_rank
    do i= 1, local_calc_count
        local_index(i) = j
        j = j + image_count
    end do

!    print *, image_id, local_index(:local_calc_count)

    do i = 1, local_calc_count
       write(*, '(i0, a, i0, a, i0)'), image_id, ' running ', i, ' of ' , local_calc_count
        call find_ZZ_polynomial(pah_array(local_index(i))%node%pah, level)
    end do

    print *, image_id, 'local jobs finished'
    
    if ( image_id == 0 ) then
        do j = 1, image_count-1
            call recv_polynomial(pah_array, j)
        end do
        call sum_up(root)
        ! ###########################
        ! # print the ZZ polynomial #
        ! ###########################
        call print_ZZ_polynomial(pah)
        call close_file()

    else
        call send_polynomial(pah_array, local_index, local_calc_count)
        call clear_tree(root)
    end if


    ! #############################################################
    ! # find recursively the ZZ polynomial of the given structure #
    ! #############################################################


    call finalize_temp_space()

    call mpi_global_finalize()
    call MPI_Finalize ( error )

end
!####################################################################################
!###################### end of program zhang_polynomial #############################
