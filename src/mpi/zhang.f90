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
    use database_m
    use tree
    use mpi
    use mpi_global

    implicit none
    integer(kint) :: i, level
    type(structure), target :: pah

!========================================
    integer(kint) :: argc
    character(len=200) :: input_fname
    character :: okey
!========================================
    
    type(tree_node), pointer :: root
    integer :: max_tree_size
    logical :: reach_limit
    type(database_entry), pointer :: head
!=======================================
! MPI specific 
    type(tree_node_ptr), dimension(:), allocatable :: pah_array
    integer :: local_pah_count, local_ext_count, local_rank, local_calc_count
    integer, dimension(:), allocatable :: local_index
    integer :: error, free_node, final_size
    integer :: j

    integer, dimension(:) ,allocatable :: dis
!=======================================
    call MPI_Init ( error )
    call mpi_global_init()

    max_tree_size = 5000

    argc = command_argument_count()

    if (argc < 1) then
        stop 'abnormal termination: wrong arguments'
    end if

    call initialize_options()

    do
        okey = getopt('Pfl:b:d:')
        if(okey == '>') exit
        if(okey == '!') then
            write(*,*) 'unknown option: ', trim(optarg)
            stop
        end if
        if(okey == 'd') then
            read(optarg, *) max_tree_size
        end if
        if(okey == 'P') then
            options%print_intermediate_structures = .true.
        end if
        if(okey == 'f') then
            options%force_print_structures = .true.
        end if

        if(okey == 'l') then
            read(optarg, *) options%print_order
            if (options%print_order < 0) then 
                options%print_order = 0
            end if
        end if
        if(okey == 'b') then
            read(optarg, '(a)') options%bondlistfile
            print *, 'Load bondlist file:', trim(options%bondlistfile)
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
    if ( options%print_intermediate_structures .and. pah%nat > 50 .and. .not. options%force_print_structures ) then
        print *, 'warning: # of atoms > 50, disabling intermediate structure printing'
        print *, '  use -P -f option to print the intermediate structures forcely'
        options%print_intermediate_structures = .false.
    end if


    call initialize_temp_space(pah%nat)
    
    level = 0

    if ( .not. options%print_intermediate_structures ) then
        allocate(root)
        root%pah => pah

        call set_max_size(int(maxval(pah%indexmapping)))

        call build_tree(root, max_tree_size, reach_limit)
 
        final_size = 0
        head => database_head
        do while(associated(head))
            if (.not. head%node%hasChild) then
                final_size = final_size + 1
            end if
            head => head%next
        end do
      
        allocate(pah_array(final_size)) 
        allocate(local_index(final_size)) 


        head => database_head
        i = 0
        do while(associated(head))
            if (.not. head%node%hasChild) then
                i = i + 1
                pah_array(i)%node => head%node
                
!                write(optarg, '(i0,a)') i, '.xyz'
!                print *, trim(optarg)
!                call write_xyz(head%node%pah, optarg)

!                print *, 'running (', i, '/' , total_required_strs , ')...'
!                call find_ZZ_polynomial(head%node%pah, level)
                
                write(*, '(a,i3,i5)') char(head%key), head%hits, head%node%pah%nat
!                write(*, '(i5,i3)') head%node%pah%nat, head%hits
            end if
            head => head%next
        end do

!        stop

!        stop
!        print *, pah_array(1)%node%pah%nat, pah_array(final_size)%node%pah%nat

!        allocate(dis(pah_array(1)%node%pah%nat))
!        dis = 0
!        head => database_head
!        do while(associated(head))
!            if (.not. head%node%hasChild) then
!                dis(head%node%pah%nat) = dis(head%node%pah%nat) + head%hits
!            end if
!            head => head%next
!        end do

!        do i=1, pah_array(1)%node%pah%nat
!            if ( dis(i) == 0 ) cycle
!            if (mod(i,2) == 1 ) cycle
!            print *, i, dis(i)
!        end do   
     
!        stop
        
        local_pah_count = final_size / image_count
        local_ext_count = mod(final_size, image_count)
        local_rank = image_count - image_id
        local_calc_count = local_pah_count
        
        if ( local_ext_count > 0 .and. (local_rank <= local_ext_count) ) then   
            local_calc_count = local_calc_count + 1
        end if


        local_index = 0
        j = local_rank
        do i= 1, local_calc_count
            local_index(i) = j
            j = j + image_count
        end do

!    print *, image_id, local_index(:local_calc_count)

        do i = local_calc_count, 1, -1
!            write(*, '(i0, a, i0, a, i0, a, i0)'), image_id, ' running ', i, ' of ' , local_calc_count, 'nat :', pah_array(local_index(i))%node%pah%nat
            call find_ZZ_polynomial(pah_array(local_index(i))%node%pah, level)
        end do

!        print *, image_id, 'local jobs finished'

        
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
!            call clear_tree(root)
        else
            call send_polynomial(pah_array, local_index, local_calc_count)
!            call clear_tree(root)
        end if
        deallocate(pah_array) 
        deallocate(local_index)
    else
        call find_ZZ_polynomial(pah, level)
        call print_ZZ_polynomial(pah)
        call close_file()
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
