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
    use hash_m

    implicit none
    integer(kint) :: i, level
    type(structure), target :: pah

!========================================
    integer(kint) :: argc
    character(len=200) :: input_fname
    character :: okey
!========================================
    
    type(tree_node), pointer :: root
    integer :: max_tree_size = 2000
    logical :: reach_limit
    type(database_entry), pointer :: head
    integer :: total_required_strs
!=======================================

!=======================================

    argc = command_argument_count()

    if (argc < 1) then
        stop 'abnormal termination: wrong arguments'
    end if

    call initialize_options()

    do
        okey = getopt('Pfl:b:n:v')
        if(okey == '>') exit
        if(okey == '!') then
            write(*,*) 'unknown option: ', trim(optarg)
            stop
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
            print *, 'bondlist file: ', trim(options%bondlistfile)
            options%has_bondlistfile = .true.
        end if
        if(okey == '.') then
            input_fname = optarg
        end if
        if (okey == 'n') then
            read(optarg, *) max_tree_size
        end if
        if (okey == 'v') then
            options%verbose = .true.
        end if
!        if(okey == 'd') then
!            call load_database_from_file()
!        end if
    end do


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

        total_required_strs = 0
        head => database_head
        do while(associated(head))
            if (.not. head%node%hasChild) then
                total_required_strs = total_required_strs + 1
            end if
            head => head%next
        end do
        head => database_head
        i = 0
        do while(associated(head))
            if (.not. head%node%hasChild) then
                i = i + 1
                if ( options%verbose ) then
                    print *, 'running (', i, '/' , total_required_strs , ')...'
                end if
                call find_ZZ_polynomial(head%node%pah, level)               
            end if
            head => head%next
        end do

        call sum_up(root)


    ! ###########################
    ! # print the ZZ polynomial #
    ! ###########################
        call print_ZZ_polynomial(pah)
        call clear_tree(root)
    else
        call find_ZZ_polynomial(pah, level)
        call print_ZZ_polynomial(pah)
        call close_file()
    end if

!    print *, 'hash: ', pah%hash_key
    ! #############################################################
    ! # find recursively the ZZ polynomial of the given structure #
    ! #############################################################


!    call print_all_database_entry()
!    call save_database_to_file()
!    print *, 'db hit: ', stat_hit
!    print *, 'db no_hit: ', stat_no_hit

    call finalize_temp_space()


end
!####################################################################################
!###################### end of program zhang_polynomial #############################
