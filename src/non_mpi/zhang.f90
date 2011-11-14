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
    use accuracy_m
    use structure_m
    use options_m
    use getopt_m
    use output_m
    use database_m
    use tree_m
    use input_m
    use decompose_print_m

    implicit none
    integer :: i, level
    type(structure), target :: pah

    real :: start, finish

!========================================
    integer :: argc
    character(len=200) :: input_fname
    character :: okey
    character(len=200) :: dummy
!========================================
    
    type(pah_tree_node), pointer :: root
    integer :: max_tree_size = 0
    logical :: reach_limit
    integer :: total_required_strs
    type(rb_treenode_ptr) :: node_ptr
    type(structure_ptr), allocatable, dimension(:) :: pah_ptr_array
!=======================================

    integer :: max_nat
    real :: avg_nat
    character(len=200) :: mem_used
    
!=======================================

    argc = command_argument_count()

    if (argc < 1) then
        stop 'abnormal termination: wrong arguments'
    end if

    call initialize_options()

    do
        okey = getopt('Pfl:b:n:vtp:c:')
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
            options%has_bondlistfile = .true.
        end if
        if(okey == '.') then
            input_fname = optarg
        end if
        if (okey == 'n') then
!            read(optarg, *) i
!            max_tree_size = max(max_tree_size, i)
            read(optarg, *) max_tree_size
        end if
        if (okey == 'v') then
            options%verbose = .true.
        end if

        if (okey == 't') then
            options%testrun = .true.
        end if
        if (okey == 'p') then
            read(optarg, *) i
            options%decompose_print = .true.
            call init_decompose_print('substructures.cml', i)
        end if

        if (okey == 'c') then
            read(optarg, '(a)') dummy
            options%use_connection_file = .true.
            options%connection_file = trim(dummy)
        end if

!        if(okey == 'd') then
!            read(optarg, '(a)') dummy
!            call load_database_from_file(len(trim(dummy)), dummy)
!            options%use_database = .true.
!        end if
!        if(okey == 'c') then
!            read(optarg, '(a)') options%databasefile
!            if (len(trim(options%databasefile)) == 0) then
!                print *, 'database filename is not valid'
!                stop
!            end if
!            options%create_database = .true.
!            options%use_database = .true.
!        end if

    end do


    call cpu_time(start)

    ! ############################################################
    ! # read initial geometry data and create topological matrix #
    ! ############################################################
    if ( options%has_bondlistfile .and. options%verbose ) then
        print *, 'bondlist file: ', trim(options%bondlistfile)
    end if

    call read_input(input_fname, pah)
    call cut_dangling_bonds(pah)



    if ( options%print_intermediate_structures .and. pah%nat > 50 .and. .not. options%force_print_structures ) then
        print *, 'warning: # of atoms > 50, disabling intermediate structure printing'
        print *, '  use -P -f option to print the intermediate structures forcely'
        options%print_intermediate_structures = .false.
    end if


    level = 0

    if ( options%decompose_print ) then
        allocate(root)
        root%pah => pah
        call write_decomposed_substructures(root)
        call clear_tree(root)
        deallocate(root)
    else if ( options%print_intermediate_structures .or. max_tree_size <= 1) then
        call find_ZZ_polynomial(pah, level)
        call print_ZZ_polynomial(pah)
        call close_file()
    else
        allocate(root)
        root%pah => pah
        call set_max_size(int(maxval(pah%indexmapping)))
        call build_tree(root, max_tree_size, reach_limit)

        total_required_strs = database_size

        if (options%verbose) then
            write(*,'(a, 1x, i0)') 'Number of decomposed unique substructures:', total_required_strs
        end if

        allocate(pah_ptr_array(total_required_strs))
        
        avg_nat = 0.0
        max_nat = 0

        i = 0
        do while( database_size > 0)
            i = i + 1
            call getMin(tree_child, node_ptr)
            if ( .not. associated(node_ptr%ptr) ) then 
                total_required_strs = i - 1
                exit
            end if

            pah_ptr_array(i)%ptr => node_ptr%ptr%keyvalue%value%pah
            if ( pah_ptr_array(i)%ptr%nat > max_nat ) then
                max_nat = pah_ptr_array(i)%ptr%nat
            end if
            avg_nat = avg_nat + real(pah_ptr_array(i)%ptr%nat)
!            write(100,'(i3,1x,a)') node_ptr%ptr%keyvalue%ckey%key1, char(node_ptr%ptr%keyvalue%ckey%key2)
            call erase_node(tree_child, node_ptr%ptr)
        end do

        if ( options%testrun ) then
            write(*,'(a, 1x, i0)') 'Number of atoms in structure:', root%pah%nat
            write(*,'(a, 1x, i0)') 'Size of Decomposed structure database:', total_required_strs
            write(*,'(a, 1x, i0, 1x, i0)') 'Intermediate structures and number of hits:', total_str, total_hit
            write(*,'(a, 1x, i0, 1x, i0)') 'Maximum and average number of atoms in decomposed database:', max_nat, &
            int(avg_nat/real(total_required_strs))
        else
!$OMP parallel default(shared) private(i, level) if ( total_required_strs > 500 .or. pah%nat > 300 )
!$OMP DO SCHEDULE(GUIDED)
            do i = 1, total_required_strs
                if ( options%verbose ) then
                    print *, 'running (', i, '/' , total_required_strs , ')...  (', pah_ptr_array(i)%ptr%nat , ')'
                end if
                call find_ZZ_polynomial(pah_ptr_array(i)%ptr, level)
                if ( options%verbose ) then
                    call print_ZZ_polynomial(pah_ptr_array(i)%ptr)
                end if
            end do
!$OMP END DO
!$OMP END parallel

            call sum_up(root)

            ! ###########################
            ! # print the ZZ polynomial #
            ! ###########################
            call print_ZZ_polynomial(pah)
        end if            
        call clear_tree(root)
        deallocate(pah_ptr_array)
        deallocate(root)
    end if

    ! #############################################################
    ! # find recursively the ZZ polynomial of the given structure #
    ! #############################################################


    call cpu_time(finish)

    if ( options%verbose ) then
        write(*, '(a,1x,f12.2,a)' ) 'Total CPU time: ', finish-start, ' s'
    end if



end program
!####################################################################################
!###################### end of program zhang_polynomial #############################
