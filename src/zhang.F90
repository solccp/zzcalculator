
subroutine print_usage()
    write(*, '(1x,a)') "Usage: ZZ_polynomial input"
    write(*, '(1x,a)') "Options:"
    write(*, '(1x,10a)') "    ", "-K", "              ",  "Compute the number of Kekule structures only"
    write(*, '(1x,10a)') "    ", "-B", "              ",  "Enable built-in bondlist generator"
    write(*, '(1x,10a)') "    ", "-b file", "         ",  "Load user defined bondlist file"
    write(*, '(1x,10a)') "    ", "-n [number]", "     ",  "The number of pretreatment substructures"
    write(*, '(1x,10a)') "    ", "-Q", "              ",  "Print the ZZ polynomial in XML format"

end subroutine

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

#ifdef USE_OPENMP
!=======================================
    integer :: nthreads, max_threads
    integer, external :: OMP_GET_MAX_THREADS
    logical :: thread_set = .false.
#endif
!=======================================

    
    argc = command_argument_count()

    if (argc < 1) then
        call print_usage()
        stop
    end if

    call initialize_options()

    do
#ifdef USE_OPENMP
        okey = getopt('b:n:vBKXQt:')
#else
        okey = getopt('b:n:vBKXQ')
#endif
        if(okey == '>') exit
        if(okey == '!') then
            write(*,*) 'unknown option: ', trim(optarg)
            stop
        end if

        if(okey == 'b') then
            read(optarg, '(a)') options%bondlistfile
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

        if (okey == 'B') then
            options%use_bipartition = .true.
        end if

        if (okey == 'K') then
            options%kekule_only = .true.
        end if

        if (okey == 'X') then
            options%read_connection_table = .true.
        end if

        if (okey == 'Q') then
            options%simple_printing = .true.
        end if
#ifdef USE_OPENMP
        if(okey == 't') then
            nthreads = 0
            read(optarg, *) nthreads
            if (nthreads < 1) then
                nthreads = 1
            end if
            max_threads = OMP_GET_MAX_THREADS()
            if (nthreads > max_threads) then
                nthreads = max_threads
            end if
            thread_set = .true.
            call omp_set_num_threads(nthreads)
        end if
#endif
    end do


#ifdef USE_OPENMP
    if (.not. thread_set ) then
        nthreads = OMP_GET_MAX_THREADS()
        if (nthreads > 1) then
            call omp_set_num_threads(nthreads-1)
        end if       
    end if
#endif


    call cpu_time(start)

    ! ############################################################
    ! # read initial geometry data and create topological matrix #
    ! ############################################################
    if ( options%has_bondlistfile .and. options%verbose ) then
        print *, 'bondlist file: ', trim(options%bondlistfile)
    end if

    call read_input(input_fname, pah)
    call cut_dangling_bonds(pah)



    level = 0

    if ( max_tree_size <= 1) then
        call find_ZZ_polynomial(pah, level)
        if (options%simple_printing) then
            call print_ZZ_polynomial_XML(pah)
        else
            call print_ZZ_polynomial(pah)
        end if
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
        if (options%simple_printing) then
            call print_ZZ_polynomial_XML(pah)
        else
            call print_ZZ_polynomial(pah)
        end if
        call clear_tree(root)
        deallocate(pah_ptr_array)
        deallocate(root)
    end if

    ! #############################################################
    ! # find recursively the ZZ polynomial of the given structure #
    ! #############################################################

    deallocate(ori_geom)
    call cpu_time(finish)

    if ( options%verbose ) then
        write(*, '(a,1x,f12.2,a)' ) 'Total CPU time: ', finish-start, ' s'
    end if



end program
!####################################################################################
!###################### end of program zhang_polynomial #############################
