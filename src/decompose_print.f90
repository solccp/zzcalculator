module decompose_print_m
    use output_m
    use structure_m
    implicit none
    save
    integer :: unitnum, level_to_print

contains

    subroutine init_decompose_print(filename, level_print)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: level_print
        unitnum = getunit()
        open(file=trim(filename), unit = unitnum) 
        level_to_print = level_print
    end subroutine


    subroutine get_path_str(cur_node, str, cur_level)
        type(pah_tree_node), intent(in), pointer :: cur_node
        type(pah_tree_node), pointer :: node
        integer, intent(in) :: cur_level
        character(len=*) :: str
        character(len=cur_level) :: str_tmp
        integer :: pos, i
        pos = 1
        node => cur_node
        do while( associated(node%parent) )
            if ( associated(node%parent%child_bond, node) ) then
                str_tmp(pos:pos) = 'S'
            else if ( associated(node%parent%child_corners, node) ) then
                str_tmp(pos:pos) = 'D'
            else if ( associated(node%parent%child_ring, node) ) then
                str_tmp(pos:pos) = 'R'
            else if ( associated(node%parent%child_son1, node) ) then
                str_tmp(pos:pos) = '1'
            else if ( associated(node%parent%child_son2, node) ) then
                str_tmp(pos:pos) = '2'
            end if
            pos = pos + 1
            node => node%parent
        end do
        do i = 1, pos - 1
            str(i:i) = str_tmp(pos-i:pos-i)
        end do

    end subroutine

    subroutine finalize_decompose_print() 
        close(unitnum)
        unitnum = -1
    end subroutine

    subroutine write_decomposed_substructures(pah_node)
        
        use options_m
        use tree_m

        type(pah_tree_node), intent(inout), pointer :: pah_node

        type(pah_tree_node), pointer :: cur_node, tmp_node
        type(pah_tree_node), pointer :: bond, corners, ring


        type(pah_tree_node_ptr), dimension(:), allocatable :: queue
        type(pah_tree_node_ptr), dimension(:), allocatable :: queue_tmp


        integer :: index_q1, index_qt

        logical :: ring_exists, are_neighbors
        integer :: atom1,atom2,i,nelim
        integer, dimension(6) :: sextet
        integer, dimension(2) :: atoms
        integer :: medat, level

        logical :: hit
        character(len=level_to_print) :: path

!        NULLIFY(ring)


        allocate(queue(3)) 
        allocate(queue_tmp(10)) 
        index_q1 = 0
        index_qt = 0

        level = 0
        index_q1 = index_q1 + 1
        queue(index_q1)%ptr => pah_node

        do while( level < level_to_print)
            level = level + 1

            if ( size(queue_tmp) < index_q1*3 ) then
                deallocate(queue_tmp)
                allocate(queue_tmp(index_q1*3))
            end if
            do while( index_q1 > 0 )

                cur_node => queue(index_q1)%ptr
                index_q1 = index_q1 - 1

                if ( cur_node%hasChild ) then
                    print *, 'some logic error'
                    stop
                end if

                if ( cur_node%pah%nat /= 0 ) then
                    call check_if_connected(cur_node%pah,medat)
                    cur_node%hasChild = .true.
                    if ( medat == 0 ) then

                        call select_edge_bond(cur_node%pah,atom1,atom2)
                        call find_edge_ring(cur_node%pah,sextet,atom1,atom2,ring_exists)
               
                        allocate(bond)
                        call initialize(bond)
                        allocate(bond%pah)
                        call create_nobond_daughter(cur_node%pah,bond%pah,atom1,atom2)
                        call cut_dangling_bonds(bond%pah)

                        if ( bond%pah%nat >= 6 ) then
                            index_qt = index_qt + 1
                            queue_tmp(index_qt)%ptr => bond
                        end if
                        if ( bond%pah%nat == 2 ) then
                           call clear_tree(bond) 
                        else
                            cur_node%child_bond => bond
                            bond%parent => cur_node
                        end if

                        allocate(corners)
                        call initialize(corners)
                        allocate(corners%pah)
                        nelim = 2
                        atoms(1)=atom1
                        atoms(2)=atom2
                        call create_noatoms_daughter(cur_node%pah,corners%pah,nelim,atoms,.false.)
                        call cut_dangling_bonds(corners%pah)
                            
                        if ( corners%pah%nat >= 6 ) then
                            index_qt = index_qt + 1
                            queue_tmp(index_qt)%ptr => corners
                        end if
                        if ( corners%pah%nat == 2 ) then
                           call clear_tree(corners)
                        else
                            cur_node%child_corners => corners
                            corners%parent => cur_node
                        end if

                        if (ring_exists) then
                            allocate(ring)
                            call initialize(ring)
                            allocate(ring%pah)
                            nelim = 6
                            call create_noatoms_daughter(cur_node%pah,ring%pah,nelim,sextet,.true.)
                            call cut_dangling_bonds(ring%pah)
                            
                            if ( ring%pah%nat >= 6 ) then
                                index_qt = index_qt + 1
                                queue_tmp(index_qt)%ptr => ring
                            end if
                            if ( ring%pah%nat == 2 ) then
                                call clear_tree(ring)
                            else
                                cur_node%child_ring => ring
                                ring%parent => cur_node
                            end if
                        end if
                    else
                        allocate(bond)
                        allocate(corners)
                        call initialize(corners)
                        call initialize(bond)
                        allocate(bond%pah)
                        allocate(corners%pah)
                        call split_structure(cur_node%pah, bond%pah, corners%pah, medat)
                        if ( mod(bond%pah%nat, 2) == 1) then
                            call clear_tree(bond)
                            call clear_tree(corners)
                        else
                            if ( bond%pah%nat >= 6 ) then
                                index_qt = index_qt + 1
                                queue_tmp(index_qt)%ptr => bond
                            end if
                            if ( corners%pah%nat >= 6 ) then
                                index_qt = index_qt + 1
                                queue_tmp(index_qt)%ptr => corners
                            end if
                            cur_node%child_son1 => bond
                            cur_node%child_son2 => corners
                            bond%parent => cur_node
                            corners%parent => cur_node
                        end if
                    end if
                end if
            end do
            if ( index_qt == 0 ) then
                exit
            end if
            if ( size(queue) /= size(queue_tmp) ) then
                deallocate(queue)
                allocate(queue(size(queue_tmp)))
            end if
            index_q1 = index_qt
            index_qt = 0
            queue(:index_q1) = queue_tmp(:index_q1)
        end do

!        if (index_q1 > 0 ) then
            call write_cml_header_start(unitnum)
!            do i=1, index_q1
                !call get_path_str(queue(i)%ptr, path, level_to_print)
                call print_structures(pah_node, 0, path)
!                write_cml_unit(queue(i)%ptr%pah, unitnum, path)
!            end do
            call write_cml_header_end(unitnum)
!        end if
     
    end subroutine 

    recursive subroutine print_structures(node, level, path)
        type(pah_tree_node), intent(in) :: node
        character(len=*), intent(inout) :: path
        integer, intent(in) :: level

        if ( associated(node%child_bond) ) then
            path(level+1:level+1) = 'S'
            call print_structures(node%child_bond, level+1, path)
        end if
        if ( associated(node%child_corners) ) then
            path(level+1:level+1) = 'D'
            call print_structures(node%child_corners, level+1, path)
        end if
        if ( associated(node%child_ring) ) then
            path(level+1:level+1) = 'R'
            call print_structures(node%child_ring, level+1, path)
        end if
        if ( associated(node%child_son1) ) then
            path(level+1:level+1) = '1'
            call print_structures(node%child_son1, level+1, path)
        end if
        if ( associated(node%child_son2) ) then
            path(level+1:level+1) = '2'
            call print_structures(node%child_son2, level+1, path)
        end if
        if ( .not. associated(node%child_bond) .and. .not. & 
                associated(node%child_corners) .and. &
             .not. associated(node%child_ring) .and. .not. &
                associated(node%child_son1) .and. .not. associated(node%child_son2) ) then

!            print *, node%pah%nat
            call write_cml_unit(node%pah, unitnum, path)
        end if


    end subroutine

    subroutine write_cml_header_start(unitnum)
        integer, intent(in) :: unitnum
        write(unitnum, '(a)') '<?xml version="1.0"?>'
        write(unitnum,'(a)') '<cml xmlns="http://www.xml-cml.org/schema">'
    end subroutine
    subroutine write_cml_header_end(unitnum)
        integer, intent(in) :: unitnum
        write(unitnum,'(a)') '</cml>'
    end subroutine
    subroutine write_cml_unit(pah, unitnum, title)
        use structure_m
        type(structure), intent(in) :: pah
        integer, intent(in) :: unitnum
        character(len=*), intent(in) :: title
        integer :: i,j
!        logical :: res

!        res = .true.

        write(unitnum, '(2x,a)') '<molecule>'
        write(unitnum, '(4x,a,a,a)') '<comment>', trim(title), '</comment>'
        if ( pah%nat > 0) then
            write(unitnum, '(4x,a)') '<atomArray>'
            do i=1, pah%nat
                if ( pah%neighbornumber(i) == 0 ) then
                    cycle
                end if
                write(unitnum, '(6x,a,i0,a,F20.12,a,F20.12,a,F20.12,a)') &
                '<atom id="a', pah%indexmapping(i),'" elementType="C" x3="', &
                ori_geom(1,pah%indexmapping(i)),'" y3="', &
                ori_geom(2,pah%indexmapping(i)),'" z3="', &
                ori_geom(3,pah%indexmapping(i)), '"/>'
            end do
            write(unitnum, '(4x,a)') '</atomArray>'
            write(unitnum, '(4x,a)') '<bondArray>'
            do i=1, pah%nat
                do j = 1, pah%neighbornumber(i)
                    write(unitnum, '(6x,a,i0,a,i0,a)') '<bond atomRefs2="a', &
                    pah%indexmapping(i),' a', pah%indexmapping(pah%neighborlist(i,j)),'" order="1"/>'
                end do
            end do
            write(unitnum, '(4x,a)') '</bondArray>'
        else
            write(unitnum, '(4x,a)') '<atomArray>'
            write(unitnum, '(4x,a)') '</atomArray>'
        end if
        write(unitnum, '(2x,a)') '</molecule>'

    end subroutine

end module
