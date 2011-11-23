

module tree_m
    use accuracy_m
    use structure_m
    use polynomial_m
    use operator_m
    implicit none

    integer, parameter :: min_limit_nat = 20
    
contains
    subroutine initialize(node)
        type(pah_tree_node), intent(inout), pointer :: node
        nullify(node%pah)
        nullify(node%child_bond)
        nullify(node%child_corners)
        nullify(node%child_ring)
        nullify(node%child_son1)
        nullify(node%child_son2)
    end subroutine
recursive subroutine clear_tree(pah_node)
    type(pah_tree_node), intent(inout), pointer :: pah_node
    logical :: ab, ac, ar, as1, as2

    ab = associated(pah_node%child_bond)
    ac = associated(pah_node%child_corners)
    ar = associated(pah_node%child_ring)
    as1 = associated(pah_node%child_son1)
    as2 = associated(pah_node%child_son2)

    if ( ab .or. ac .or. ar ) then
        if ( ab ) then
            call clear_tree(pah_node%child_bond)
            if ( pah_node%child_bond%shared_count == 0 ) then
                call destory(pah_node%child_bond%pah)
                deallocate(pah_node%child_bond%pah)
                deallocate(pah_node%child_bond)
            else
                pah_node%child_bond%shared_count = pah_node%child_bond%shared_count - 1
            end if
            nullify(pah_node%child_bond)
        end if
        if ( ac ) then
            call clear_tree(pah_node%child_corners)
            if ( pah_node%child_corners%shared_count == 0 ) then
                call destory(pah_node%child_corners%pah)
                deallocate(pah_node%child_corners%pah)
                deallocate(pah_node%child_corners)
            else
                pah_node%child_corners%shared_count = pah_node%child_corners%shared_count - 1
            end if
            nullify(pah_node%child_corners)
        end if
        if (ar) then
            call clear_tree(pah_node%child_ring)
            if ( pah_node%child_ring%shared_count == 0 ) then
                call destory(pah_node%child_ring%pah)
                deallocate(pah_node%child_ring%pah)
                deallocate(pah_node%child_ring)
            else
                pah_node%child_ring%shared_count = pah_node%child_ring%shared_count - 1
            end if
            nullify(pah_node%child_ring)
        end if

    else if ( as1 .or. as2 ) then
        call clear_tree(pah_node%child_son1)
        call clear_tree(pah_node%child_son2)

        if ( pah_node%child_son1%shared_count == 0 ) then
            call destory(pah_node%child_son1%pah)
            deallocate(pah_node%child_son1%pah)
            deallocate(pah_node%child_son1)
        else
            pah_node%child_son1%shared_count = pah_node%child_son1%shared_count - 1
        end if
        nullify(pah_node%child_son1)
        if ( pah_node%child_son2%shared_count == 0 ) then
            call destory(pah_node%child_son2%pah)
            deallocate(pah_node%child_son2%pah)
            deallocate(pah_node%child_son2)
        else
            pah_node%child_son2%shared_count = pah_node%child_son2%shared_count - 1
        end if
        nullify(pah_node%child_son2)
    end if

    
end subroutine
recursive subroutine sum_up(pah_node)
    use database_m

    type(pah_tree_node), intent(inout), pointer :: pah_node
    logical :: ab, ac, ar, as1, as2

    if ( pah_node%pah%polynomial_computed ) then
        return
    end if

    ab = associated(pah_node%child_bond)
    ac = associated(pah_node%child_corners)
    ar = associated(pah_node%child_ring)
    as1 = associated(pah_node%child_son1)
    as2 = associated(pah_node%child_son2)

    if ( ab .or. ac .or. ar ) then
        if ( ab ) then
            call sum_up(pah_node%child_bond)
        end if
        if ( ac ) then
            call sum_up(pah_node%child_corners)
        end if
        if ( ar ) then
            call sum_up(pah_node%child_ring)
            call sum_polynomials(pah_node%pah,pah_node%child_bond%pah, &
            pah_node%child_corners%pah,pah_node%child_ring%pah, .true.)
        else
            call sum_polynomials(pah_node%pah,pah_node%child_bond%pah, &
            pah_node%child_corners%pah,pah_node%child_corners%pah, .false.)
        end if
    else if ( as1 .or. as2 ) then
        if ( as1 ) then
            call sum_up(pah_node%child_son1)
        end if
        if ( as2 ) then
            call sum_up(pah_node%child_son2)
        end if
        call multiply_polynomials(pah_node%pah,pah_node%child_son1%pah, &
        pah_node%child_son2%pah)
    end if

end subroutine
subroutine sort_pah(pah_array, final_tree_size, acc)
    integer, intent(in) :: final_tree_size
    type(pah_tree_node_ptr), dimension(:), intent(inout) :: pah_array
    logical :: acc

    type(pah_tree_node_ptr) :: temp
    integer :: i, j
    
    do i = 1, final_tree_size
        do j = final_tree_size, i + 1, -1
            if (acc) then
                if ( pah_array(j-1)%ptr%pah%nat < pah_array(j)%ptr%pah%nat ) then
                    temp = pah_array(j-1)
                    pah_array(j-1) = pah_array(j)
                    pah_array(j) = temp
                end if
            else
                if ( pah_array(j-1)%ptr%pah%nat > pah_array(j)%ptr%pah%nat ) then
                    temp = pah_array(j-1)
                    pah_array(j-1) = pah_array(j)
                    pah_array(j) = temp
                end if
            end if
        end do
    end do


end subroutine


subroutine build_tree(pah_node, max_tree_size, reach_limit)
    
    use database_m
    use options_m

    type(pah_tree_node), intent(inout), pointer :: pah_node
    integer, intent(in) :: max_tree_size
    logical, intent(out) :: reach_limit

    type(pah_tree_node), pointer :: cur_node, tmp_node
    type(pah_tree_node), pointer :: bond, corners, ring


    logical :: ring_exists, are_neighbors
    integer :: atom1,atom2,i,nelim
    integer, dimension(6) :: sextet
    integer, dimension(2) :: atoms
    integer :: medat

    logical :: hit


    NULLIFY(ring)

    reach_limit = .false.

    call add_database_entry(pah_node, hit, tmp_node)

outer: do
        if (reach_limit) then
            exit
        end if

        call get_database_first_entry(cur_node)
        if ( .not. associated(cur_node) ) then
            if ( options%verbose ) then
                print *, 'empty'
            end if
            exit
        end if

        if ( cur_node%hasChild ) then
            print *, 'some logic error'
            stop
        end if

        if ( cur_node%pah%nat == 0 ) then
            call set_polynomial(cur_node%pah, 1_kint)
            cycle
        else if ( mod(cur_node%pah%nat,2) /= 0 ) then
            call set_polynomial(cur_node%pah, 0_kint)
            cycle
        else
            call check_if_connected(cur_node%pah,medat)
            cur_node%hasChild = .true.
            if ( medat == 0 ) then

                call select_edge_bond(cur_node%pah,atom1,atom2)
                ring_exists = .false.
                if ( .not. options%kekule_only ) then
                    call find_edge_ring(cur_node%pah,sextet,atom1,atom2,ring_exists)
                end if
       
                allocate(bond)
                call initialize(bond)
                allocate(bond%pah)
                call create_nobond_daughter(cur_node%pah,bond%pah,atom1,atom2)
                call cut_dangling_bonds(bond%pah)

                if ( bond%pah%nat == 0 ) then
                    call set_polynomial(bond%pah, 1_kint)
                else if ( bond%pah%nat < 6 ) then
                    call set_polynomial(bond%pah, 0_kint)
                else
                    call add_database_entry(bond, hit, tmp_node)
                    if ( hit ) then
                        call clear_tree(tmp_node)
                        deallocate(tmp_node)
                    end if
                end if
                cur_node%child_bond => bond
                if ( associated(cur_node%child_bond, cur_node) ) then
                    print *, 'wtf'
                    stop
                end if

                allocate(corners)
                call initialize(corners)
                allocate(corners%pah)
                nelim = 2
                atoms(1) = atom1
                atoms(2) = atom2
                call create_noatoms_daughter(cur_node%pah,corners%pah,nelim,atoms,.false.)
                call cut_dangling_bonds(corners%pah)
                    
                if ( corners%pah%nat == 0 ) then
                    call set_polynomial(corners%pah, 1_kint)
                else if ( corners%pah%nat < 6 ) then
                    call set_polynomial(corners%pah, 0_kint)
                else
                    call add_database_entry(corners, hit, tmp_node)
                    if ( hit ) then
                        call clear_tree(tmp_node)
                        deallocate(tmp_node)
                    end if
                end if
                cur_node%child_corners => corners

                if (ring_exists) then
                    allocate(ring)
                    call initialize(ring)
                    allocate(ring%pah)
                    nelim = 6
                    call create_noatoms_daughter(cur_node%pah,ring%pah,nelim,sextet,.true.)
                    call cut_dangling_bonds(ring%pah)
                    
                    if ( ring%pah%nat == 0 ) then
                        call set_polynomial(ring%pah, 1_kint)
                    else if ( ring%pah%nat < 6 ) then
                        call set_polynomial(ring%pah, 0_kint)
                    else
                        call add_database_entry(ring, hit, tmp_node)
                        if ( hit ) then
                            call clear_tree(tmp_node)
                            deallocate(tmp_node)
                        end if
                    end if
                    cur_node%child_ring => ring
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
                    call set_polynomial(cur_node%pah, 0_kint)
                else
                    if ( bond%pah%nat == 0 ) then
                        call set_polynomial(bond%pah, 1_kint)
                    else if ( bond%pah%nat < 6 ) then
                        call set_polynomial(bond%pah, 0_kint)
                    else
                        call add_database_entry(bond, hit, tmp_node)
                        if ( hit ) then
                            call clear_tree(tmp_node)
                            deallocate(tmp_node)
                        end if
                    end if
                    if ( corners%pah%nat == 0 ) then
                        call set_polynomial(corners%pah, 1_kint)
                    else if ( corners%pah%nat < 6 ) then
                        call set_polynomial(corners%pah, 0_kint)
                    else
                        call add_database_entry(corners, hit, tmp_node)
                        if ( hit ) then
                            call clear_tree(tmp_node)
                            deallocate(tmp_node)
                        end if
                    end if
                    cur_node%child_son1 => bond
                    cur_node%child_son2 => corners
                end if
            end if
        end if
        if ( database_size >= max_tree_size ) then
            reach_limit = .true.
        end if
        if ( cur_node%pah%nat < min_limit_nat )  then
            reach_limit = .true.
        end if 
    end do outer   
end subroutine 
end module
