

module tree
    use types_module
    use structure_module
    implicit none

    integer, parameter :: min_limit_nat = 10
    
contains
recursive subroutine clear_tree(pah_node)
    use types_module
    use structure_module
    type(tree_node), intent(inout), pointer :: pah_node
    logical :: ab, ac, ar, as1, as2

    ab = associated(pah_node%child_bond)
    ac = associated(pah_node%child_corners)
    ar = associated(pah_node%child_ring)
    as1 = associated(pah_node%child_son1)
    as2 = associated(pah_node%child_son2)

    if ( ab .or. ac .or. ar ) then
        if ( ab ) then
            call clear_tree(pah_node%child_bond)
        end if
        if ( ac ) then
            call clear_tree(pah_node%child_corners)
        end if
        if (ar) then
            call clear_tree(pah_node%child_ring)
        end if
        call destory(pah_node%child_bond%pah)
        call destory(pah_node%child_corners%pah)

        deallocate(pah_node%child_bond%pah)
        deallocate(pah_node%child_corners%pah)

        deallocate(pah_node%child_bond, pah_node%child_corners)


        nullify(pah_node%child_bond)
        nullify(pah_node%child_corners)

        if (ar) then
            call destory(pah_node%child_ring%pah)
            deallocate(pah_node%child_ring%pah)
            deallocate(pah_node%child_ring)
            nullify(pah_node%child_ring)
        end if
    else if ( as1 .or. as2 ) then               
        call clear_tree(pah_node%child_son1)
        call clear_tree(pah_node%child_son2)
        call destory(pah_node%child_son1%pah)
        call destory(pah_node%child_son2%pah)
        deallocate(pah_node%child_son1%pah)
        deallocate(pah_node%child_son2%pah)
        nullify(pah_node%child_son1)
        nullify(pah_node%child_son2)
    end if

    
end subroutine
recursive subroutine sum_up(pah_node)
    use types_module
    use structure_module
    use database_m

    type(tree_node), intent(inout), pointer :: pah_node
    logical :: ab, ac, ar, as1, as2


    ab = associated(pah_node%child_bond)
    ac = associated(pah_node%child_corners)
    ar = associated(pah_node%child_ring)
    as1 = associated(pah_node%child_son1)
    as2 = associated(pah_node%child_son2)

!    print *, char(pah_node%key)


    if ( ab .or. ac .or. ar ) then
        if ( ab ) then
            call sum_up(pah_node%child_bond)
        end if
        if ( ac ) then
            call sum_up(pah_node%child_corners)
        end if
        if ( ar ) then
            call sum_up(pah_node%child_ring)
            call sum_polynomials(pah_node%pah,pah_node%child_bond%pah,pah_node%child_corners%pah,pah_node%child_ring%pah, .true.)
        else
            call sum_polynomials(pah_node%pah,pah_node%child_bond%pah,pah_node%child_corners%pah,pah_node%child_corners%pah, .false.)
        end if
        nullify(pah_node%child_bond)
        nullify(pah_node%child_corners)
        if ( ar ) then
            nullify(pah_node%child_ring)
        end if
    else if ( as1 .or. as2 ) then
        if ( as1 ) then
            call sum_up(pah_node%child_son1)
        end if
        if ( as2 ) then
            call sum_up(pah_node%child_son2)
        end if
        call multiply_polynomials(pah_node%pah,pah_node%child_son1%pah,pah_node%child_son2%pah)
        nullify(pah_node%child_son1)
        nullify(pah_node%child_son2)
    end if

    
end subroutine
recursive subroutine visit_tree(pah_node)
    use types_module
    use structure_module
    use database_m

    type(tree_node), intent(inout), pointer :: pah_node
    logical :: ab, ac, ar, as1, as2


    ab = associated(pah_node%child_bond)
    ac = associated(pah_node%child_corners)
    ar = associated(pah_node%child_ring)
    as1 = associated(pah_node%child_son1)
    as2 = associated(pah_node%child_son2)


    if ( ab .or. ac .or. ar ) then
        if ( ab ) then
            call visit_tree(pah_node%child_bond)
        end if
        if ( ac ) then
            call visit_tree(pah_node%child_corners)
        end if
        if ( ar ) then
            call visit_tree(pah_node%child_ring)
        end if
!        nullify(pah_node%child_bond)
!        nullify(pah_node%child_corners)
!        if ( ar ) then
!            nullify(pah_node%child_ring)
!        end if
    else if ( as1 .or. as2 ) then
        if ( as1 ) then
            call visit_tree(pah_node%child_son1)
        end if
        if ( as2 ) then
            call visit_tree(pah_node%child_son2)
        end if
!        nullify(pah_node%child_son1)
!        nullify(pah_node%child_son2)
    end if

    
end subroutine
subroutine sort_pah(pah_array, final_tree_size, acc)
    integer, intent(in) :: final_tree_size
    type(tree_node_ptr), dimension(:), intent(inout) :: pah_array
    logical :: acc

    type(tree_node_ptr) :: temp
    integer :: i, j
    
    do i = 1, final_tree_size
        do j = final_tree_size, i + 1, -1
            if (acc) then
                if ( pah_array(j-1)%node%pah%nat < pah_array(j)%node%pah%nat ) then
                    temp = pah_array(j-1)
                    pah_array(j-1) = pah_array(j)
                    pah_array(j) = temp
                end if
            else
                if ( pah_array(j-1)%node%pah%nat > pah_array(j)%node%pah%nat ) then
                    temp = pah_array(j-1)
                    pah_array(j-1) = pah_array(j)
                    pah_array(j) = temp
                end if
            end if
        end do
    end do


end subroutine

subroutine set_polynomial(pah, value)
    type(structure), intent(inout) :: pah
    integer(kint), intent(in) :: value
    
    pah%order = 0
    allocate(pah%polynomial(1))
    pah%polynomial(1)=setvli(value)
    pah%polynomial_computed = .true.
end subroutine

subroutine initialize(node)
    type(tree_node), intent(inout), pointer :: node
    nullify(node%pah)
    nullify(node%child_bond)
    nullify(node%child_corners)
    nullify(node%child_ring)
    nullify(node%child_son1)
    nullify(node%child_son2)
end subroutine
subroutine build_tree(pah_node, max_tree_size, reach_limit)
    
    use hash_m
    use database_m

    type(tree_node), intent(inout), pointer :: pah_node
    integer, intent(in) :: max_tree_size
    logical, intent(out) :: reach_limit

    type(tree_node), pointer :: cur_node 
    type(tree_node), pointer :: bond, corners, ring


    logical :: ring_exists, are_neighbors
    integer(kint) :: atom1,atom2,i,nelim
    integer(kint), dimension(6) :: sextet
    integer(kint), dimension(2) :: atoms
    integer(kint) :: medat

    logical :: hit
    type(database_entry), pointer :: head


    NULLIFY(ring)

    reach_limit = .false.

!    call get_hash(pah_node%pah, pah_node%pah%hash_key) 
    call add_database_entry(pah_node, hit)

    outer: do
        if (reach_limit) then
            exit
        end if

        call get_database_first_entry(cur_node)
        if ( .not. associated(cur_node) ) then
!            print *, 'all str'
            exit
        end if

!        call get_polynomial(cur_node%pah, hit)
!        if (hit) then
!            print *, cur_node%pah%hash_key
!            print *, 'hit'
!            cycle
!        end if


        if ( cur_node%hasChild ) then
            print *, 'some logic error'
        end if

        if ( cur_node%pah%nat == 0 ) then
            call set_polynomial(cur_node%pah, 1_kint)
        else
            call check_if_connected(cur_node%pah,medat)
            cur_node%hasChild = .true.
            database_size = database_size -1
            if ( medat == 0 ) then

                call select_edge_bond(cur_node%pah,atom1,atom2)
                call find_edge_ring(cur_node%pah,sextet,atom1,atom2,ring_exists)
       
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
                    
!                    call get_hash(bond%pah, bond%pah%hash_key) 
                    call add_database_entry(bond, hit)
                end if
                cur_node%child_bond => bond

                allocate(corners)
                call initialize(corners)
                allocate(corners%pah)
                nelim=2
                atoms(1)=atom1
                atoms(2)=atom2
                call create_noatoms_daughter(cur_node%pah,corners%pah,nelim,atoms,.false.)
                call cut_dangling_bonds(corners%pah)
                    
                if ( corners%pah%nat == 0 ) then
                    call set_polynomial(corners%pah, 1_kint)
                else if ( corners%pah%nat < 6 ) then
                    call set_polynomial(corners%pah, 0_kint)
                else
!                    call get_hash(corners%pah, corners%pah%hash_key) 
                    call add_database_entry(corners, hit)
                end if
                cur_node%child_corners => corners

                if (ring_exists) then
                    allocate(ring)
                    call initialize(ring)
                    allocate(ring%pah)
                    nelim=6
                    call create_noatoms_daughter(cur_node%pah,ring%pah,nelim,sextet,.true.)
                    call cut_dangling_bonds(ring%pah)

                    
                    if ( ring%pah%nat == 0 ) then
                        call set_polynomial(ring%pah, 1_kint)
                    else if ( ring%pah%nat < 6 ) then
                        call set_polynomial(ring%pah, 0_kint)
                    else
!                        call get_hash(ring%pah, ring%pah%hash_key) 
                        call add_database_entry(ring, hit)
                    end if
                    cur_node%child_ring => ring
                end if
            else
!                print *, 'found a disconnected pah, splitting into sons'
                allocate(bond)
                allocate(corners)
                call initialize(corners)
                call initialize(bond)
                allocate(bond%pah)
                allocate(corners%pah)
                call split_structure(cur_node%pah, bond%pah, corners%pah, medat)
                if ( mod(bond%pah%nat, 2_kint) == 1) then
!                    print *, 'found invalid splitting to odd nats'
                    call destory(bond%pah)
                    call destory(corners%pah)
                    deallocate(bond%pah)
                    deallocate(corners%pah)
                    deallocate(bond)
                    deallocate(corners)
                    call set_polynomial(cur_node%pah, 0_kint)
                else
!                    print *, 'found valid splitting to even nats'
                    if ( bond%pah%nat == 0 ) then
                        call set_polynomial(bond%pah, 1_kint)
                    else
!                        call get_hash(bond%pah, bond%pah%hash_key) 
                        call add_database_entry(bond, hit)
                    end if
                    if ( corners%pah%nat == 0 ) then
                        call set_polynomial(corners%pah, 1_kint)
                    else
!                        call get_hash(corners%pah, corners%pah%hash_key) 
                        call add_database_entry(corners, hit)
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
!            print *, 'reach limit'
!            cycle
        end if 
!        write (*, '(a)') "=================================="
!        head => database_head
!        do while(associated(head))
!            write(*, '(a,i3)') char(head%key), head%hits
!            head => head%next
!        end do
!        write (*, '(a)') "=================================="
    end do outer   
end subroutine 
end module
