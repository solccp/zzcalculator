

module tree
    use types_module
    use structure_module
    implicit none

    integer, parameter :: min_limit_nat = 20
    integer, parameter :: str_per_core = 25
    
contains
recursive subroutine clear_tree(pah)
    use types_module
    use structure_module
    type(structure), intent(inout) :: pah
    logical :: ab, ac, ar, as1, as2

    ab = associated(pah%child_bond)
    ac = associated(pah%child_corners)
    ar = associated(pah%child_ring)
    as1 = associated(pah%child_son1)
    as2 = associated(pah%child_son2)

    if ( ab .and. ac ) then
        call clear_tree(pah%child_bond)
        call clear_tree(pah%child_corners)
        if (ar) then
            call clear_tree(pah%child_ring)
        end if

        call destory(pah%child_bond)
        call destory(pah%child_corners)

        deallocate(pah%child_bond)
        deallocate(pah%child_corners)

        nullify(pah%child_bond)
        nullify(pah%child_corners)

        if (ar) then
            call destory(pah%child_ring)
            deallocate(pah%child_ring)
            nullify(pah%child_ring)
        end if
    else if ( as1 .and. as2 ) then
                
        call clear_tree(pah%child_son1)
        call clear_tree(pah%child_son2)
        call destory(pah%child_son1)
        call destory(pah%child_son2)
        nullify(pah%child_son1)
        nullify(pah%child_son2)
    end if

    
end subroutine
recursive subroutine sum_up(pah)
    use types_module
    use structure_module
    type(structure), intent(inout) :: pah
    logical :: ab, ac, ar, as1, as2

    ab = associated(pah%child_bond)
    ac = associated(pah%child_corners)
    ar = associated(pah%child_ring)
    as1 = associated(pah%child_son1)
    as2 = associated(pah%child_son2)



    if ( ab .and. ac ) then
        call sum_up(pah%child_bond)
        call sum_up(pah%child_corners)
        if (ar) then
            call sum_up(pah%child_ring)
        end if
        call sum_polynomials(pah,pah%child_bond,pah%child_corners,pah%child_ring, ar)

        call destory(pah%child_bond)
        call destory(pah%child_corners)

        deallocate(pah%child_bond)
        deallocate(pah%child_corners)

        nullify(pah%child_bond)
        nullify(pah%child_corners)

        if (ar) then
            call destory(pah%child_ring)
            deallocate(pah%child_ring)
            nullify(pah%child_ring)
        end if
    else if ( as1 .and. as2 ) then
        call sum_up(pah%child_son1)
        call sum_up(pah%child_son2)
        call multiply_polynomials(pah,pah%child_son1,pah%child_son2)
        call destory(pah%child_son1)
        call destory(pah%child_son2)
        deallocate(pah%child_son1)
        deallocate(pah%child_son2)

        nullify(pah%child_son1)
        nullify(pah%child_son2)
    end if

    
end subroutine
subroutine sort_pah(pah_array, final_tree_size, acc)
    integer, intent(in) :: final_tree_size
    type(pah_ptr), dimension(:), intent(inout) :: pah_array
    logical :: acc

    type(pah_ptr) :: temp
    integer :: i, j
    
    do i = 1, final_tree_size
        do j = final_tree_size, i + 1, -1
            if (acc) then
                if ( pah_array(j-1)%ptr%nat < pah_array(j)%ptr%nat ) then
                    temp = pah_array(j-1)
                    pah_array(j-1) = pah_array(j)
                    pah_array(j) = temp
                end if
            else
                if ( pah_array(j-1)%ptr%nat > pah_array(j)%ptr%nat ) then
                    temp = pah_array(j-1)
                    pah_array(j-1) = pah_array(j)
                    pah_array(j) = temp
                end if
            end if
        end do
    end do


end subroutine
subroutine build_tree(pah, image_count, max_tree_size, pah_array, final_tree_size, reach_limit)
    
    type(structure), target, intent(inout) :: pah
    integer, intent(in) :: image_count
    integer, intent(in) :: max_tree_size
    integer, intent(out) :: final_tree_size
    type(pah_ptr), dimension(max_tree_size), intent(inout) :: pah_array
    logical, intent(out) :: reach_limit

    
    integer :: counter1, counter2
    type(structure), pointer :: cur_node 
    type(structure), pointer :: bond, corners, ring
    type(pah_ptr), dimension(max_tree_size) :: t_pah_array


    logical :: ring_exists, are_neighbors
    integer(kint) :: atom1,atom2,i,nelim
    integer(kint), dimension(6) :: sextet
    integer(kint), dimension(2) :: atoms
    integer(kint) :: medat

    NULLIFY(ring)

    reach_limit = .false.

    counter1 = 1
    pah_array(1)%ptr => pah   

    if ( image_count == 1 ) then
        final_tree_size = 1
        return
    end if
 
!    if ( counter1 > 0 .and. mod(counter1, image_count) == 0 ) then
!        final_tree_size = counter1 
!        return
!    end if
   
    outer: do
        if (reach_limit) then
            exit
        end if
        counter2 = 0
        do while(counter1 > 0)
            cur_node => pah_array(counter1)%ptr
            if ( cur_node%nat < min_limit_nat )  then
                reach_limit = .true.
            end if 
            counter1 = counter1 - 1

            if ( cur_node%nat > 0 ) then

                call check_if_connected(cur_node,medat)
                if ( medat == 0 ) then

                    call select_edge_bond(cur_node,atom1,atom2)
                    call find_edge_ring(cur_node,sextet,atom1,atom2,ring_exists)
       
                    allocate(bond)
                    call create_nobond_daughter(cur_node,bond,atom1,atom2)
                    call cut_dangling_bonds(bond)
                    counter2 = counter2 + 1
                    t_pah_array(counter2)%ptr => bond
                    cur_node%child_bond => bond


                    allocate(corners)
                    nelim=2
                    atoms(1)=atom1
                    atoms(2)=atom2
                    call create_noatoms_daughter(cur_node,corners,nelim,atoms,.false.)
                    call cut_dangling_bonds(corners)
                    counter2 = counter2 + 1
                    t_pah_array(counter2)%ptr => corners
                    cur_node%child_corners => corners

                    if (ring_exists) then
                        allocate(ring)
                        nelim=6
                        call create_noatoms_daughter(cur_node,ring,nelim,sextet,.true.)
                        call cut_dangling_bonds(ring)
                        counter2 = counter2 + 1
                        t_pah_array(counter2)%ptr => ring
                        cur_node%child_ring => ring
                    end if
                else
!                    print *, 'found a disconnected pah, splitting into sons'
                    allocate(bond)
                    allocate(corners)
                    call split_structure(cur_node, bond, corners, medat)
                    if ( mod(bond%nat, 2_kint) == 1) then
!                        print *, 'found invalid splitting to odd nats'
                        call destory(bond)
                        call destory(corners)
                        deallocate(bond)
                        deallocate(corners)
                        cur_node%order = 0
                        allocate(cur_node%polynomial(1))
                        cur_node%polynomial(1)=setvli(0_kint)
                    else
!                        print *, 'found valid splitting to odd nats'
                        counter2 = counter2 + 1
                        t_pah_array(counter2)%ptr => bond
                        cur_node%child_son1 => bond
                        counter2 = counter2 + 1
                        t_pah_array(counter2)%ptr => corners
                        cur_node%child_son2 => corners
                    end if
                end if
            else
                reach_limit = .true.
                counter2 = counter2 + 1
                t_pah_array(counter2)%ptr => cur_node
            end if
        end do
        pah_array = t_pah_array
        counter1 = counter2
!        if ( mod(counter1, image_count) == 0 .or. counter1>str_per_core*image_count ) then
        if ( counter1>str_per_core*image_count ) then
            exit
        end if
    end do outer   
    final_tree_size = counter1    
    call sort_pah(pah_array, final_tree_size, .true.)
end subroutine 
end module
