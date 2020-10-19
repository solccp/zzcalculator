
!
!
!module for build the simple filter for structures
!
!using red-black tree as internal storage
!


module database_m

    use iso_varying_string
    use structure_m
    use rbtree_m

    implicit none

    integer :: max_size = 1
    logical :: setted_max_size = .false.

    save
    type(RedBlackTree) :: tree_total
    type(RedBlackTree) :: tree_child

    integer :: database_size = 0 

    integer :: total_str = 0
    integer :: total_hit = 0

contains

subroutine set_max_size(max_len)
    integer, intent(in) :: max_len
    max_size = max_len    
    setted_max_size = .true.
    call init(tree_total)
    call init(tree_child)
end subroutine


subroutine get_key(pah, key)
    type(structure), intent(inout) :: pah
    integer, dimension(max_size*2) :: temp_key
    character(len=max_size*2) :: ch_key
    type(varying_string), intent(out) :: key
    integer :: i

    if ( .not. setted_max_size ) then
        print *, 'error database::set_max_size was not called'
        stop
    end if

    temp_key = 0
    do i=1, pah%nat
        temp_key(pah%indexmapping(i)) = 1
    end do
    do i=1, pah%nat
        temp_key(pah%indexmapping(i)+max_size) = pah%neighbornumber(i)
    end do

    do i = 1, size(temp_key)
        write(ch_key(i:i),'(i1)') temp_key(i)
    end do

    key = ch_key
end subroutine


subroutine add_database_entry(node, hit, oldnode)
    type(pah_tree_node), intent(inout), pointer :: node
    type(pah_tree_node), intent(inout), pointer :: oldnode
    
    logical, intent(out) :: hit
    type(varying_string) :: str_key

    type(nodedata) :: data_
    logical :: inserted

    call get_key(node%pah, str_key)
    data_%ckey%key1 = node%pah%nat
    data_%ckey%key2 = str_key
    data_%value => node

    total_str = total_str + 1
    
!    if ( data_%ckey%key1 == 94 .and. data_%ckey%key2 == "111111111111111111111111111111111111111111111111111110101111111111111111111111111111111111111111" ) then
!        print *, "WTF"
!        call print_tree(tree_total%root, 1)
!        print *, '----'
!    end if

    call insert_node(tree_total, data_, inserted)

!    if ( data_%ckey%key1 == 94 .and. data_%ckey%key2 == "111111111111111111111111111111111111111111111111111110101111111111111111111111111111111111111111" ) then
!        call print_tree(tree_total%root, 1)
!    end if

    if ( inserted ) then
        hit = .false.
        call insert_node(tree_child, data_, inserted)
        database_size = tree_child%nitems
    else
        hit = .true.
        oldnode => node
        node => data_%value
        node%shared_count = node%shared_count + 1
        total_hit = total_hit + 1
    end if

end subroutine

subroutine get_database_first_entry(node)
    type(pah_tree_node), intent(inout), pointer :: node
    type(rb_treenode_ptr) :: node_ptr

    call getMax(tree_child, node_ptr)
    node => node_ptr%ptr%keyvalue%value 
    node_ptr%ptr%keyvalue%value => null()
!    call find(tree_total 
    call erase_node(tree_child, node_ptr%ptr)
    database_size = tree_child%nitems
end subroutine

subroutine free_database
    call finalize(tree_total)
    call finalize(tree_child)
end subroutine
    

end module
