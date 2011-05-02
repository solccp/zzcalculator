module database_m

    use iso_varying_string
    use structure_module

    type :: database_entry
        type(varying_string) :: key
        integer :: hits
        type(tree_node), pointer :: node => NULL()
        type(database_entry), pointer :: next => NULL()
!        type(database_entry), pointer :: previous => NULL()
    end type
    

    save
    type(database_entry), pointer :: database_head => NULL()
    integer :: max_size = 1
    logical :: set_maxsize = .false.
    integer :: database_size = 0

contains

subroutine set_max_size(nat)
    integer, intent(in) :: nat
    max_size = nat
    set_maxsize = .true.
end subroutine

subroutine add_database_entry(node, hit)
    type(tree_node), intent(inout), pointer :: node
    type(structure), pointer :: pah
    logical, intent(out) :: hit
    character(len=max_size) :: ch_key
    type(varying_string) :: str_key
    integer, dimension(max_size) :: temp_key
    type(database_entry), pointer :: new_entry
    type(database_entry), pointer :: head 

    if ( .not. set_maxsize ) then
        stop 'you forgot to call set_max_size'
    end if

    pah => node%pah

!============================================
! create key
    temp_key = 0
    do i=1, pah%nat
        temp_key(pah%indexmapping(i)) = 1
    end do

    do i = 1, size(temp_key)
        write(ch_key(i:i),'(i1)') temp_key(i)
    end do

    str_key = ch_key
!==========================================

    hit = .false.

    if ( .not. associated(database_head) ) then
!        print *, 'this should be printed only once'
        allocate(database_head)
        database_head%key = str_key
        database_head%node => node
        database_head%node%key = str_key
        database_head%hits = 1
        database_size = database_size + 1
        return
    else
        head => database_head   
        do while(associated(head%next))
            if ( str_key == head%key ) then
                hit = .true.
                head%hits = head%hits + 1
                call destory(node%pah)
                deallocate(node%pah)
                deallocate(node)
                node => head%node
                return
            end if
            if ( head%next%node%pah%nat < node%pah%nat ) then
                allocate(new_entry)
                new_entry%key = str_key
                new_entry%hits = 1
                new_entry%next => head%next 
                new_entry%node => node
                new_entry%node%key = str_key
                head%next => new_entry
                database_size = database_size + 1
                return
            end if
            head => head%next
        end do
        if ( str_key /= head%key ) then
            allocate(new_entry)
            new_entry%hits = 1
            new_entry%key = str_key
            new_entry%node => node
            new_entry%node%key = str_key
            head%next => new_entry
            database_size = database_size + 1
        else
            hit = .true.
            head%hits = head%hits + 1
            call destory(node%pah)
            deallocate(node%pah)
            deallocate(node)
            node => head%node
        end if
        return
    end if 

end subroutine

subroutine get_database_first_entry(node)
    type(tree_node), intent(inout), pointer :: node
    type(database_entry), pointer :: d_entry

    d_entry => database_head
    if ( associated(d_entry) ) then
        node => d_entry%node
        do while(associated(d_entry))
            node => d_entry%node
            if ( .not. node%hasChild ) then
                return
            end if
            d_entry => d_entry%next
        end do
        node => NULL()
    else
        node => NULL()
        return
    end if
end subroutine
subroutine remove_database_first_entry()
    if ( associated(database_head%next) ) then
        database_head => database_head%next
        database_size = database_size - 1
    else
        deallocate(database_head)
        nullify(database_head)
        database_size = 0
    end if
end subroutine

subroutine free_database
    type(database_entry), pointer :: head
    do while(associated(database_head))
        head => database_head
        database_head => database_head%next
        deallocate(head)
    end do
    nullify(database_head)
end subroutine
    

end module
