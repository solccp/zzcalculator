module hash_m
    use types_module
    use structure_module
contains

subroutine get_hash(pah, key)
    type(structure), intent(in) :: pah
    character(len=32), intent(out) :: key
    
    integer :: i, j
    integer :: numEdges
    integer, dimension(:,:), allocatable :: edges

    numEdges = 0
    do i = 1, pah%nat
        do j = 1, pah%neighbornumber(i)
            if ( pah%neighborlist(i, j) > i ) then
                numEdges = numEdges + 1
            end if
        end do
    end do

    allocate(edges(2,numEdges))
    numEdges = 0
    do i = 1, pah%nat
        do j = 1, pah%neighbornumber(i)
            if ( pah%neighborlist(i, j) > i ) then
                numEdges = numEdges + 1
                edges(1, numEdges) = i
                edges(2, numEdges) = pah%neighborlist(i, j)
            end if
        end do
    end do

    call make_hash(pah%nat, numEdges, edges, key)
    print *, 'hash: ', key

    deallocate(edges)
end subroutine
end module
