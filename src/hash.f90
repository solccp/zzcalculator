module hash_m
    use types_module
    use structure_module
    save 
    integer(kint) :: stat_hit = 0
    integer(kint) :: stat_no_hit = 0
    
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
!    print *, 'hash: ', key

    deallocate(edges)
end subroutine

subroutine get_polynomial(pah, hit)
    type(structure), intent(inout) :: pah
    logical, intent(inout) :: hit
    integer :: order

    hit = .false.

    call get_polynomial_order_kernel(pah%hash_key, order)
    if ( order == -1) then
        stat_no_hit = stat_no_hit + 1
        return 
    end if
    
    hit = .true.
    stat_hit = stat_hit + 1

    if ( allocated(pah%polynomial) .and. pah%order /= order ) then
        deallocate(pah%polynomial)
        pah%order = order
        allocate(pah%polynomial(pah%order+1))
    else if (.not. allocated(pah%polynomial) ) then
        pah%order = order
        allocate(pah%polynomial(pah%order+1))
    end if

    call get_polynomial_kernel(pah%hash_key, block_size, pah%polynomial)

end subroutine

subroutine add_polynomial(pah)
    use types_module
    type(structure), intent(in) :: pah
    integer :: i,j
    
!    call print_ZZ_polynomial(pah)
    
!    print *, pah%order+1
!    do i = 1, pah%order+1
!        print *, 'term ', i
!        print *, '  leadpow', pah%polynomial(i)%leadpow
!        do j=1, block_size
!           print *, '  coeffs[', j ,']:', pah%polynomial(i)%tabl(j) 
!        end do
!    end do

!    print *, '=============add========='
    call add_polynomial_kernel(pah%hash_key, pah%order+1, block_size, pah%polynomial)

!    print *, 'add done'
    
end subroutine

end module
