
module partition_m

contains
subroutine bipartition(pah)
    use accuracy_m
    use structure_m
    use lapack

    implicit none
    type(structure), intent(inout) :: pah


    logical, dimension(pah%nat) :: g1, g2

    double precision , dimension(pah%nat, pah%nat) :: LP_mat
    double precision , dimension(pah%nat) :: ev
    integer, dimension(2,pah%nat) :: bondlist

    integer :: i, j, k

    LP_mat = 0.0d0

    do i = 1, pah%nat
        LP_mat(i,i) = pah%neighbornumber(i)
        do k = 1, pah%neighbornumber(i)
           LP_mat(i, pah%neighborlist(i, k)) = -1
        end do
    end do

    call syev('V', 'L', pah%nat, LP_mat, pah%nat, ev)
    g1 = .false.
    g2 = .false.
    do i=1, pah%nat
        if (LP_mat(i,2)> 0) then
            g1(i) = .true.
        else
            g2(i) = .true.
        end if
    end do

    k = 0
    do i = 1, pah%nat
        do j = 1, pah%neighbornumber(i)
            if ( g1(i) == .true. .and. g2(pah%neighborlist(i, j)) == .true. ) then
                k = k + 1
                bondlist(1, k) = i
                bondlist(2, k) = pah%neighborlist(i, j)
            end if
        end do      
    end do

!    do i = 1, k 
!        print *, pah%indexmapping(bondlist(1,i)), pah%indexmapping(bondlist(2,i))
!    end do
    pah%nbondlistentries = k
    allocate(pah%bondlist(2, k))
    pah%bondlist = bondlist(:, :k)

end subroutine

end module
