

module utils_m
    implicit none

contains
!############################### function dist ######################################
!####################################################################################
function dist(nat,atom1,atom2,geom) result(ret)
!
! compute the distance between atom1 & atom2
!
    use accuracy_m
    integer, intent(in) :: nat                           ! number of atoms
    integer, intent(in) :: atom1                         ! position of atom 1
    integer, intent(in) :: atom2                         ! position of atom 2
    real(kreal), intent(in), dimension(3,nat) :: geom    ! geometry table

    integer :: i                                         ! local counter
    real(kreal) :: r,x                                  ! local variables
    real(kreal) :: ret

    r = 0.0d0
    do i = 1, 3
        x = geom(i,atom1) - geom(i,atom2)
        r = r + x*x
    end do
    ret = sqrt(r)
    return

end function dist
!####################################################################################
!############################### end of function dist ###############################

subroutine swap_int(int1, int2)
    integer, intent(inout) :: int1, int2
    integer :: temp
    temp = int1
    int1 = int2
    int2 = temp
end subroutine

!########################## function are_neighbors ##################################
!####################################################################################
logical function are_neighbors(pah,atom1,atom2)
!
! return .true. if atoms atom1 and atom2 are neighbors in structure pah;
! otherwise, return .false.
!
    use structure_m
    implicit none

    type(structure), intent(in) :: pah
    integer, intent(in) :: atom1, atom2
    
    integer :: i
  
    are_neighbors = .false.
    do i = 1, pah%neighbornumber(atom1)
        if (pah%neighborlist(atom1,i) == atom2) then
            are_neighbors = .true.
            exit
        end if
    end do
    return

end function are_neighbors
!####################################################################################
!###################### end of function are_neighbors ###############################

subroutine print_bondlist(pah)
    use structure_m
    implicit none

    type(structure), intent(in) :: pah
    integer :: i 

    write(*,*) '----------------------------'
    if ( pah%nbondlistentries > 0 ) then
        do i = 1, pah%nbondlistentries
            write(*, '(i0,1x,a,1x,i0)') pah%bondlist(1, i), '-', pah%bondlist(2, i)
        end do
    end if
    write(*,*) '----------------------------'
end subroutine

end module utils_m
