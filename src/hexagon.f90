!####################### subroutine find_all_hexagons ###############################
!####################################################################################
subroutine find_all_hexagons(nat,pah,nhex,lista)
! 
! find a list of all hexagons in a polycyclic structure pah
!
    use types_module
    use structure_module
    implicit none
    integer(kint), intent(in) :: nat
    type(structure), intent(in) :: pah
    integer(kint), intent(out) :: nhex ! number of hexagons in structure pah
    integer(kint), intent(inout), dimension(6,nat) :: lista
    
    integer(kint) :: i,j,k,l,atom2,atom3
    integer(kint),dimension(6) :: sextet
    logical :: ring_exists

    nhex=0
! #######################
! # loop over all atoms #
! #######################
    atomloop: do i=1,pah%nat

!   ##################################
!   # if atom i has only 2 neighbors #
!   ##################################
        if (pah%neighbornumber(i) == 2) then
            if (pah%neighborlist(i,1) < i) cycle atomloop
            if (pah%neighborlist(i,2) < i) cycle atomloop
            call find_aromatic_sextet(pah,sextet,i,pah%neighborlist(i,2),pah%neighborlist(i,1),ring_exists)
            if (ring_exists) then
                do j=4,6
                    if (sextet(j) < i) cycle atomloop
                end do
                nhex = nhex+1
                do l=1,6
                    lista(l,nhex)=sextet(l)
                end do
            end if

!   #############################
!   # if atom i has 3 neighbors #
!   #############################
        else if (pah%neighbornumber(i) == 3) then

!     ###################################################
!     # loop over all 2-combinations of three neighbors #
!     ###################################################
            innerloop: do j=1,3
                atom2=pah%neighborlist(i,mod(j,3)+1)
                if (atom2 < i) cycle innerloop
                atom3=pah%neighborlist(i,mod(j+1,3)+1)
                if (atom3 < i) cycle innerloop
                call find_aromatic_sextet(pah,sextet,i,atom2,atom3,ring_exists)
                if (ring_exists) then
                    do k=4,6
                        if (sextet(k) < i) cycle innerloop
                    end do
                    nhex = nhex+1
                    do l=1,6
                        lista(l,nhex)=sextet(l)
                    end do
                end if
            end do innerloop
        end if
    end do atomloop
    return

end subroutine find_all_hexagons
!####################################################################################
!#################### end of subroutine find_all_hexagons ###########################
