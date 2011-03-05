!############################ subroutine decompose ##################################
!####################################################################################
subroutine decompose(pah,level)
!
! decompose the original parent structure (pah) into three daughter structures:
!   1. with one edge bond deleted (between atom1 and atom2)
!   2. with two edge atoms (atom1 & atom2) deleted
!   3. with one ring containing atom1 & atom2 deleted
! 
    use types_module
    use structure_module

    implicit none
    type(structure), intent(in) :: pah
    integer(kint), intent(in) :: level

    type(structure) :: bond, corners, ring
    integer(kint) :: atom1,atom2,i,nelim
    integer(kint), dimension(6) :: sextet
    integer(kint), dimension(2) :: atoms
    logical :: ring_exists, are_neighbors

! #############################################
! # select an edge bond between atom1 & atom2 #
! #############################################
    call select_edge_bond(pah,atom1,atom2)

! ##########################################
! # find the ring containing atom1 & atom2 #
! ##########################################
    call find_edge_ring(pah,sextet,atom1,atom2,ring_exists)

! ##################################
! # create the daughter structures #
! ##################################
    call create_nobond_daughter(pah,bond,atom1,atom2)
    nelim=2
    atoms(1)=atom1
    atoms(2)=atom2
    call create_noatoms_daughter(pah,corners,nelim,atoms,.false.)
    if (ring_exists) then
        nelim=6
        call create_noatoms_daughter(pah,ring,nelim,sextet,.true.)
    end if

! ###################################################
! # eliminate dangling bonds in daughter structures #
! ###################################################
    call cut_dangling_bonds(bond)
    call cut_dangling_bonds(corners)
    if (ring_exists) then
        call cut_dangling_bonds(ring)
    end if

! ###############################################
! # find ZZ polynomials for daughter structures #
! ###############################################
    call find_ZZ_polynomial(bond,level+1)
    call find_ZZ_polynomial(corners,level+1)
    if (ring_exists) then
        call find_ZZ_polynomial(ring,level+1)
    end if

! ###############################################
! # find ZZ polynomial for the parent structure #
! ###############################################
    call sum_polynomials(pah,bond,corners,ring,ring_exists)

! ##################################
! # deallocate daughter structures #
! ##################################
    deallocate(bond%neighbornumber)
    deallocate(corners%neighbornumber)
    deallocate(bond%neighborlist)
    deallocate(corners%neighborlist)
    deallocate(bond%polynomial)
    deallocate(corners%polynomial)
    if (ring_exists) then
        deallocate(ring%neighbornumber)
        deallocate(ring%neighborlist)
        deallocate(ring%polynomial)
        if (ring%nbondlistentries > 0) then
            deallocate(ring%bondlist)
        end if
    end if
    if (bond%nbondlistentries > 0) then
        deallocate(bond%bondlist)
    end if
    if (corners%nbondlistentries > 0) then
        deallocate(corners%bondlist)
    end if

end subroutine decompose
!####################################################################################
!######################## end of subroutine decompose ###############################
