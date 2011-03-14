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
    use mpi_global

    implicit none
    type(structure), intent(in) :: pah
    integer(kint), intent(in) :: level

    type(structure) :: bond, corners, ring
    integer(kint) :: atom1,atom2,i,nelim
    integer(kint), dimension(6) :: sextet
    integer(kint), dimension(2) :: atoms
    logical :: ring_exists, are_neighbors
    integer :: free_node, ring_node, bond_node

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
    if (ring_exists) then
        if (image_id == 0) then
            ring_node = get_free_node()    
            if ( ring_node /= -1 ) then
                call send_structure(ring, level+1, ring_node)
            else
                call find_ZZ_polynomial(ring,level+1)
            end if
        else
            call find_ZZ_polynomial(ring,level+1)
        end if
    end if
    
    if (image_id == 0) then
        bond_node = get_free_node()    
        if ( bond_node /= -1 ) then
            call send_structure(bond, level+1, bond_node)
        else
            call find_ZZ_polynomial(bond,level+1)
        end if
    else
        call find_ZZ_polynomial(bond,level+1)
    end if

    call find_ZZ_polynomial(corners,level+1)
    if ( image_id == 0 ) then
        if ( ring_exists ) then
            if ( ring_node /= -1 ) then
                call recv_polynomial(ring, ring_node) 
            end if
        end if
        if ( bond_node /= -1 ) then
            call recv_polynomial(bond, bond_node) 
        end if
    end if


! ###############################################
! # find ZZ polynomial for the parent structure #
! ###############################################
    call sum_polynomials(pah,bond,corners,ring,ring_exists)

! ##################################
! # deallocate daughter structures #
! ##################################
    call destory(bond)
    call destory(corners)
    if (ring_exists) then
        call destory(ring)
    end if

end subroutine decompose
!####################################################################################
!######################## end of subroutine decompose ###############################
