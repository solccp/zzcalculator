module decompose_m

    use accuracy_m
    use structure_m
    use operator_m
    use polynomial_m
!    use output_m
    use options_m
    implicit none
    
contains
!############################ subroutine decompose ##################################
!####################################################################################
subroutine decompose(pah, level)
!
! decompose the original parent structure (pah) into three daughter structures:
!   1. with one edge bond deleted (between atom1 and atom2)
!   2. with two edge atoms (atom1 & atom2) deleted
!   3. with one ring containing atom1 & atom2 deleted
! 

    implicit none
    type(structure), intent(inout) :: pah
    integer, intent(in) :: level

    type(structure) :: bond, corners, ring
    integer :: atom1, atom2, i, nelim
    integer, dimension(6) :: sextet
    integer, dimension(2) :: atoms
    logical :: ring_exists, are_neighbors
    

    ring_exists = .false.
! #############################################
! # select an edge bond between atom1 & atom2 #
! #############################################
    call select_edge_bond(pah, atom1, atom2)


! ##################################
! # create the daughter structures #
! ##################################
    call create_nobond_daughter(pah, bond, atom1, atom2)
    nelim = 2
    atoms(1) = atom1
    atoms(2) = atom2
    call create_noatoms_daughter(pah, corners, nelim, atoms, .false. )


! ###################################################
! # eliminate dangling bonds in daughter structures #
! ###################################################
    call cut_dangling_bonds(bond)
    call cut_dangling_bonds(corners)


    
! ###############################################
! # find ZZ polynomials for daughter structures #
! ###############################################
    call find_ZZ_polynomial(bond, level+1)
    call find_ZZ_polynomial(corners, level+1)

! ##########################################
! # find the ring containing atom1 & atom2 #
! ##########################################
    if ( .not. options%kekule_only ) then
        call find_edge_ring(pah, sextet, atom1, atom2, ring_exists)
        if (ring_exists) then
            nelim = 6
            call create_noatoms_daughter(pah, ring, nelim, sextet, .true. )
            call cut_dangling_bonds(ring)
            call find_ZZ_polynomial(ring, level+1)
        end if
    end if

! ###############################################
! # find ZZ polynomial for the parent structure #
! ###############################################
    call sum_polynomials(pah, bond, corners, ring, ring_exists)

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



!###################### subroutine split_and_decompose ##############################
!####################################################################################
subroutine split_and_decompose(pah, medat, level)
!
! split a disconnected polycyclic benzenoid structure pah into two substructures
! * son1 which is connected and contains (medat-1) atoms
! * son2 which can be connected or disconnected and contains (pah%nat-medat+1) atoms
! decompose both structures further and multiply their resulting ZZ polynomials
!
!        ZZ(pah) = ZZ(son1) * ZZ(son2)
!
    use output_m
    use options_m
    implicit none
    
    type(structure), intent(inout) :: pah
    integer, intent(in) :: medat, level

    integer :: i,j
    type(structure) :: son1, son2


    call split_structure(pah, son1, son2, medat)

    if (options%print_intermediate_structures) then    

    !##########################################
    !# allocate storage for connection output #
    !##########################################
        allocate(son1%doublebondlist(2,size(pah%doublebondlist,2)))
        allocate(son1%ringlist(6,size(pah%ringlist,2)))
        allocate(son2%doublebondlist(2,size(pah%doublebondlist,2)))
        allocate(son2%ringlist(6,size(pah%ringlist,2)))

    !##########################################################
    ! # fragments only contain the rest connection infomation #
    !##########################################################
        son1%doublebondnumber = 0
        son1%ringnumber = 0
        son2%doublebondnumber = 0
        son2%ringnumber = 0


    !######################
    !# open scratch file for sons
    !######################
        i = getunit()
        open(unit=i, status='scratch')
        son1%storage_unit = i
        i = getunit()
        open(unit=i, status='scratch')
        son2%storage_unit = i
    
        son1%hasDisconnectedParent = .true.
        son2%hasDisconnectedParent = .true.
    end if

! ###################################################
! # find the ZZ polynomials for both son structures #
! ###################################################
    call find_ZZ_polynomial(son1,level+1)
    call find_ZZ_polynomial(son2,level+1)

! ######################################################
! # multiply the ZZ polynomials of both son structures #
! ######################################################
    call multiply_polynomials(pah,son1,son2)

    if (options%print_intermediate_structures) then
    !########################################################
    !# combine all connection info from pah, son1, and son2 #
    !########################################################
        call combine_connection_output(pah,son1,son2)


        close(son1%storage_unit,status='delete')
        close(son2%storage_unit,status='delete')
    end if

! #################################
! # deallocate the son structures #
! #################################
    call destory(son1)
    call destory(son2)
    return

end subroutine split_and_decompose
!####################################################################################
!################### end of subroutine split_and_decompose ##########################

end module
