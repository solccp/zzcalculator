!######################## subroutine cut_dangling_bonds #############################
subroutine remove_atom(pah, atom, r1)
    use types_module
    use structure_module
    implicit none
    type(structure), intent(inout) :: pah
    integer(kint), intent(in) :: atom
    integer(kint), intent(in) :: r1

    integer(kint) :: i, j, k
    integer(kint) :: n1
    do j=1, pah%neighbornumber(atom)
        n1 = pah%neighborlist(atom, j)
        do k = 1, pah%neighbornumber(n1)
            if (pah%neighborlist(n1,k) == atom ) then
                pah%neighborlist(n1,k:3-1) = pah%neighborlist(n1,k+1:3)
                pah%neighbornumber(n1) = pah%neighbornumber(n1) - 1
                exit
            end if
        end do
    end do
    if ( r1 /= atom ) then
        pah%neighbornumber(atom) = pah%neighbornumber(r1)
        pah%neighborlist(atom,:) = pah%neighborlist(r1,:)
        do j=1, pah%neighbornumber(r1)
            n1 = pah%neighborlist(r1, j)
            do k = 1, pah%neighbornumber(n1)
                if (pah%neighborlist(n1,k) == r1 ) then
                    pah%neighborlist(n1,k) = atom
                    exit
                end if
            end do
        end do
    end if
    
    pah%neighbornumber(r1) = 0
    pah%nat = pah%nat - 1


end subroutine

subroutine get_remove_indexes(pah, atom1, atom2, r1, r2)
    use types_module
    use structure_module
    implicit none
    type(structure), intent(in) :: pah
    integer(kint), intent(out) :: r1, r2
    integer(kint), intent(in) :: atom1, atom2
    integer(kint) :: temp

    r1 = pah%nat 
    if (r1 == atom2) then
        r2 = r1
        r1 = r1 - 1
    else
        r2 = r1 - 1
    end if

    if (r2 == atom1) then
        temp = r2
        r2 = r1
        r1 = temp
    end if

        
end subroutine

!####################################################################################
subroutine cut_dangling_bonds(pah)
!
! removes all dangling bonds from a given polycyclic benzenoid structure (pah)
! the dangling bond is defined via a single connected atom
! a structure without singly connected atoms does not have dangling bonds
!
    use types_module
    use structure_module
    use options_m
    implicit none
    type(structure), intent(inout) :: pah

    integer(kint) :: atom1, atom2, r1, r2
    type(structure) :: pah1
    logical :: has_dangling_bonds
    integer(kint) :: i, j


! ############################
! # eliminate dangling bonds #
! ############################
    has_dangling_bonds = .true.
    do while (has_dangling_bonds)
        call check_for_dangling_bonds(pah,has_dangling_bonds,atom1)
        if (has_dangling_bonds) then
            atom2 = pah%neighborlist(atom1,1)
            call get_remove_indexes(pah, atom1, atom2, r1, r2)
            if ( options%print_intermediate_structures) then
                pah%doublebondnumber = pah%doublebondnumber + 1
                pah%doublebondlist(1, pah%doublebondnumber) = pah%indexmapping(atom1)
                pah%doublebondlist(2, pah%doublebondnumber) = pah%indexmapping(atom2)
                pah%indexmapping(atom1) = pah%indexmapping(r1)
                pah%indexmapping(atom2) = pah%indexmapping(r2)
            endif

            call remove_atom(pah, atom1, r1)
            call remove_atom(pah, atom2, r2)
            ! #######################################
            ! # translate to the new atom numbering #
            ! #######################################
            do i = 1, pah%nbondlistentries
                do j = 1, 2
                    if ( pah%bondlist(j,i) == r1 ) pah%bondlist(j,i) = atom1
                    if ( pah%bondlist(j,i) == r2 ) pah%bondlist(j,i) = atom2
                    if ( pah%bondlist(j,i) == atom1 ) pah%bondlist(j,i) = 0
                    if ( pah%bondlist(j,i) == atom2 ) pah%bondlist(j,i) = 0
                end do
            end do
        end if
    end do
  
    call clean_bond_list(pah)

end subroutine cut_dangling_bonds
!####################################################################################
!##################### end of subroutine cut_dangling_bonds #########################



!####################### subroutine check_for_dangling_bonds ########################
!####################################################################################
subroutine check_for_dangling_bonds(pah,has_dangling_bonds,atom1)
!
! check if a given polycyclic benzenoid structure
! has dangling bonds
!
    use types_module
    use structure_module
    implicit none
    type(structure), intent(in) :: pah
    logical, intent(out) :: has_dangling_bonds
    integer(kint), intent(out) :: atom1

    integer(kint) :: i

    has_dangling_bonds = .false.
    do i=1, pah%nat
        if (pah%neighbornumber(i) == 1) then
            atom1 = i
            has_dangling_bonds = .true.
            exit
        end if
    end do
    return

end subroutine check_for_dangling_bonds
!####################################################################################
!#################### end of subroutine check_for_dangling_bonds ####################


!######################## subroutine find_aromatic_sextet ###########################
!####################################################################################
subroutine find_aromatic_sextet(pah,sextet,atom1,atom2,atom3,ring_exists)
!
! for a given polycyclic benzenoid structure, the routine
! finds an aromatic ring containing atoms: atom1, atom2, and atom3
!
    use types_module
    use structure_module
    implicit none
    type(structure), intent(in) :: pah
    integer(kint), intent(out), dimension(6) :: sextet
    logical, intent(out) :: ring_exists
    integer(kint), intent(in) :: atom1,atom2,atom3

    integer(kint) :: i,j,k,l
    integer(kint) :: atom4,atom5,atom6

! ##################
! # initialization #
! ##################
    ring_exists=.false.
    sextet(1)=atom1
    sextet(2)=atom2
    sextet(3)=atom3

! ######################################
! # look for the remaining three atoms #
! ######################################
    do i=1,pah%neighbornumber(atom3)
        atom4=pah%neighborlist(atom3,i)
        if (atom4 /= atom1) then
            do j=1,pah%neighbornumber(atom2)
                atom5=pah%neighborlist(atom2,j)
                if (atom5 /= atom1) then
                    do k=1,pah%neighbornumber(atom4)
                        atom6=pah%neighborlist(atom4,k)
                        do l=1,pah%neighbornumber(atom5)
                            if (atom6 == pah%neighborlist(atom5,l)) then
                                sextet(4) = atom4
                                sextet(5) = atom5
                                sextet(6) = atom6
                                ring_exists = .true.
                            end if
                        end do
                    end do
                end if
            end do
        end if
    end do
    return

end subroutine find_aromatic_sextet
!####################################################################################
!##################### end of subroutine find_aromatic_sextet #######################




!########################### subroutine find_edge_ring ##############################
!####################################################################################
subroutine find_edge_ring(pah,sextet,atom1,atom2,ring_exists)
!
! for a given polycyclic benzenoid structure, the routine
! finds an aromatic ring containing atoms: atom1 and atom2
!
    use types_module
    use structure_module
    implicit none
    type(structure), intent(in) :: pah
    integer(kint), intent(out), dimension(6) :: sextet
    integer(kint), intent(in) :: atom1,atom2
    logical, intent(out) :: ring_exists
    
    integer(kint) :: i,j,k,l,m
    integer(kint) :: atom3,atom4,atom5,atom6

! ##################
! # initialization #
! ##################
    ring_exists = .false.
    sextet(1) = atom1
    sextet(2) = atom2

! #####################################
! # look for the remaining four atoms #
! #####################################
    outer: do i=1,pah%neighbornumber(atom1)
        atom3=pah%neighborlist(atom1,i)
        if (atom3 /= atom2) then
            do j=1,pah%neighbornumber(atom2)
                atom4=pah%neighborlist(atom2,j)
                if (atom4 /= atom1) then
                    do k=1,pah%neighbornumber(atom3)
                        atom5=pah%neighborlist(atom3,k)
                        if (atom5 /= atom1) then
                            do l=1,pah%neighbornumber(atom4)
                                atom6=pah%neighborlist(atom4,l)
                                if (atom6 /= atom2) then
                                    do m=1,pah%neighbornumber(atom5)
                                        if (atom6 == pah%neighborlist(atom5,m)) then
                                            sextet(3) = atom4
                                            sextet(4) = atom6
                                            sextet(5) = atom5
                                            sextet(6) = atom3
                                            ring_exists = .true.
                                            exit outer
                                        end if
                                    end do
                                end if
                            end do
                        end if
                    end do
                end if
            end do
        end if
    end do outer
    return

end subroutine find_edge_ring
!####################################################################################
!######################## end of subroutine find_edge_ring ##########################

subroutine swap_int(int1, int2)
    use types_module
    integer(kint), intent(inout) :: int1, int2
    integer(kint) :: temp
    temp = int1
    int1 = int2
    int2 = temp
end subroutine

!######################## subroutine check_if_connected #############################
!####################################################################################
subroutine check_if_connected(pah,medat)
!
! check if the given polycyclic benzenoid structure is connected
! i.e., if there exists a bonded path between any two randomly chosen atoms
! * if pah is connected, medat is returned as 0
! * if pah is disconnected, the structure is reordered in such a way that
!   first connected substructure, containing (medat-1) atoms, is returned
!   in first (medat-1) positions of pah, and the remaining substructure (possibly
!   disconnected) in the remaining positions
!
    use types_module
    use structure_module
    use options_m
    use temp_space
    implicit none
    type(structure), intent(inout) :: pah
    integer(kint), intent(out) :: medat

    type(structure) :: pah1
    integer(kint) :: i,j,k,lnat,start, l, n1

!    map => int_1darray_1
!    visit_list => bool_1darray_1

! #########################################################
! # find the connected cluster of atoms containing atom 1 #
! #########################################################
    lnat = 0
    start = 1
    bool_1darray_1(1:pah%nat) = .false.
    call dfs(pah, pah%nat, bool_1darray_1, lnat)

! ##########################################
! # return if all atoms are in the cluster #
! ##########################################
    if (lnat == pah%nat) then
        medat = 0
        return
    else
        medat = lnat+1
    end if

    i = 1
    j = pah%nat
    do while( i < j )
        do while( bool_1darray_1(i) )
            int_1darray_1(i) = i
            i = i + 1
        end do
        do while( .not. bool_1darray_1(j) )
            int_1darray_1(j) = j
            j = j - 1
        end do
        if ( i >= j ) then
            exit
        end if
        if (options%print_intermediate_structures) then
            call swap_int(pah%indexmapping(i), pah%indexmapping(j))
        end if
        call swap_int(pah%neighbornumber(i), pah%neighbornumber(j))
        do k = 1, 3
            call swap_int(pah%neighborlist(i,k), pah%neighborlist(j,k))
        end do
        do k = 1, pah%neighbornumber(i)
            n1 = pah%neighborlist(i, k)
            do l = 1, pah%neighbornumber(n1) 
                if ( pah%neighborlist(n1,l) == j ) then
                    pah%neighborlist(n1,l) = i
                    exit
                end if
            end do
        end do
        do k = 1, pah%neighbornumber(j)
            n1 = pah%neighborlist(j, k)
            do l = 1, pah%neighbornumber(n1) 
                if ( pah%neighborlist(n1,l) == i ) then
                    pah%neighborlist(n1,l) = j
                    exit
                end if
            end do
        end do
        int_1darray_1(i) = j
        int_1darray_1(j) = i
        i = i + 1
        j = j - 1
        
    end do

! ######################################
! # map the bond list to the new order #
! ######################################
    if (pah%nbondlistentries > 0) then
        forall (i=1:pah%nbondlistentries, j=1:2) !, pah%bondlist(j,i) /= 0)
            pah%bondlist(j,i)=int_1darray_1(pah%bondlist(j,i))
        end forall
    end if
    return

end subroutine check_if_connected
!####################################################################################
!###################### end of subroutine check_if_connected ########################



!################################# subroutine dfs ###################################
!####################################################################################
subroutine dfs(pah,nat,visit_list,lnat)
!
! find all atoms in structure pah that are connected to atom 1; if atom k is connected 
! to atom 1 via a sequence of bonds, then visit_list(k)=.true.; otherwise, visit_list(k)=.false.
!
    use types_module
    use structure_module
    implicit none
    type(structure), intent(in) :: pah
    integer(kint), intent(in) :: nat
    logical, dimension(nat), intent(inout) :: visit_list
    integer(kint), intent(inout) :: lnat
    
    integer(kint) :: i
    integer(kint), dimension(0:nat) :: stack
    integer(kint) :: cur_index

    stack(0) = 1
    stack(1) = 1
    do while(stack(0) > 0)
        cur_index = stack(stack(0))
        stack(0) = stack(0) - 1
        if (.not. visit_list(cur_index)) then
            visit_list(cur_index) = .true.
            lnat = lnat + 1
            do i = 1, pah%neighbornumber(cur_index)
                if (.not. visit_list(pah%neighborlist(cur_index, i))) then
                    stack(0) = stack(0) + 1
                    stack(stack(0)) = pah%neighborlist(cur_index, i)
                end if
            end do
        end if
    end do
    return

end subroutine dfs
!####################################################################################
!############################# end of subroutine dfs ################################


!###################### subroutine split_and_decompose ##############################
!####################################################################################
subroutine split_and_decompose(pah,medat,level)
!
! split a disconnected polycyclic benzenoid structure pah into two substructures
! * son1 which is connected and contains (medat-1) atoms
! * son2 which can be connected or disconnected and contains (pah%nat-medat+1) atoms
! decompose both structures further and multiply their resulting ZZ polynomials
!
!        ZZ(pah) = ZZ(son1) * ZZ(son2)
!
    use types_module
    use output
    use structure_module
    use options_m
    implicit none
    
    type(structure), intent(in) :: pah
    integer(kint), intent(in) :: medat, level

    integer(kint) :: i,j
    type(structure) :: son1, son2

! ###############################
! # allocate the son structures #
! ###############################
    son1%nat=medat-1
    son1%nbondlistentries=pah%nbondlistentries
    allocate(son1%neighbornumber(son1%nat))
    allocate(son1%neighborlist(son1%nat,3))
    son2%nat=pah%nat-medat+1
    son2%nbondlistentries=pah%nbondlistentries
    allocate(son2%neighbornumber(son2%nat))
    allocate(son2%neighborlist(son2%nat,3))
    if (pah%nbondlistentries > 0) then 
        allocate(son1%bondlist(2,son1%nbondlistentries))
        allocate(son2%bondlist(2,son2%nbondlistentries))
        son1%bondlist=pah%bondlist
        son2%bondlist=pah%bondlist
    end if

    if (options%print_intermediate_structures) then
    ! ########################
    ! # update index mapping #
    ! ########################
        allocate(son1%indexmapping(son1%nat))
        son1%indexmapping = pah%indexmapping(1:medat-1)
        allocate(son2%indexmapping(son2%nat))
        son2%indexmapping = pah%indexmapping(medat:)

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

! #################################
! # initialize the son structures #
! #################################
    son1%neighbornumber=pah%neighbornumber(1:medat-1)
    son1%neighborlist=pah%neighborlist(1:medat-1,1:3)
    son2%neighbornumber=pah%neighbornumber(medat:pah%nat)
    son2%neighborlist=pah%neighborlist(medat:pah%nat,1:3)
    forall (i=1:son2%nat, j=1:3, son2%neighborlist(i,j) /= 0)
        son2%neighborlist(i,j)=son2%neighborlist(i,j)-son1%nat
    end forall
    forall (i=1:son1%nbondlistentries, j=1:2, son1%bondlist(j,i) > son1%nat)
        son1%bondlist(j,i)=0
    end forall
    forall (i=1:son2%nbondlistentries, j=1:2, son2%bondlist(j,i) <= son1%nat)
        son2%bondlist(j,i)=0
    end forall
    forall (i=1:son2%nbondlistentries, j=1:2, son2%bondlist(j,i) /= 0)
        son2%bondlist(j,i)=son2%bondlist(j,i)-son1%nat
    end forall

    call clean_bond_list(son1)
    call clean_bond_list(son2)


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







!######################### subroutine select_edge_bond ##############################
!####################################################################################
subroutine select_edge_bond(pah,atom1,atom2)
!
! selects two atoms (atom1 & atom2) located on the edge of a given
! polycyclic benzenoid structure pah; it uses the fact that the edge
! carbon atoms are surounded maximally by 2 hexagons
!
    use types_module
    use structure_module
    implicit none
    type(structure), intent(in) :: pah
    integer(kint), intent(out) :: atom1, atom2
    
    integer(kint) :: i,j,atom3,sextet(6)
    logical :: ring_exists,selected
    
    atom2 = 0
    selected = .false.

! #################################################
! # try to select the bond from the provided list #
! #################################################
    outer1: do i=1,pah%nbondlistentries
        atom1=pah%bondlist(1,i)
        atom2=pah%bondlist(2,i)
        if (pah%neighbornumber(atom1) == 2 .or. pah%neighbornumber(atom2) == 2) then
            selected=.true.
            exit outer1
        end if
        do  j=1,pah%neighbornumber(atom1)
            if (pah%neighborlist(atom1,j) == atom2) cycle
            atom3=pah%neighborlist(atom1,j)
            call find_aromatic_sextet(pah,sextet,atom1,atom2,atom3,ring_exists)
            if (.not. ring_exists) then
                selected=.true.
                exit outer1
            end if
        end do
    end do outer1

! #############################################################
! # otherwise choose an edge bond from the topological matrix #
! #############################################################
    if (.not. selected) then
        atom2 = 0
        outer: do i=1, pah%nat
            if (pah%neighbornumber(i) == 2) then
                atom1 = i
                atom2 = pah%neighborlist(i,1)
                exit outer
            else if (pah%neighbornumber(i) == 3) then
                atom1 = i
                do j=1, 3
                    atom2=pah%neighborlist(atom1,min(j,mod(j,3_kint)+1))
                    atom3=pah%neighborlist(atom1,max(j,mod(j,3_kint)+1))
                    call find_aromatic_sextet(pah,sextet,atom1,atom2,atom3,ring_exists)
                    if (.not. ring_exists) then
                        exit outer
                    end if
                end do
                atom2=0
            end if
        end do outer
    end if

    if (atom2 == 0) then
        write(*,*)"No atom is located on the edge"
        write(*,*)"despite of non-zero number of atoms:",pah%nat
        write(*,*)"Logical error - program is terminated in select_edge_bond"
        stop
    end if
    return

end subroutine select_edge_bond
!####################################################################################
!###################### end of subroutine select_edge_bond ##########################




!########################## subroutine clean_bond_list ##############################
!####################################################################################
subroutine clean_bond_list(pah)
!
! clean the list of bonds by removing not longer
! valid entries, i.e. not existing atoms & not existing bonds
!
    use types_module
    use structure_module
    implicit none
    type(structure), intent(inout) :: pah
    
    integer(kint) :: i,j
    type(structure) :: pah1
    logical :: are_neighbors

    if (pah%nbondlistentries > 0) then
        allocate(pah1%bondlist(2,pah%nbondlistentries))
        pah1%bondlist=0
        j=0
        do i=1,pah%nbondlistentries
            if (pah%bondlist(1,i) /= 0 .and. pah%bondlist(1,i)<=pah%nat) then
                if (pah%bondlist(2,i) /= 0 .and. pah%bondlist(2,i)<=pah%nat) then
                    if (are_neighbors(pah,pah%bondlist(1,i),pah%bondlist(2,i))) then
                        j=j+1
                        pah1%bondlist(1,j)=pah%bondlist(1,i)
                        pah1%bondlist(2,j)=pah%bondlist(2,i)
                    end if
                end if
            end if
        end do
        pah%nbondlistentries=j

        if (pah%nbondlistentries == 0) then
            deallocate(pah%bondlist)
        else
            deallocate(pah%bondlist)
            allocate(pah%bondlist(2,pah%nbondlistentries))
            pah%bondlist=pah1%bondlist(:,1:pah%nbondlistentries)
            deallocate(pah1%bondlist)
        end if
    end if
    return

! #################################################
! # try to select the bond from the provided list #
! #################################################


end subroutine clean_bond_list
!####################################################################################
!####################### end of subroutine clean_bond_list ##########################


!########################## function are_neighbors ##################################
!####################################################################################
logical function are_neighbors(pah,atom1,atom2)
!
! return .true. if atoms atom1 and atom2 are neighbors in structure pah;
! otherwise, return .false.
!
    use types_module
    use structure_module
    implicit none
    type(structure), intent(in) :: pah
    integer(kint), intent(in) :: atom1,atom2
    
    integer(kint) :: i
  
    are_neighbors=.false.
    do i=1,pah%neighbornumber(atom1)
        if (pah%neighborlist(atom1,i) == atom2) then
            are_neighbors=.true.
            exit
        end if
    end do
    return

end function are_neighbors
!####################################################################################
!###################### end of function are_neighbors ###############################
