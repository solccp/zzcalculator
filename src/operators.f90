module operator_m

    public :: check_if_connected, select_edge_bond, find_edge_ring, create_nobond_daughter, create_noatoms_daughter
    public :: cut_dangling_bonds, split_structure, clean_bond_list
    private

contains
!########################## subroutine clean_bond_list ##############################
!####################################################################################
subroutine clean_bond_list(pah)
!
! clean the list of bonds by removing not longer
! valid entries, i.e. not existing atoms & not existing bonds
!
    use structure_m
    use utils_m
    implicit none
    type(structure), intent(inout) :: pah
    
    integer :: i,j
    integer, allocatable, dimension(:,:) :: bl
!    logical :: are_neighbors

    if (pah%nbondlistentries > 0) then
!        print *, 'before clean bondlist'
!        call print_bondlist(pah)
        allocate(bl(2,pah%nbondlistentries))
        bl = 0
        j = 0
        do i = 1, pah%nbondlistentries
            if (pah%bondlist(1,i) /= 0 .and. pah%bondlist(1,i)<=pah%nat) then
                if (pah%bondlist(2,i) /= 0 .and. pah%bondlist(2,i)<=pah%nat) then
                    if (are_neighbors(pah,pah%bondlist(1,i),pah%bondlist(2,i))) then
                        j = j + 1
                        bl(1,j) = pah%bondlist(1,i)
                        bl(2,j) = pah%bondlist(2,i)
                    end if
                end if
            end if
        end do
        pah%nbondlistentries = j

        if (pah%nbondlistentries == 0) then
            deallocate(pah%bondlist)
            deallocate(bl)
        else
            deallocate(pah%bondlist)
            allocate(pah%bondlist(2,pah%nbondlistentries))
            pah%bondlist = bl(:,1:pah%nbondlistentries)
            deallocate(bl)
        end if
    end if

!    if (pah%nbondlistentries > 0) then
!        print *, 'after clean bondlist'
!        call print_bondlist(pah)
!    end if

    return

! #################################################
! # try to select the bond from the provided list #
! #################################################


end subroutine clean_bond_list
!####################################################################################
!####################### end of subroutine clean_bond_list ##########################

!##################### subroutine create_nobond_daughter ############################
!####################################################################################
subroutine create_nobond_daughter(pah,bond,atom1,atom2)
!
! creates a daughter structure (bond) from a parent structure (pah)
! by deleting a bond between atoms: atom1 and atom2
! 
    use accuracy_m
    use structure_m
    use options_m

    implicit none
    type(structure), intent(in) :: pah
    type(structure), intent(inout) :: bond
    integer, intent(in) :: atom1, atom2

    integer :: i, j

! #####################################
! # initialize the daughter structure #
! #####################################
    bond%nat = pah%nat
    bond%order = 0
    bond%nbondlistentries = pah%nbondlistentries
    allocate(bond%neighbornumber(bond%nat))
    allocate(bond%neighborlist(bond%nat,3))
    if (bond%nbondlistentries > 0) then
        allocate(bond%bondlist(2,bond%nbondlistentries))
        bond%bondlist=pah%bondlist
    end if
    bond%neighbornumber = 0
    bond%neighborlist = 0
 

    allocate(bond%indexmapping(bond%nat))
    bond%indexmapping = pah%indexmapping
    
    if (options%print_intermediate_structures) then

        allocate(bond%doublebondlist(2,size(pah%doublebondlist,2))) 
        allocate(bond%ringlist(6,size(pah%ringlist,2))) 
        bond%doublebondnumber = pah%doublebondnumber
        bond%ringnumber = pah%ringnumber
        bond%doublebondlist = pah%doublebondlist
        bond%ringlist = pah%ringlist
        bond%hasDisconnectedParent = pah%hasDisconnectedParent
        bond%storage_unit = pah%storage_unit
    end if


! ###############################
! # fill the daughter structure #
! ###############################
    bond%neighbornumber = pah%neighbornumber
    bond%neighbornumber(atom1) = bond%neighbornumber(atom1)-1
    bond%neighbornumber(atom2) = bond%neighbornumber(atom2)-1
    bond%neighborlist = pah%neighborlist
    bond%neighborlist(atom1,1:3) = 0
    bond%neighborlist(atom2,1:3) = 0
    j = 0
    do i = 1, 3
        if (pah%neighborlist(atom1,i) /= atom2) then
            j = j+1
            bond%neighborlist(atom1,j) = pah%neighborlist(atom1,i)
        end if
    end do
    j = 0
    do i = 1, 3
        if (pah%neighborlist(atom2,i) /= atom1) then
            j = j+1
            bond%neighborlist(atom2,j) = pah%neighborlist(atom2,i)
        end if
    end do
    call clean_bond_list(bond)
    return

end subroutine create_nobond_daughter
!####################################################################################
!################## end of subroutine create_nobond_daughter ########################




!######################## subroutine create_noatoms_daughter ########################
!####################################################################################
subroutine create_noatoms_daughter(pah,pah1,nelim,delatoms,ring_exist)
!
! creates a daughter structure (pah1) from a parent structure (pah)
! by deleting nelim atoms: delatom(1),...,delatom(nelim)
! 
    use accuracy_m
    use structure_m
    use options_m
    implicit none
    type(structure), intent(in) :: pah
    type(structure), intent(inout) :: pah1
    integer, intent(in) :: nelim
    integer, intent(in) :: delatoms(nelim)
    logical, intent(in) :: ring_exist

    integer :: j,i,k,l,m

    integer, dimension(pah%nat) :: mapping
    logical, dimension(pah%nat) :: offlist


    mapping = 0
    offlist = .true.

! ########################################  
! # logical vector of atoms for deleting #
! ########################################  
    forall (i=1:nelim)
        offlist(delatoms(i))=.false.
    end forall

! #####################################
! # initialize the daughter structure #
! #####################################
    pah1%nat = pah%nat-nelim
    pah1%order = 0
    pah1%nbondlistentries = pah%nbondlistentries
    allocate(pah1%neighbornumber(pah1%nat))
    allocate(pah1%neighborlist(pah1%nat,3))
    if (pah1%nbondlistentries > 0) then
        allocate(pah1%bondlist(2,pah1%nbondlistentries))
    end if
    pah1%neighbornumber = 0
    pah1%neighborlist = 0

    allocate(pah1%indexmapping(pah%nat))
    pah1%indexmapping = pah%indexmapping

    if (options%print_intermediate_structures) then
        allocate(pah1%doublebondlist(2,size(pah%doublebondlist,2))) 
        allocate(pah1%ringlist(6,size(pah%ringlist,2))) 
        pah1%doublebondnumber = pah%doublebondnumber
        pah1%doublebondlist = pah%doublebondlist
        pah1%ringnumber = pah%ringnumber
        pah1%ringlist = pah%ringlist
        pah1%hasDisconnectedParent = pah%hasDisconnectedParent
        pah1%storage_unit = pah%storage_unit
    !###########################################
    !# add deleted atom to the happy atom list #
    !###########################################
        if (ring_exist) then 
            pah1%ringnumber = pah1%ringnumber+1
            pah1%ringlist(:,pah1%ringnumber) = pah1%indexmapping(delatoms(1:6))
        else
            j = nelim
            k = 1
            do while(j>0)
                pah1%doublebondnumber = pah1%doublebondnumber+1
                do i=1, 2
                    pah1%doublebondlist(i,pah1%doublebondnumber) = pah1%indexmapping(delatoms(k))
                    k = k + 1
                end do
                j = j - 2
            end do
        end if

    end if

! #############################
! # create the transition map #
! #############################
    j=0
    do i=1,pah%nat
        if (offlist(i)) then
            j=j+1
            mapping(i)=j
        end if
    end do

! ###########################################
! # create the structure with deleted atoms #
! ###########################################
    do i=1,pah%nat
        if (mapping(i) /= 0) then 
            k=0
            do l=1,pah%neighbornumber(i)
                m=mapping(pah%neighborlist(i,l))
                if (m /= 0) then
                    k = k + 1
                    pah1%neighborlist(mapping(i),k) = m
                end if
            end do
            pah1%neighbornumber(mapping(i)) = k
        end if
    end do

!    if (options%print_intermediate_structures) then
        do i=1,pah%nat
            if (mapping(i) /= 0) then 
                pah1%indexmapping(mapping(i)) = pah%indexmapping(i)
            end if
        end do
!    end if
        

! #######################################
! # translate to the new atom numbering #
! #######################################
    forall (i=1:pah1%nbondlistentries, j=1:2)
        pah1%bondlist(j,i)=mapping(pah%bondlist(j,i))
    end forall

    call clean_bond_list(pah1)
    return

end subroutine create_noatoms_daughter
!####################################################################################
!#################### end of subroutine create_noatoms_daughter ######################

subroutine remove_atom(pah, atom, r1)
    use structure_m

    implicit none
    type(structure), intent(inout) :: pah
    integer, intent(in) :: atom
    integer, intent(in) :: r1

    integer :: i, j, k
    integer :: n1

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
    use structure_m
    implicit none
    type(structure), intent(in) :: pah
    integer, intent(out) :: r1, r2
    integer, intent(in) :: atom1, atom2
    integer :: temp

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

!######################## subroutine cut_dangling_bonds #############################
!####################################################################################
subroutine cut_dangling_bonds(pah)
!
! removes all dangling bonds from a given polycyclic benzenoid structure (pah)
! the dangling bond is defined via a single connected atom
! a structure without singly connected atoms does not have dangling bonds
!
    use structure_m
    use options_m
    implicit none
    type(structure), intent(inout) :: pah

    integer :: atom1, atom2, r1, r2
    logical :: has_dangling_bonds
    integer :: i, j


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
            endif
            pah%indexmapping(atom1) = pah%indexmapping(r1)
            pah%indexmapping(atom2) = pah%indexmapping(r2)

            call remove_atom(pah, atom1, r1)
            call remove_atom(pah, atom2, r2)
            ! #######################################
            ! # translate to the new atom numbering #
            ! #######################################
            do i = 1, pah%nbondlistentries
                do j = 1, 2
                    if ( pah%bondlist(j,i) == atom1 .or. pah%bondlist(j,i) == atom2 ) then
                        pah%bondlist(j,i) = 0
                    else if ( pah%bondlist(j,i) == r1 ) then 
                        pah%bondlist(j,i) = atom1
                    else if ( pah%bondlist(j,i) == r2 ) then 
                        pah%bondlist(j,i) = atom2
                    end if
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
    use structure_m
    implicit none
    type(structure), intent(in) :: pah
    logical, intent(out) :: has_dangling_bonds
    integer, intent(out) :: atom1

    integer :: i

    has_dangling_bonds = .false.
    do i = 1, pah%nat
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





!########################### subroutine find_edge_ring ##############################
!####################################################################################
subroutine find_edge_ring(pah,sextet,atom1,atom2,ring_exists)
!
! for a given polycyclic benzenoid structure, the routine
! finds an aromatic ring containing atoms: atom1 and atom2
!
    use structure_m
    implicit none
    type(structure), intent(in) :: pah
    integer, intent(out), dimension(6) :: sextet
    integer, intent(in) :: atom1,atom2
    logical, intent(out) :: ring_exists
    
    integer :: i,j,k,l,m
    integer :: atom3,atom4,atom5,atom6

! ##################
! # initialization #
! ##################
    ring_exists = .false.
    sextet(1) = atom1
    sextet(2) = atom2

! #####################################
! # look for the remaining four atoms #
! #####################################
    outer: do i = 1, pah%neighbornumber(atom1)
        atom3 = pah%neighborlist(atom1,i)
        if (atom3 /= atom2) then
            do j = 1, pah%neighbornumber(atom2)
                atom4 = pah%neighborlist(atom2,j)
                if (atom4 /= atom1) then
                    do k = 1, pah%neighbornumber(atom3)
                        atom5 = pah%neighborlist(atom3,k)
                        if (atom5 /= atom1) then
                            do l = 1, pah%neighbornumber(atom4)
                                atom6 = pah%neighborlist(atom4,l)
                                if (atom6 /= atom2) then
                                    do m = 1, pah%neighbornumber(atom5)
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
    use structure_m
    use options_m
    use utils_m

    implicit none
    type(structure), intent(inout) :: pah
    integer, intent(out) :: medat

    type(structure) :: pah1
    integer :: i,j,k,lnat,start, l, n1

    logical, dimension(pah%nat+1) :: visited_list
    integer, dimension(pah%nat+1) :: map


! #########################################################
! # find the connected cluster of atoms containing atom 1 #
! #########################################################
    lnat = 0
    start = 1
    visited_list = .false.
    map = 0
    call dfs(pah, pah%nat, visited_list, lnat)

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
        do while( visited_list(i) )
            map(i) = i
            i = i + 1
        end do
        do while( .not. visited_list(j) )
            map(j) = j
            j = j - 1
        end do
        if ( i >= j ) then
            exit
        end if
!        if (options%print_intermediate_structures) then
            call swap_int(pah%indexmapping(i), pah%indexmapping(j))
!        end if
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
        map(i) = j
        map(j) = i
        i = i + 1
        j = j - 1
        if ( i >= j ) then
            exit
        end if
        
    end do

! ######################################
! # map the bond list to the new order #
! ######################################
    if (pah%nbondlistentries > 0) then
        forall (i=1:pah%nbondlistentries, j=1:2) !, pah%bondlist(j,i) /= 0)
            pah%bondlist(j,i) = map(pah%bondlist(j,i))
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
    use structure_m
    implicit none
    type(structure), intent(in) :: pah
    integer, intent(in) :: nat
    logical, dimension(nat), intent(inout) :: visit_list
    integer, intent(inout) :: lnat
    
    integer :: i
    integer, dimension(0:nat) :: stack
    integer :: cur_index

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
subroutine split_structure(pah, son1, son2, medat)
!
! split a disconnected polycyclic benzenoid structure pah into two substructures
! * son1 which is connected and contains (medat-1) atoms
! * son2 which can be connected or disconnected and contains (pah%nat-medat+1) atoms
!
! * note: only use for splitting purpose.
!

    use accuracy_m
    use structure_m
    implicit none
    
    type(structure), intent(in) :: pah
    type(structure), intent(inout) :: son1, son2
    integer, intent(in) :: medat

    integer :: i,j

! ###############################
! # allocate the son structures #
! ###############################
    son1%nat = medat-1
    son1%nbondlistentries = pah%nbondlistentries
    allocate(son1%neighbornumber(son1%nat))
    allocate(son1%neighborlist(son1%nat,3))
    son2%nat = pah%nat-medat+1
    son2%nbondlistentries = pah%nbondlistentries
    allocate(son2%neighbornumber(son2%nat))
    allocate(son2%neighborlist(son2%nat,3))
    if (pah%nbondlistentries > 0) then 
        allocate(son1%bondlist(2,son1%nbondlistentries))
        allocate(son2%bondlist(2,son2%nbondlistentries))
        son1%bondlist = pah%bondlist
        son2%bondlist = pah%bondlist
    end if

    ! ########################
    ! # update index mapping #
    ! ########################
    allocate(son1%indexmapping(son1%nat))
    son1%indexmapping = pah%indexmapping(1:medat-1)
    allocate(son2%indexmapping(son2%nat))
    son2%indexmapping = pah%indexmapping(medat:)

! #################################
! # initialize the son structures #
! #################################
    son1%neighbornumber = pah%neighbornumber(1:medat-1)
    son1%neighborlist = pah%neighborlist(1:medat-1,1:3)
    son2%neighbornumber = pah%neighbornumber(medat:pah%nat)
    son2%neighborlist = pah%neighborlist(medat:pah%nat,1:3)
    forall (i=1:son2%nat, j=1:3, son2%neighborlist(i,j) /= 0)
        son2%neighborlist(i,j) = son2%neighborlist(i,j)-son1%nat
    end forall
    forall (i=1:son1%nbondlistentries, j=1:2, son1%bondlist(j,i) > son1%nat)
        son1%bondlist(j,i) = 0
    end forall
    forall (i=1:son2%nbondlistentries, j=1:2, son2%bondlist(j,i) <= son1%nat)
        son2%bondlist(j,i) = 0
    end forall
    forall (i=1:son2%nbondlistentries, j=1:2, son2%bondlist(j,i) /= 0)
        son2%bondlist(j,i) = son2%bondlist(j,i)-son1%nat
    end forall

    call clean_bond_list(son1)
    call clean_bond_list(son2)
    
end subroutine


!######################### subroutine select_edge_bond ##############################
!####################################################################################
subroutine select_edge_bond(pah,atom1,atom2)
!
! selects two atoms (atom1 & atom2) located on the edge of a given
! polycyclic benzenoid structure pah; it uses the fact that the edge
! carbon atoms are surounded maximally by 2 hexagons
!
    use structure_m
    implicit none
    type(structure), intent(inout) :: pah
    integer, intent(out) :: atom1, atom2
    
    integer :: i,j,atom3,sextet(6)
    logical :: ring_exists,selected
    
    atom2 = 0
    selected = .false.

! #################################################
! # try to select the bond from the provided list #
! #################################################
    outer1: do i = 1, pah%nbondlistentries
        atom1 = pah%bondlist(1,i)
        atom2 = pah%bondlist(2,i)
        if (pah%neighbornumber(atom1) == 2 .or. pah%neighbornumber(atom2) == 2) then
            selected = .true.
            exit outer1
        end if
!        if (pah%neighbornumber(atom1) == 3 .and. pah%neighbornumber(atom2) == 3 )
!             cycle because filter key
!            cycle
!        end if
        do  j = 1, pah%neighbornumber(atom1)
            if (pah%neighborlist(atom1,j) == atom2) cycle
            atom3 = pah%neighborlist(atom1,j)
            call find_aromatic_sextet(pah,sextet,atom1,atom2,atom3,ring_exists)
            if (.not. ring_exists) then
                selected = .true.
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
                do j = 1, 3
                    atom2 = pah%neighborlist(atom1,min(j,mod(j,3)+1))
                    if ( pah%neighbornumber(atom2) == 3 ) then
                        cycle
                    end if
                    atom3 = pah%neighborlist(atom1,max(j,mod(j,3)+1))
                    call find_aromatic_sextet(pah,sextet,atom1,atom2,atom3,ring_exists)
                    if (.not. ring_exists) then
                        exit outer
                    end if
                end do
                atom2 = 0
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





!######################## subroutine find_aromatic_sextet ###########################
!####################################################################################
subroutine find_aromatic_sextet(pah,sextet,atom1,atom2,atom3,ring_exists)
!
! for a given polycyclic benzenoid structure, the routine
! finds an aromatic ring containing atoms: atom1, atom2, and atom3
!
    use structure_m
    implicit none
    type(structure), intent(in) :: pah
    integer, intent(out), dimension(6) :: sextet
    logical, intent(out) :: ring_exists
    integer, intent(in) :: atom1, atom2, atom3

    integer(kint) :: i,j,k,l
    integer(kint) :: atom4, atom5, atom6

! ##################
! # initialization #
! ##################
    ring_exists = .false.
    sextet(1) = atom1
    sextet(2) = atom2
    sextet(3) = atom3

! ######################################
! # look for the remaining three atoms #
! ######################################
    do i = 1, pah%neighbornumber(atom3)
        atom4 = pah%neighborlist(atom3,i)
        if (atom4 /= atom1) then
            do j = 1, pah%neighbornumber(atom2)
                atom5 = pah%neighborlist(atom2,j)
                if (atom5 /= atom1) then
                    do k = 1, pah%neighbornumber(atom4)
                        atom6 = pah%neighborlist(atom4,k)
                        do l = 1, pah%neighbornumber(atom5)
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


!####################### subroutine find_all_hexagons ###############################
!####################################################################################
subroutine find_all_hexagons(nat,pah,nhex,lista)
! 
! find a list of all hexagons in a polycyclic structure pah
!
    use accuracy_m
    use structure_m

    implicit none
    integer, intent(in) :: nat
    type(structure), intent(in) :: pah
    integer, intent(out) :: nhex ! number of hexagons in structure pah
    integer, intent(inout), dimension(6,nat) :: lista
    
    integer :: i,j,k,l,atom2,atom3
    integer, dimension(6) :: sextet
    logical :: ring_exists

    nhex=0
! #######################
! # loop over all atoms #
! #######################
    atomloop: do i = 1, pah%nat

!   ##################################
!   # if atom i has only 2 neighbors #
!   ##################################
        if (pah%neighbornumber(i) == 2) then
            if (pah%neighborlist(i,1) < i) cycle atomloop
            if (pah%neighborlist(i,2) < i) cycle atomloop
            call find_aromatic_sextet(pah,sextet,i,pah%neighborlist(i,2),pah%neighborlist(i,1),ring_exists)
            if (ring_exists) then
                do j = 4, 6
                    if (sextet(j) < i) cycle atomloop
                end do
                nhex = nhex+1
                do l = 1, 6
                    lista(l,nhex) = sextet(l)
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
                atom2 = pah%neighborlist(i,mod(j,3)+1)
                if (atom2 < i) cycle innerloop
                atom3 = pah%neighborlist(i,mod(j+1,3)+1)
                if (atom3 < i) cycle innerloop
                call find_aromatic_sextet(pah, sextet, i, atom2, atom3, ring_exists)
                if (ring_exists) then
                    do k = 4, 6
                        if (sextet(k) < i) cycle innerloop
                    end do
                    nhex = nhex + 1
                    do l = 1, 6
                        lista(l,nhex) = sextet(l)
                    end do
                end if
            end do innerloop
        end if
    end do atomloop
    return

end subroutine find_all_hexagons
!####################################################################################
!#################### end of subroutine find_all_hexagons ###########################

end module operator_m
