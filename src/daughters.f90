!##################### subroutine create_nobond_daughter ############################
!####################################################################################
subroutine create_nobond_daughter(pah,bond,atom1,atom2)
!
! creates a daughter structure (bond) from a parent structure (pah)
! by deleting a bond between atoms: atom1 and atom2
! 
    use types_module
    use structure_module
    use options_m

    implicit none
    type(structure), intent(in) :: pah
    type(structure), intent(inout) :: bond
    integer(kint), intent(in) :: atom1, atom2

    integer(kint) :: i, j

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
 

    if (options%print_intermediate_structures) then
        allocate(bond%indexmapping(bond%nat))
        bond%indexmapping = pah%indexmapping

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
    bond%neighbornumber=pah%neighbornumber
    bond%neighbornumber(atom1)=bond%neighbornumber(atom1)-1
    bond%neighbornumber(atom2)=bond%neighbornumber(atom2)-1
    bond%neighborlist=pah%neighborlist
    bond%neighborlist(atom1,1:3)=0
    bond%neighborlist(atom2,1:3)=0
    j=0
    do i=1,3
        if (pah%neighborlist(atom1,i) /= atom2) then
            j=j+1
            bond%neighborlist(atom1,j)=pah%neighborlist(atom1,i)
        end if
    end do
    j=0
    do i=1,3
        if (pah%neighborlist(atom2,i) /= atom1) then
            j=j+1
        bond%neighborlist(atom2,j)=pah%neighborlist(atom2,i)
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
    use types_module
    use structure_module
    use options_m
    use temp_space
    implicit none
    type(structure), intent(in) :: pah
    type(structure), intent(inout) :: pah1
    integer(kint), intent(in) :: nelim
    integer(kint), intent(in) :: delatoms(nelim)
    logical, intent(in) :: ring_exist

    integer(kint) :: j,i,k,l,m

!####################################
!# int_1darray_1 use as mapping     #
!# bool_1darray_1 use as offlist    #
!####################################

    int_1darray_1(1:pah%nat)=0
    bool_1darray_1(1:pah%nat)=.true.

! ########################################  
! # logical vector of atoms for deleting #
! ########################################  
    forall (i=1:nelim)
        bool_1darray_1(delatoms(i))=.false.
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


    if (options%print_intermediate_structures) then
        allocate(pah1%indexmapping(pah%nat))
        pah1%indexmapping = pah%indexmapping
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
        if (bool_1darray_1(i)) then
            j=j+1
            int_1darray_1(i)=j
        end if
    end do

! ###########################################
! # create the structure with deleted atoms #
! ###########################################
    do i=1,pah%nat
        if (int_1darray_1(i) /= 0) then 
            k=0
            do l=1,pah%neighbornumber(i)
                m=int_1darray_1(pah%neighborlist(i,l))
                if (m /= 0) then
                    k = k + 1
                    pah1%neighborlist(int_1darray_1(i),k) = m
                end if
            end do
            pah1%neighbornumber(int_1darray_1(i)) = k
        end if
    end do
    if (options%print_intermediate_structures) then
        do i=1,pah%nat
            if (int_1darray_1(i) /= 0) then 
                pah1%indexmapping(int_1darray_1(i)) = pah%indexmapping(i)
            end if
        end do
    end if
        

! #######################################
! # translate to the new atom numbering #
! #######################################
    forall (i=1:pah1%nbondlistentries, j=1:2)
        pah1%bondlist(j,i)=int_1darray_1(pah%bondlist(j,i))
    end forall

    call clean_bond_list(pah1)
    return

end subroutine create_noatoms_daughter
!####################################################################################
!#################### end of subroutine create_noatoms_daughter ######################
