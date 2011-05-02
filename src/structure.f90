module structure_module
    
    use iso_varying_string
    use types_module
    type, public :: structure
        integer(kint) :: nat
        integer(kint),allocatable,dimension(:) :: neighbornumber
        integer(kint),allocatable,dimension(:,:) :: neighborlist
        integer(kint) :: order
        type(vlonginteger),allocatable,dimension(:) :: polynomial
        integer(kint), allocatable,dimension(:,:) :: bondlist
        integer(kint) :: nbondlistentries

        integer(kint),allocatable,dimension(:) :: indexmapping
        integer(kint) :: doublebondnumber
        integer(kint), allocatable, dimension(:,:) :: doublebondlist
        integer(kint) :: ringnumber
        integer(kint), allocatable, dimension(:,:) :: ringlist
        logical :: hasDisconnectedParent
        integer :: storage_unit
    end type structure

    type ::  tree_node
        type(varying_string) :: key
        logical :: hasChild = .false.
        type(structure), pointer :: pah => NULL()
        
        type(tree_node), pointer :: child_corners => NULL()
        type(tree_node), pointer :: child_bond => NULL()   
        type(tree_node), pointer :: child_ring => NULL()   
        type(tree_node), pointer :: child_son1 => NULL()   
        type(tree_node), pointer :: child_son2 => NULL()   
    end type

    type :: tree_node_ptr
        type(tree_node), pointer :: node
    end type
contains 
    subroutine destory(pah)
        use options_m
        type(structure), intent(inout) :: pah
        if (allocated(pah%neighbornumber)) then
            deallocate(pah%neighbornumber)
        end if
        if (allocated(pah%neighborlist)) then
            deallocate(pah%neighborlist)
        end if
        if (allocated(pah%indexmapping)) then
            deallocate(pah%indexmapping)
        end if
        if (allocated(pah%polynomial)) then
            deallocate(pah%polynomial)
        end if
        if (allocated(pah%bondlist)) then
            deallocate(pah%bondlist)
        end if

        if (options%print_intermediate_structures) then

            if (allocated(pah%doublebondlist)) then
                deallocate(pah%doublebondlist)
            end if
            if (allocated(pah%ringlist)) then
                deallocate(pah%ringlist)
            end if
        end if
    end subroutine

end module
