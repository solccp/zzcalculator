module structure_m
    use accuracy_m
    use big_integer_m
    use iso_varying_string

    type, public :: structure
        integer :: nat
        integer, allocatable, dimension(:) :: neighbornumber
        integer, allocatable, dimension(:,:) :: neighborlist
        integer :: order
        type(big_integer), allocatable, dimension(:) :: polynomial
        integer, allocatable, dimension(:,:) :: bondlist
        integer :: nbondlistentries = 0

        logical :: polynomial_computed = .false.

        integer, allocatable, dimension(:) :: indexmapping
        integer :: doublebondnumber
        integer, allocatable, dimension(:,:) :: doublebondlist
        integer :: ringnumber
        integer, allocatable, dimension(:,:) :: ringlist
        logical :: hasDisconnectedParent = .false.
        integer :: storage_unit = 0

        logical :: hasDanglingBond = .false.
        integer :: nextDanglingBond = 0


    end type structure

    type structure_ptr
        type(structure), pointer :: ptr
    end type

    type ::  pah_tree_node
        integer :: shared_count = 0
        logical :: hasChild = .false.
        type(structure), pointer :: pah => NULL()
        
        type(pah_tree_node), pointer :: parent => NULL()
        type(pah_tree_node), pointer :: child_corners => NULL()
        type(pah_tree_node), pointer :: child_bond => NULL()   
        type(pah_tree_node), pointer :: child_ring => NULL()   
        type(pah_tree_node), pointer :: child_son1 => NULL()   
        type(pah_tree_node), pointer :: child_son2 => NULL()   
    end type

    type :: pah_tree_node_ptr
        type(pah_tree_node), pointer :: ptr
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
