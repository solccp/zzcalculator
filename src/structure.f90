module structure_module
    
    use types_module
    type,public :: structure
        integer(kint) :: nat
        integer(kint),allocatable,dimension(:) :: initiallabel
        integer(kint),allocatable,dimension(:) :: neighbornumber
        integer(kint),allocatable,dimension(:,:) :: neighborlist
        integer(kint),allocatable,dimension(:) :: indexmapping
        integer(kint) :: order
        type(vlonginteger),allocatable,dimension(:) :: polynomial
        integer(kint), allocatable,dimension(:,:) :: bondlist
        integer(kint) :: nbondlistentries
        integer(kint) :: doublebondnumber
        integer(kint), allocatable, dimension(:,:) :: doublebondlist
        integer(kint) :: ringnumber
        integer(kint), allocatable, dimension(:,:) :: ringlist
        logical :: hasDisconnectedParent
        integer :: storage_unit
    end type structure

contains 
    subroutine create(pah, pah1)
        type(structure), intent(in) :: pah
        type(structure), intent(inout) :: pah1
        allocate(pah1%neighbornumber(pah1%nat))
        allocate(pah1%neighborlist(pah1%nat,3))
        if (pah1%nbondlistentries > 0) then
            allocate(pah1%bondlist(2,pah1%nbondlistentries))
        end if
        pah1%neighbornumber = 0
        pah1%neighborlist = 0
        allocate(pah1%indexmapping(size(pah%indexmapping)))
        pah1%indexmapping = pah%indexmapping
        allocate(pah1%doublebondlist(2,size(pah%doublebondlist,2))) 
        allocate(pah1%ringlist(6,size(pah%ringlist,2))) 
        pah1%hasDisconnectedParent = pah%hasDisconnectedParent
        pah1%storage_unit = pah%storage_unit
    end subroutine
    subroutine destory(pah)
        type(structure), intent(inout) :: pah
        if (allocated(pah%initiallabel)) then
            deallocate(pah%initiallabel)
        end if
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
        if (allocated(pah%doublebondlist)) then
            deallocate(pah%doublebondlist)
        end if
        if (allocated(pah%ringlist)) then
            deallocate(pah%ringlist)
        end if
    end subroutine

end module
