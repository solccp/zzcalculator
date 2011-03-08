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