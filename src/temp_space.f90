module temp_space
    use types_module
    implicit none
    save
    integer(kint), dimension(:), allocatable, target :: int_1darray_1
    integer(kint), dimension(:), allocatable, target :: int_1darray_2
    logical, dimension(:), allocatable, target :: bool_1darray_1
contains
    subroutine initialize_temp_space(nat)
        integer(kint), intent(in) :: nat
        allocate(int_1darray_1(nat))
        allocate(int_1darray_2(nat))
        allocate(bool_1darray_1(nat))
    end subroutine
    subroutine finalize_temp_space
        if (allocated(int_1darray_1)) then
            deallocate(int_1darray_1)
        end if
        if (allocated(int_1darray_2)) then
            deallocate(int_1darray_2)
        end if
        if (allocated(bool_1darray_1)) then    
            deallocate(bool_1darray_1)
        end if
    end subroutine
end module
