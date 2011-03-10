module options_m
    implicit none
    save

    type :: options_t
        logical :: print_intermediate_structures
        integer :: print_order
    end type

    type(options_t) :: options    
contains
    subroutine initialize_options()
        options%print_intermediate_structures = .false.
        options%print_order = 0
    end subroutine

end module
