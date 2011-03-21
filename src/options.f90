module options_m
    implicit none
    save

    type :: options_t
        logical :: print_intermediate_structures
        integer :: print_order
        character(len=80) :: bondlistfile
        logical :: has_bondlistfile
    end type

    type(options_t) :: options    
contains
    subroutine initialize_options()
        options%print_intermediate_structures = .false.
        options%print_order = 0
        options%has_bondlistfile = .false.
    end subroutine

end module
