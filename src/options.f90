module options_m
    implicit none
    save

    type :: options_t
        logical :: print_intermediate_structures
        logical :: force_print_structures
        integer :: print_order
        character(len=80) :: bondlistfile
        logical :: has_bondlistfile
        logical :: verbose 
        logical :: decompose_print
        logical :: testrun
!        logical :: use_database
!        logical :: create_database
!        character(len=80) :: databasefile

    end type

    type(options_t) :: options    
contains
    subroutine initialize_options()
        options%print_intermediate_structures = .false.
        options%print_order = 0
        options%has_bondlistfile = .false.
        options%force_print_structures = .false.
        options%verbose = .false.
        options%testrun = .false.               
        options%decompose_print = .false.
!        options%use_database = .false.
!        options%create_database = .false.
    end subroutine

end module
