module options_m
    implicit none
    save

    type :: options_t
        logical :: print_intermediate_structures
        logical :: force_print_structures
        integer :: print_order
        character(len=100) :: bondlistfile
        logical :: has_bondlistfile
        logical :: verbose 
        logical :: decompose_print
        logical :: testrun
        logical :: use_connection_file
        character(len=100) :: connection_file
        logical :: use_bipartition
        logical :: kekule_only
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
        options%use_connection_file = .false.
        options%use_bipartition = .false.
        options%kekule_only = .false.
!        options%use_database = .false.
!        options%create_database = .false.
    end subroutine

end module
