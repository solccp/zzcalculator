module options_m
    implicit none
    save

    type :: options_t
        character(len=100) :: bondlistfile
        logical :: has_bondlistfile
        logical :: verbose
        logical :: use_bipartition
        logical :: kekule_only
		logical :: read_connection_table
		logical :: simple_printing

    end type

    type(options_t) :: options
contains
    subroutine initialize_options()
        options%has_bondlistfile = .false.
        options%verbose = .false.
        options%use_bipartition = .false.
        options%kekule_only = .false.
        options%read_connection_table = .false.
        options%simple_printing = .false.
    end subroutine

end module
