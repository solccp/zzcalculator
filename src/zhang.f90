!######################### program zhang_polynomial #################################
!####################################################################################
program zhang_polynomial
!
! This program calculates the Zhang-Zhang polynomial for benzenoid structures;
! the ZZ polynomial contains beside other quantities also the Clar number,
! Clar count, and Kekule number of a given structure
!   
! Reference:
!   I. Gutman, B. Furtula, and A. T. Balaban
!   Polycyclic Aromatic Compounds 26 pp.17-35, 2006
!
    use types_module
    use structure_module
    implicit none
    integer(kint) :: i,nhex,level=0
    type(structure) :: pah


    integer(kint) :: argc
    integer(kint) :: input_unit
    character(len=256) :: input_fname

    argc = command_argument_count()

    if (argc /= 1) then
        stop 'abnormal termination: wrong arguments'
    end if

    call get_command_argument(1, input_fname)

! ############################################################
! # read initial geometry data and create topological matrix #
! ############################################################
    call read_input(input_fname, pah)

! #############################################################
! # find recursively the ZZ polynomial of the given structure #
! #############################################################
    call find_ZZ_polynomial(pah,level)

! ###########################
! # print the ZZ polynomial #
! ###########################
    call print_ZZ_polynomial(pah)

end
!####################################################################################
!###################### end of program zhang_polynomial #############################
