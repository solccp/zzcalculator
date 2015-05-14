

!###################### subroutine find_ZZ_polynomial ###############################
!####################################################################################
subroutine find_ZZ_polynomial(pah,level)
! 
! find resursively the ZZ polynomial for the structure pah
!
    use accuracy_m
    use structure_m
    use output_m
    use options_m
    use polynomial_m
    use decompose_m

    implicit none
    type(structure), intent(inout) :: pah
    integer, intent(in) :: level

    integer :: medat


!    if ( options%decompose_print ) then
!        if ( level == level_to_print ) then
!            call write_xyz_unit(pah, unitnum, level, 'Find_ZZ')
!            return
!        end if
!    end if


! ###########################
! # if pah contains 0 atoms #
! ###########################
    if (pah%nat == 0) then
        call set_polynomial(pah, 1_kint)
        if (options%print_intermediate_structures) then
            if (pah%hasDisconnectedParent) then
                call write_connections_partial(pah)
            else
                call write_connections(pah)
            end if
        end if

! ##########################################
! # if pah contains an odd number of atoms #
! ##########################################
    else if ( mod(pah%nat,2) == 1) then
        call set_polynomial(pah, 0_kint)

! ##########################################
! # check if the graph of pah is connected #
! ##########################################
    else
        call check_if_connected(pah,medat)

!   #############################################
!   # if pah is connected, decompose it further #
!   #############################################
        if (medat == 0) then
            call decompose(pah, level)

!   ########################################################################
!   # if pah is disconnected, split it and decompose the fragments further #
!   ########################################################################
        else
            call split_and_decompose(pah, medat, level)
        end if

        pah%polynomial_computed = .true.
    end if

    return

end subroutine find_ZZ_polynomial
!####################################################################################
!################## end of subroutine find_ZZ_polynomial ############################


