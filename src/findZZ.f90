!###################### subroutine find_ZZ_polynomial ###############################
!####################################################################################
subroutine find_ZZ_polynomial(pah,level)
! 
! find resursively the ZZ polynomial for the structure pah
!
    use types_module
    use structure_module
    use output
    use options_m
    use hash_m

    implicit none
    type(structure), intent(inout) :: pah
    integer(kint), intent(in) :: level

    integer(kint) :: medat
    logical :: hit

    if ( pah%nat >= 6 ) then
        call get_hash(pah, pah%hash_key)
        call get_polynomial(pah, hit)
        if (hit) then
!            print *, 'got a hit in database, nat=', pah%nat
            return
        end if
    end if

    


! ###########################
! # if pah contains 0 atoms #
! ###########################
    if (pah%nat == 0) then
        pah%order = 0
        allocate(pah%polynomial(pah%order+1))
        pah%polynomial(1) = setvli(1_kint)

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
    else if (mod(pah%nat,2_kint) == 1) then
        pah%order = 0
        allocate(pah%polynomial(pah%order+1))
        pah%polynomial(1)=setvli(0_kint)

! ##########################################
! # check if the graph of pah is connected #
! ##########################################
    else
        call check_if_connected(pah,medat)

!   #############################################
!   # if pah is connected, decompose it further #
!   #############################################
        if (medat == 0) then
            call decompose(pah,level)

!   ########################################################################
!   # if pah is disconnected, split it and decompose the fragments further #
!   ########################################################################
        else
            call split_and_decompose(pah,medat,level)
        end if
        if ( pah%nat >= 6 .and. medat == 0) then
            call add_polynomial(pah)
        end if
    end if

    return

end subroutine find_ZZ_polynomial
!####################################################################################
!################## end of subroutine find_ZZ_polynomial ############################


