!###################### subroutine find_ZZ_polynomial ###############################
!####################################################################################
recursive subroutine find_ZZ_polynomial(pah,level)
! 
! find resursively the ZZ polynomial for the structure pah
!
  use types_module
  use output
  implicit none
  integer(kint) :: medat,level
  type(structure) :: pah

! ###########################
! # if pah contains 0 atoms #
! ###########################
  if (pah%nat == 0) then
    pah%order=0
    allocate(pah%polynomial(pah%order+1))
    pah%polynomial(1)=setvli(1)
    if (pah%hasDisconnectedParent) then
        call write_connections_partial(pah)
    else
        call write_connections(pah)
    end if

! ##########################################
! # if pah contains an odd number of atoms #
! ##########################################
  else if (mod(pah%nat,2) == 1) then
    pah%order=0
    allocate(pah%polynomial(pah%order+1))
    pah%polynomial(1)=setvli(0)

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
      print *, '*****************disconnected'
      call split_and_decompose(pah,medat,level)
    end if
  end if
  return

end subroutine find_ZZ_polynomial
!####################################################################################
!################## end of subroutine find_ZZ_polynomial ############################
