!################### subroutine find_all_kekule structures ##########################
!####################################################################################
recursive subroutine find_all_kekule_structures(pah,nstr)
! 
! find a number of all possible Kekule structures for a polycyclic structure pah
!
  use types_module
  implicit none
  integer(kint) :: atom1,atom2,dnstr,i,medat,j,k
  integer(kint) :: nstr ! number of Kekule structures for pah
  type(structure) :: pah,daughter

! ##############################
! # if no atoms left, return 1 #
! ##############################
  if (pah%nat == 0) then
    nstr=1

! ############################################
! # if an odd number of atoms left, return 0 #
! ############################################
  else if (mod(pah%nat,2) == 1) then
    nstr=0
  else 

!   ###################################################
!   # if disconnected, treat each fragment separately #
!   ###################################################
    call check_if_connected(pah,medat)
    if (medat /=0) then
      call split_and_kekule(pah,medat,nstr)      

!   ##############################################
!   # if connected, follow the usual kekule road #
!   ##############################################
    else

!     #############################################
!     # select an edge bond between atom1 & atom2 #
!     #############################################
      call select_edge_bond(pah,atom1,atom2)

!     #######################################################
!     # loop over two possible locations of the double bond #
!     #######################################################
      nstr=0
      do i=1,2
        atom2=pah%neighborlist(atom1,i)
        call create_nobond_daughter(pah,daughter,atom1,atom2)
        call cut_dangling_bonds(daughter)
        call find_all_kekule_structures(daughter,dnstr)
        nstr=nstr+dnstr
        deallocate(daughter%neighborlist)
        deallocate(daughter%neighbornumber)
      end do
    end if
  end if
  return

end subroutine find_all_kekule_structures
!####################################################################################
!################## end of subroutine find_all_kekule structures ####################


!######################## subroutine split_and_kekule ###############################
!####################################################################################
recursive subroutine split_and_kekule(pah,medat,nstr)
!
! split a disconnected polycyclic benzenoid structure pah into two substructures
! * son1 which is connected and contains (medat-1) atoms
! * son2 which can be connected or disconnected and contains (pah%nat-medat+1) atoms
! find kekule number for both structures and multiply it into the kekule number of the parent
!
  use types_module
  implicit none
  integer(kint) :: medat,i,j,nstr,sstr1,sstr2
  type(structure) :: pah,son1,son2

! ###############################
! # allocate the son structures #
! ###############################
  son1%nat=medat-1
  allocate(son1%neighbornumber(son1%nat))
  allocate(son1%neighborlist(son1%nat,3))
  son2%nat=pah%nat-medat+1
  allocate(son2%neighbornumber(son2%nat))
  allocate(son2%neighborlist(son2%nat,3))

! #################################
! # initialize the son structures #
! #################################
  son1%neighbornumber=pah%neighbornumber(1:medat-1)
  son1%neighborlist=pah%neighborlist(1:medat-1,1:3)
  son2%neighbornumber=pah%neighbornumber(medat:pah%nat)
  son2%neighborlist=pah%neighborlist(medat:pah%nat,1:3)
  forall (i=1:son2%nat, j=1:3, son2%neighborlist(i,j) /= 0)
    son2%neighborlist(i,j)=son2%neighborlist(i,j)-son1%nat
  end forall

! ##############################################
! # find Kekule number for both son structures #
! ##############################################
  call find_all_kekule_structures(son1,sstr1)
  call find_all_kekule_structures(son2,sstr2)

! ##################################################
! # multiply Kekule numbers of both son structures #
! ##################################################
  nstr=sstr1*sstr2

! #################################
! # deallocate the son structures #
! #################################
  deallocate(son1%neighbornumber)
  deallocate(son1%neighborlist)
  deallocate(son2%neighbornumber)
  deallocate(son2%neighborlist)
  return

end subroutine split_and_kekule
!####################################################################################
!##################### end of subroutine split_and_kekule ###########################
