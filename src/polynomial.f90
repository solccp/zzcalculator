!######################### subroutine sum_polynomials ###############################
!####################################################################################
subroutine sum_polynomials(pah,daughter1,daughter2,daughter3,ring_exists)
! 
! obtain the ZZ polynomial of the parent structure (pah) 
! by summing  the ZZ polynomials for three daughter structures
!
!  ZZ(pah) = ZZ(daughter1) + ZZ(daughter2) + x * ZZ(daughter3)
!
  use types_module
    use structure_module
  implicit none
  integer(kint) :: i
  type(structure) :: pah,daughter1,daughter2,daughter3
  logical :: ring_exists

! ###################################
! # initialize parent ZZ polynomial #
! ###################################
  if (ring_exists) then
    pah%order=max0(daughter1%order,daughter2%order,daughter3%order+1)
  else
    pah%order=max0(daughter1%order,daughter2%order)
  end if
  allocate(pah%polynomial(pah%order+1))
  pah%polynomial=setvli(0)

! #########################################################
! # add the contribution from daughter structure 1 (bond) #
! #########################################################
  do i=0,daughter1%order
!    pah%polynomial(i+1)=pah%polynomial(i+1)+daughter1%polynomial(i+1)
    pah%polynomial(i+1)=addvli(pah%polynomial(i+1),daughter1%polynomial(i+1))
  end do

! ############################################################
! # add the contribution from daughter structure 2 (corners) #
! ############################################################
  do i=0,daughter2%order
!    pah%polynomial(i+1)=pah%polynomial(i+1)+daughter2%polynomial(i+1)
    pah%polynomial(i+1)=addvli(pah%polynomial(i+1),daughter2%polynomial(i+1))
  end do

! #########################################################
! # add the contribution from daughter structure 3 (ring) #
! #########################################################
  if (ring_exists) then
    do i=0,daughter3%order
!      pah%polynomial(i+2)=pah%polynomial(i+2)+daughter3%polynomial(i+1)
      pah%polynomial(i+2)=addvli(pah%polynomial(i+2),daughter3%polynomial(i+1))
    end do
  end if
  return

end subroutine sum_polynomials
!####################################################################################
!###################### end of subroutine sum_polynomials ###########################



!######################### subroutine multiply_polynomials ##########################
!####################################################################################
subroutine multiply_polynomials(pah,son1,son2)
!
! compute the ZZ polynomial of the parent structure pah 
! by multiplying the ZZ polynomial of the son structures: son1 & son2
!
  use types_module
    use structure_module
  implicit none
  integer(kint) :: deg1,deg2,i,j
  type(structure) :: pah,son1,son2

! #####################################################################
! # find the highest non-vanishing power of the ZZ polynomial of son1 #
! #####################################################################
  deg1=-1
  do i=son1%order,0,-1
!    if (son1%polynomial(i+1) == 0) cycle
    if (son1%polynomial(i+1)%leadpow == 0) cycle
    deg1=i
    exit
  end do

! #####################################################################
! # find the highest non-vanishing power of the ZZ polynomial of son2 #
! #####################################################################
  deg2=-1
  do i=son2%order,0,-1
!    if (son2%polynomial(i+1) == 0) cycle
    if (son2%polynomial(i+1)%leadpow == 0) cycle
    deg2=i
    exit
  end do

! ##################################################################
! # check if any of the ZZ polynomials for son structures vanished #
! ##################################################################
  if (deg1 == -1 .or. deg2 == -1 )then
    pah%order=0
    allocate(pah%polynomial(pah%order+1))
    pah%polynomial(1)=setvli(0)
!    pah%polynomial=0

! #######################################################
! # allocate the ZZ polynomial for the parent structure #
! #######################################################
  else
    pah%order=deg1+deg2
    allocate(pah%polynomial(pah%order+1))
    pah%polynomial=setvli(0)

!   #####################################################
!   # multiply the ZZ polynomials of the son structures # 
!   #####################################################
    do i=0,deg1
      do j=0,deg2
        pah%polynomial(i+j+1)=addvli(pah%polynomial(i+j+1),multvli(son1%polynomial(i+1),son2%polynomial(j+1)))
      end do
    end do
  end if
  return

end subroutine multiply_polynomials
!####################################################################################
!###################### end of subroutine multiply_polynomials ######################
