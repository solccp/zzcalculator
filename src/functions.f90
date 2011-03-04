!############################### function dist ######################################
!####################################################################################
real(kreal) function dist(nat,atom1,atom2,geom)
!
! compute the distance between atom1 & atom2
!
  use types_module
  implicit none
  integer(kint),intent(in) :: nat                    ! number of atoms
  integer(kint),intent(in) :: atom1                  ! position of atom 1
  integer(kint),intent(in) :: atom2                  ! position of atom 2
  integer(kint) :: i                                 ! local counter
  real(kreal),intent(in),dimension(3,nat) :: geom    ! geometry table
  real(kreal) :: r,x                                 ! local variables

  r=0.0d0
  do i=1,3
    x=geom(i,atom1)-geom(i,atom2)
    r=r+x*x
  end do
  dist=sqrt(r)
  return

end function dist
!####################################################################################
!############################### end of function dist ###############################


