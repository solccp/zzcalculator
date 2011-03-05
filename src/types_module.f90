!################################ module types_module ###############################
!####################################################################################
module types_module

    integer, parameter :: kint = 8
    integer, parameter :: kreal = kind(0.0d0)

    type,public :: vlonginteger
        integer,public :: leadpow
        integer(kint), private :: tabl(40)
    end type vlonginteger


contains

function setvli(a) result(c)   
    implicit none
    integer(kint),intent(in) :: a
    integer(kint) :: b,i
    type(vlonginteger) :: c

    c%tabl=0
    b=a
    do i=1,40
        c%tabl(i)=mod(b,10000)
        b=(b-c%tabl(i))/10000
        if (b == 0) then
            c%leadpow=i
            exit
        end if
    end do

end function setvli


function addvli(a,b) result(c)
    implicit none
    type(vlonginteger),intent(in) :: a,b
    type(vlonginteger) :: c
    integer(kint) :: i,val

    c%tabl=0
    c%leadpow=max(a%leadpow,b%leadpow)
    do i=1,c%leadpow
        c%tabl(i)=c%tabl(i)+a%tabl(i)+b%tabl(i)
    end do
    do i=1,min(39,c%leadpow)
        if (c%tabl(i) >= 10000) then
            val=mod(c%tabl(i),10000)
            c%tabl(i+1)=c%tabl(i+1)+(c%tabl(i)-val)/10000
            c%tabl(i)=val
        end if
    end do

    if (c%tabl(40) >= 10000) then
        print*,"overflow in addvli, enlarge size of tabl"
        stop
    end if
 
    if (c%leadpow < 40) then  
        if (c%tabl(c%leadpow+1) /= 0) c%leadpow=c%leadpow+1
    end if

end function addvli
 
function multvli(a,b) result(c)
    implicit none
    type(vlonginteger),intent(in) :: a,b
    type(vlonginteger) :: c
    integer(kint) :: i,j,val

    c%tabl=0
    do i=1,a%leadpow
        do j=1,min(40-i,b%leadpow)
            c%tabl(i+j-1)=c%tabl(i+j-1)+a%tabl(i)*b%tabl(j)
        end do
    end do
    do i=1,min(39,a%leadpow+b%leadpow)
        if (c%tabl(i) >= 10000) then
            val=mod(c%tabl(i),10000)
            c%tabl(i+1)=c%tabl(i+1)+(c%tabl(i)-val)/10000
            c%tabl(i)=val
        end if
    end do

    if (c%tabl(40) >= 10000) then
        print*,"overflow in multvli, enlarge size of tabl"
        stop
    end if

    c%leadpow=0
    do i=min(40,a%leadpow+b%leadpow+1),1,-1
        if (c%tabl(i) /=0) then
            c%leadpow=i
            exit
        end if
    end do
 
end function multvli
 
subroutine printvli(a)
    implicit none
    type(vlonginteger),intent(in) :: a
    integer(kint) :: i,j,val

    if (a%leadpow == 0) then
        write(*,*)"0"
    else
        write(*,'(1X,I4,40I4.4)')(a%tabl(i),i=a%leadpow,1,-1)
    end if

end subroutine printvli


subroutine print_vli_in_string(pos,string,val)
!   prints integer val in the string at position pos
    implicit none
    integer(kint) :: i
    integer(kint), intent(inout) :: pos
    type(vlonginteger) :: val
    character(len=500) :: string

    do i=val%leadpow,val%leadpow
        select case (val%tabl(i))
        case (0:9)
            write(string(pos:),'(i1)')val%tabl(i)
            pos=pos+1
        case (10:99)
            write(string(pos:),'(i2)')val%tabl(i)
            pos=pos+2
        case (100:999)
            write(string(pos:),'(i3)')val%tabl(i)
            pos=pos+3
        case default
            write(string(pos:),'(i4)')val%tabl(i)
            pos=pos+4
        end select
    end do
    do i=val%leadpow-1,1,-1
        write(string(pos:),'(i4.4)')val%tabl(i)
        pos=pos+4
    end do
    return

end subroutine print_vli_in_string

end module types_module

!####################################################################################
!############################ end of module types_module ############################
