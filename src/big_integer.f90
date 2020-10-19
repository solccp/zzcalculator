!################################ module types_module ###############################
!####################################################################################
module big_integer_m
    use accuracy_m
    implicit none


    integer, parameter :: entry_power = (range(1_kint)+1)/2
    integer, parameter :: block_size = ceiling( real(desired_digits) / real(entry_power))

    integer(kint), parameter :: block_max = 10**entry_power

    type, public :: big_integer
        sequence
        integer(kint) :: leadpow
        integer(kint) :: tabl(block_size)
    end type big_integer


contains

function setvli(a) result(c)   
    implicit none
    integer(kint),intent(in) :: a
    integer(kint) :: b,i
    type(big_integer) :: c

    c%tabl = 0
    b = a
    do i = 1, block_size
        c%tabl(i) = mod(b,block_max)
        b = (b-c%tabl(i)) / block_max
        if (b == 0) then
            c%leadpow = i
            exit
        end if
    end do

end function setvli


function addvli(a,b) result(c)
    implicit none
    type(big_integer),intent(in) :: a,b
    type(big_integer) :: c
    integer(kint) :: i,val

    c%tabl = 0
    c%leadpow = max(a%leadpow, b%leadpow)
    do i = 1, c%leadpow
        c%tabl(i) = c%tabl(i)+a%tabl(i)+b%tabl(i)
    end do
    do i = 1, min(block_size-1, c%leadpow)
        if (c%tabl(i) >= block_max) then
            val = mod(c%tabl(i),block_max)
            c%tabl(i+1) = c%tabl(i+1) + (c%tabl(i))/block_max
            c%tabl(i) = val
        end if
    end do

    if (c%tabl(block_size) >= block_max ) then
        print*,"big integer overflow in addvli, increase [desired_digits] in accuracy.f90"
        stop
    end if
 
    if (c%leadpow < block_size) then  
        if (c%tabl(c%leadpow+1) /= 0) c%leadpow = c%leadpow+1
    end if

end function addvli
 
function multvli(a,b) result(c)
    implicit none
    type(big_integer),intent(in) :: a,b
    type(big_integer) :: c
    integer(kint) :: i,j,val

    c%tabl = 0
    do i=1, a%leadpow
        do j=1, min(block_size-i,b%leadpow)
            c%tabl(i+j-1) = c%tabl(i+j-1) + a%tabl(i)*b%tabl(j)
        end do
    end do
    do i=1, min(block_size,a%leadpow+b%leadpow)
        if (c%tabl(i) >= block_max) then
            val = mod(c%tabl(i),block_max)
            c%tabl(i+1) = c%tabl(i+1) + (c%tabl(i))/block_max
            c%tabl(i) = val
        end if
    end do

    if (c%tabl(block_size) >= block_max) then
        print*,"big integer overflow in multvli, increase [desired_digits] in accuracy.f90"
        stop
    end if

    c%leadpow = 0
    do i = min(block_size, a%leadpow+b%leadpow+1), 1, -1
        if (c%tabl(i) /=0) then
            c%leadpow = i
            exit
        end if
    end do
 
end function multvli
 
subroutine print_vli_in_string(pos,string,val)
!   prints integer val in the string at position pos
    use accuracy_m
    implicit none
    integer :: i
    integer, intent(inout) :: pos
    type(big_integer) :: val
    character(len=str_len_long) :: string
    character(len=10) :: fmt1

    do i=val%leadpow, 1, -1
        if (i == val%leadpow) then
            write(string(pos:),'(I0)') val%tabl(i)
        else
            write(fmt1, '(a,i0,a,i0,a)') '(I0', entry_power, '.', entry_power, ')'
            write(string(pos:),fmt1) val%tabl(i)
        end if
        pos = pos + len(trim(string(pos:pos+entry_power)))
    end do
    return

end subroutine print_vli_in_string

end module big_integer_m

!####################################################################################
!############################ end of module types_module ############################
