!######################### subroutine print_ZZ_polynomial ###########################
!####################################################################################
subroutine print_ZZ_polynomial(pah)
!
! prints the computed ZZ polynomial of a given structure pah;
! the polynomial is first printed in a formated manner to a string
! and later the chunks of the string longer than 80 characters are
! displayed at stdout
!
    use types_module
    use structure_module
    implicit none
    type(structure), intent(in) :: pah
    integer(kint) :: i,cpos
    type(vlonginteger) :: total
    character(len=500) :: finalZZpolynomial

! #########################
! # initialize the string #
! #########################
    finalZZpolynomial=''
    cpos=1
    call print_vli_in_string(cpos,finalZZpolynomial,pah%polynomial(1))
    if (pah%order > 0) then
        write(finalZZpolynomial(cpos:),*)'+ '
        cpos=cpos+3
        call print_vli_in_string(cpos,finalZZpolynomial,pah%polynomial(2))
        write(finalZZpolynomial(cpos:),*)'x'
        cpos=cpos+2
        total=addvli(pah%polynomial(1),pah%polynomial(2))
    end if
    if (pah%order <= 1) then
        write(*,*)finalZZpolynomial(1:cpos-1)
    end if

! ##############################################
! # loop over all degrees of the ZZ polynomial #
! ##############################################
    do i=2,pah%order
        total=addvli(total,pah%polynomial(i+1))
        write(finalZZpolynomial(cpos:),*)'+ '
        cpos=cpos+3
        call print_vli_in_string(cpos,finalZZpolynomial,pah%polynomial(i+1))
        write(finalZZpolynomial(cpos:),*)'x^'
        cpos=cpos+3
        call print_int_in_string(cpos,finalZZpolynomial,i)

!   ####################################################
!   # flush out the chunks of the polynomial to stdout #
!   ####################################################
        if (i == pah%order) then
            write(*,'(1x,a)')finalZZpolynomial(1:cpos-1)
        else if (cpos > 450) then
            write(*,'(1x,a)')finalZZpolynomial(1:cpos-1)
            finalZZpolynomial=''
            cpos=1
        end if
    end do
    finalZZpolynomial=''
    cpos = 1
    call print_vli_in_string(cpos,finalZZpolynomial,total)
    write(*,'(1x,2a)')"total: ",trim(finalZZpolynomial)

    return

end subroutine print_ZZ_polynomial
!####################################################################################
!################### end of subroutine print_ZZ_polynomial ##########################



!##################### subroutine print_int_in_string ###############################
!####################################################################################
subroutine print_int_in_string(pos,string,val)
!
! prints integer val in the string at position pos
!
    use types_module
    use structure_module
    implicit none
    integer(kint) :: val,pos
    character(len=500) :: string
    integer, parameter :: int_len = range(val)

    write(string(pos:),'(i0)')val
    pos = pos + len(trim(string(pos:pos+int_len)))
    return
end subroutine print_int_in_string
!####################################################################################
!################## end of subroutine print_int_in_string ###########################
