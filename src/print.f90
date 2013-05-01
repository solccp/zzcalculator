!######################### subroutine print_ZZ_polynomial ###########################
!####################################################################################
subroutine print_ZZ_polynomial(pah)
!
! prints the computed ZZ polynomial of a given structure pah;
! the polynomial is first printed in a formated manner to a string
! and later the chunks of the string longer than 80 characters are
! displayed at stdout
!
    use accuracy_m
    use structure_m
    use big_integer_m
    use options_m
    implicit none
    type(structure), intent(in) :: pah
    integer :: i,cpos
    type(big_integer) :: total
    character(len=40000) :: finalZZpolynomial


    if ( options%verbose ) then
        write(*, '(1x,a)') 'Final ZZ_polynomial:'
    end if

! #########################
! # initialize the string #
! #########################
    finalZZpolynomial = ''
    cpos = 1
    total = pah%polynomial(1)
    call print_vli_in_string(cpos,finalZZpolynomial,pah%polynomial(1))
    if (pah%order > 0) then
        write(finalZZpolynomial(cpos:),*)'+ '
        cpos = cpos+3
        call print_vli_in_string(cpos,finalZZpolynomial,pah%polynomial(2))
        write(finalZZpolynomial(cpos:),*)'x'
        cpos = cpos+2
        total = addvli(pah%polynomial(1),pah%polynomial(2))
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
subroutine print_ZZ_polynomial_simple(pah)
    use accuracy_m
    use structure_m
    use big_integer_m
    use options_m
    implicit none
    type(structure), intent(in) :: pah
    integer :: i,cpos
    type(big_integer) :: total
    character(len=40000) :: finalZZpolynomial


! #########################
! # initialize the string #
! #########################


    total = setvli(0_kint)

    do i = 0, pah%order
        finalZZpolynomial = ''
        cpos = 1
        call print_vli_in_string(cpos,finalZZpolynomial,pah%polynomial(i+1))
        write(*,'(i0,1x,a)') i, trim(finalZZpolynomial)
        total = addvli(total,pah%polynomial(i+1))
    end do

    finalZZpolynomial=''
    cpos = 1
    call print_vli_in_string(cpos,finalZZpolynomial,total)
    write(*,'(2a)')"t:",trim(finalZZpolynomial)

    return

end subroutine print_ZZ_polynomial_simple
subroutine print_ZZ_polynomial_XML(pah)

    use accuracy_m
    use structure_m
    use big_integer_m
    use options_m
    implicit none
    type(structure), intent(in) :: pah
    integer :: i,cpos
    type(big_integer) :: total
    character(len=40000) :: finalZZpolynomial


! #########################
! # initialize the string #
! #########################


    total = setvli(0_kint)
    write(*, '(a)') '<zzpolynomial>'

    do i = 0, pah%order
        write(*,'(2x,a)') '<term>'
        finalZZpolynomial = ''
        cpos = 1
        call print_vli_in_string(cpos,finalZZpolynomial,pah%polynomial(i+1))
        write(*,'(4x,a,i0,a)') '<order>', i, '</order>'
        write(*,'(4x,a,a,a)')  '<coefficient>', trim(finalZZpolynomial), '</coefficient>'
        total = addvli(total,pah%polynomial(i+1))
        write(*,'(2x,a)') '</term>'
    end do

    finalZZpolynomial=''
    cpos = 1
    call print_vli_in_string(cpos,finalZZpolynomial,total)
    write(*,'(2x,a,a,a)') '<total>',trim(finalZZpolynomial), '</total>'
    write(*, '(a)') '</zzpolynomial>'
    return

end subroutine print_ZZ_polynomial_XML
!####################################################################################





!##################### subroutine print_int_in_string ###############################
!####################################################################################
subroutine print_int_in_string(pos,string,val)
!
! prints integer val in the string at position pos
!
    use accuracy_m
    use structure_m
    implicit none
    integer, intent(in) :: val
    integer, intent(inout) :: pos
    character(len=*), intent(inout) :: string
    integer, parameter :: int_len = range(val)

    write(string(pos:),'(i0)')val
    pos = pos + len(trim(string(pos:pos+int_len)))
    return
end subroutine print_int_in_string
!####################################################################################
!################## end of subroutine print_int_in_string ###########################
