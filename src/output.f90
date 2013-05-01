module output_m
    use input_m
    implicit none
    integer, save :: output_unit = -1

contains
    subroutine close_file()
        if ( output_unit /= -1 ) then
            close(output_unit) 
        end if
    end subroutine
    function getunit() 
        integer :: getunit 
        logical :: inuse 
        inuse = .true. 
        getunit = 10
        do while (inuse) 
            getunit = getunit + 1 
            INQUIRE(UNIT = getunit, OPENED = inuse) 
        end do 
    end function

    subroutine write_xyz(pah, filename)
        use structure_m
        type(structure), intent(in) :: pah
        character(len=*), intent(in) :: filename
        integer :: i

        open(unit=99, file=trim(filename))
        call write_xyz_unit(pah, 99, '')
        close(99) 
    end subroutine
    subroutine write_xyz_unit(pah, unitnum, title)
        use structure_m
        type(structure), intent(in) :: pah
        integer, intent(in) :: unitnum
        character(len=*), intent(in) :: title
        integer :: i

        if ( pah%nat == 0 ) then
            return 
        end if
        write(unitnum, '(i0)') pah%nat
        write(unitnum, '(a)') trim(title)
        do i=1, pah%nat
            write(unitnum, '(a, 3(2x, F12.6))') 'C', ori_geom(:,pah%indexmapping(i))
        end do

    end subroutine



end module
