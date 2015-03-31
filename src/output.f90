module output_m
    use input_m
    implicit none
    character(len=100), save :: output_filename = 'intermediate_strs.yaml'

    integer, save :: output_unit = -1

contains
    subroutine close_file()
        if ( output_unit /= -1 ) then
            close(output_unit) 
        end if
    end subroutine
    subroutine open_file()
        output_unit = getunit()
        open(file=trim(output_filename), unit=output_unit) 
    end subroutine
    subroutine write_connections(pah)
    
        use structure_m
        use options_m
        type(structure), intent(in) :: pah
        logical, parameter :: flow = .true.
        integer :: i,j
        integer, save :: output_index = 0


        if (pah%ringnumber < options%print_order) then
            return 
        end if

        if (output_unit < 0) then
            call open_file()
        end if

        output_index = output_index +1

        write(output_unit,'(a)') '---'
        write(output_unit,'(a,I0)') 'structure_id: ', output_index
        write(output_unit,'(a,I0)') 'number_of_atoms: ', pah%doublebondnumber*2+pah%ringnumber*6
        if (pah%doublebondnumber > 0) then
            write(output_unit,'(a)') 'double_bonds:'
            if (flow) then
                do i=1, pah%doublebondnumber
                    write(output_unit, '(a, 2(I0, a))') '- [', pah%doublebondlist(1,i), ', ', pah%doublebondlist(2,i), ']'
                end do
            else
                do i=1, pah%doublebondnumber
                    write(output_unit, '(a, I0)') '- - ', pah%doublebondlist(1,i)
                    write(output_unit, '(a, I0)') '  - ', pah%doublebondlist(2,i)
                end do
            end if
        end if
        if (pah%ringnumber > 0) then
            write(output_unit,'(a)') 'rings:'
            if (flow) then
                do i=1, pah%ringnumber
                    write(output_unit, '(a, 6(I0, a))') '- [', (pah%ringlist(j,i), ', ', j=1,5), pah%ringlist(6,i), ']'
                end do
            else
                do i=1, pah%ringnumber
                    write(output_unit,'(a, I0)') '- - ', pah%ringlist(1,i)
                    write(output_unit,'(a, I0)') '  - ', pah%ringlist(2,i)
                    write(output_unit,'(a, I0)') '  - ', pah%ringlist(3,i)
                    write(output_unit,'(a, I0)') '  - ', pah%ringlist(4,i)
                    write(output_unit,'(a, I0)') '  - ', pah%ringlist(5,i)
                    write(output_unit,'(a, I0)') '  - ', pah%ringlist(6,i)
                end do
            end if
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
    subroutine write_connections_partial(pah)
        use structure_m
        type(structure), intent(in) :: pah
        logical :: inuse
        integer :: i


        INQUIRE(UNIT=pah%storage_unit, OPENED = inuse)
        if (.not. inuse) then
            print *, 'something wrong in write_connections_partial'
            stop
        end if

        write(pah%storage_unit, '(2(I6,1X))') pah%doublebondnumber, pah%ringnumber
        do i = 1, pah%doublebondnumber
            write(pah%storage_unit, '(2(I6,1X))') pah%doublebondlist(:,i)
        end do
        do i = 1, pah%ringnumber
            write(pah%storage_unit, '(6(I6,1X))') pah%ringlist(:,i)
        end do
    end subroutine

    subroutine combine_connection_output(pah, son1, son2)
        use structure_m
        type(structure), intent(in) :: pah, son1, son2
        type(structure) :: temp
        integer :: ndb1, nring1
        integer :: ndb2, nring2
        integer :: i,j
        integer :: indexes(6)
        integer :: errorcode

        allocate(temp%doublebondlist(2,size(pah%doublebondlist,2)))
        allocate(temp%ringlist(6,size(pah%ringlist,2)))
        temp%storage_unit = pah%storage_unit
        temp%doublebondlist = pah%doublebondlist
        temp%ringlist = pah%ringlist
        temp%ringnumber = pah%ringnumber
        temp%doublebondnumber = pah%doublebondnumber

        rewind(son1%storage_unit)
        do
            temp%ringnumber = pah%ringnumber
            temp%doublebondnumber = pah%doublebondnumber
            read(son1%storage_unit,*, iostat=errorcode) ndb1, nring1
            if (errorcode /= 0) exit
            do i=1, ndb1
                read(son1%storage_unit,*) indexes(1:2)
                temp%doublebondnumber = temp%doublebondnumber+1
                temp%doublebondlist(:,temp%doublebondnumber) = indexes(1:2)
            end do
            do i=1, nring1
                read(son1%storage_unit,*) indexes(1:6)
                temp%ringnumber = temp%ringnumber+1
                temp%ringlist(:,temp%ringnumber) = indexes(1:6)
            end do
            rewind(son2%storage_unit)
            do
                temp%doublebondnumber = pah%doublebondnumber + ndb1
                temp%ringnumber = pah%ringnumber + nring1
                read(son2%storage_unit,*, iostat=errorcode) ndb2, nring2
                if (errorcode /= 0) exit
                do i=1, ndb2
                    read(son2%storage_unit,*) indexes(1:2)
                    temp%doublebondnumber = temp%doublebondnumber+1
                    temp%doublebondlist(:,temp%doublebondnumber) = indexes(1:2)
                end do
                do i=1, nring2
                    read(son2%storage_unit,*) indexes(1:6)
                    temp%ringnumber = temp%ringnumber+1
                    temp%ringlist(:,temp%ringnumber) = indexes(1:6)
                end do
                if (pah%hasDisconnectedParent) then
                    call write_connections_partial(temp)
                else
                    call write_connections(temp)
                end if        
            end do
            100 continue
        end do

        deallocate(temp%doublebondlist)
        deallocate(temp%ringlist)

    end subroutine

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
