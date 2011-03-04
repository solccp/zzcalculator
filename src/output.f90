module output

contains
    subroutine write_connections(pah)
        use types_module
        type(structure), intent(in) :: pah
        integer :: i
        integer, save :: output_index = 0
        
        output_index = output_index +1
        write(*,*) 'writing structure ', output_index
        write(*,*) 'number of atoms: ', pah%doublebondnumber*2+pah%ringnumber*6
        if (pah%doublebondnumber > 0) then
            write(*,*) 'double bond:'
            do i=1, pah%doublebondnumber
                write(*,*) pah%doublebondlist(:,i)
            end do
        end if
        if (pah%ringnumber > 0) then
            write(*,*) 'ring:'
            do i=1, pah%ringnumber
                write(*,'(999(I,X))') pah%ringlist(:,i)
            end do
        end if
        write(*,*) "===================================================="
    end subroutine
    function getunit() 
        integer :: getunit 
        logical :: inuse 
        inuse = .true. 
        getunit = 0
        do while (inuse == .true.) 
            getunit = getunit + 1 
            INQUIRE(UNIT = getunit, OPENED = inuse) 
        end do 
    end function
    subroutine write_connections_partial(pah)
        use types_module
        type(structure), intent(in) :: pah
        logical :: inuse
        integer :: i


        INQUIRE(UNIT=pah%storage_unit, OPENED = inuse)
        if (.not. inuse) then
            print *, 'something wrong in write_connections_partial'
            stop
        end if

        write(pah%storage_unit, '(2(I6,X))') pah%doublebondnumber, pah%ringnumber
        do i=1, pah%doublebondnumber
            write(pah%storage_unit, '(2(I6,X))') pah%doublebondlist(:,i)
        end do
        do i=1, pah%ringnumber
            write(pah%storage_unit, '(6(I6,X))') pah%ringlist(:,i)
        end do
    end subroutine

    subroutine combine_connection_output(pah, son1, son2)
        use types_module
        type(structure), intent(in) :: pah, son1, son2
        type(structure) :: temp
        integer :: ndb, nring
        integer :: i,j, k1, k2
        integer :: indexes(6)
        integer :: errorcode

        allocate(temp%doublebondlist(2,size(pah%doublebondlist,2)))
        allocate(temp%ringlist(2,size(pah%ringlist,2)))
        temp%storage_unit = pah%storage_unit
        temp%doublebondlist = pah%doublebondlist
        temp%ringlist = pah%ringlist
        temp%ringnumber = pah%ringnumber
        temp%doublebondnumber = pah%doublebondnumber

        k1 = 1
        rewind(son1%storage_unit)
        do
            temp%ringnumber = pah%ringnumber
            temp%doublebondnumber = pah%doublebondnumber
            read(son1%storage_unit,*, iostat=errorcode) ndb, nring
            if (errorcode /= 0) exit
            do i=1, ndb
                read(son1%storage_unit,*) indexes(1:2)
                temp%doublebondnumber = temp%doublebondnumber+1
                temp%doublebondlist(:,temp%doublebondnumber) = indexes(1:2)
            end do
            do i=1, nring
                read(son1%storage_unit,*) indexes(1:6)
                temp%ringnumber = temp%ringnumber+1
                temp%ringlist(:,temp%ringnumber) = indexes(1:6)
            end do
            k2 = 1
            rewind(son2%storage_unit)
            do
                temp%doublebondnumber = pah%doublebondnumber + ndb
                temp%ringnumber = pah%ringnumber + nring
                read(son2%storage_unit,*,end=100) ndb, nring
                if (errorcode /= 0) exit
                print *,'k1,k2', k1,k2
                do i=1, ndb
                    read(son2%storage_unit,*) indexes(1:2)
                    temp%doublebondnumber = temp%doublebondnumber+1
                    temp%doublebondlist(:,temp%doublebondnumber) = indexes(1:2)
                end do
                do i=1, nring
                    read(son2%storage_unit,*) indexes(1:6)
                    temp%ringnumber = temp%ringnumber+1
                    temp%ringlist(:,temp%ringnumber) = indexes(1:6)
                end do
                if (pah%hasDisconnectedParent) then
                    call write_connections_partial(temp)
                else
                    call write_connections(temp)
                end if        
                k2 = k2 + 1
            end do
            100 continue
            k1 = k1 + 1
        end do

        deallocate(temp%doublebondlist)
        deallocate(temp%ringlist)

    end subroutine

end module
