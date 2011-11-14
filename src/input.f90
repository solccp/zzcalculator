
module input_m
    use accuracy_m
    use structure_m
    use options_m
    use utils_m
    use operator_m
    use partition_m
 
    implicit none

    real(kreal), allocatable, dimension(:,:), save :: ori_geom

contains

!############################ subroutine read_input #################################
!####################################################################################
subroutine read_input(input_fname,pah)
!
! read the geometry of a polycyclic benzenoid structure pah from the file 'geometry',
! filter out only the carbon atoms, and create the topological matrix for it
!
    character(len=*), intent(in) :: input_fname
    type(structure), intent(inout) :: pah

    integer :: info
    integer :: cnat, status, bnat, i, j, k, nhex, l, m, errorcode, a1, a2
    integer, allocatable, dimension(:,:) :: lista
    integer, allocatable, dimension(:,:) :: localbondlist
    character(len=2) :: atname
    real(kreal) :: inertia(3,3), eival(3), work(100)
    real(kreal), allocatable, dimension(:,:) :: geom
    integer, allocatable, dimension(:) :: map
    integer, allocatable, dimension(:) :: rmap
    real(kreal), dimension(3) :: x
    logical :: inlist, bondfileexists

    ! ######################
    ! # read geometry file #
    ! ######################
    open(20,file=trim(input_fname),status='old')
    read(20,*) bnat
    read(20,*)
    allocate(geom(3,bnat))
    allocate(ori_geom(3,bnat))
    allocate(map(bnat))
    allocate(rmap(bnat))
    cnat = 0
    map = 0
    rmap = 0
    do i = 1, bnat
        read(20,*) atname,(x(j),j=1,3)
        if (atname == 'C' .or. atname == 'c') then
            cnat = cnat+1
            map(i) = cnat
            rmap(cnat) = i
            geom(:,cnat) = x(:)
        end if
        ori_geom(:,i) = x(:)
    end do
    close(20)


    ! ######################
    ! # allocate structure #
    ! ######################
    pah%nat = cnat
    pah%order = 0
    allocate(pah%neighbornumber(pah%nat))
    allocate(pah%neighborlist(pah%nat,3))
    pah%neighbornumber=0

    allocate(pah%indexmapping(pah%nat))
    do i=1, pah%nat
        pah%indexmapping(i) = rmap(i)
    end do

    if (options%print_intermediate_structures) then
        pah%doublebondnumber = 0
        allocate(pah%doublebondlist(2,pah%nat))
        pah%ringnumber = 0
        allocate(pah%ringlist(6,pah%nat))
    end if


    bondfileexists = .false.
    if ( options%use_connection_file) then
        !use connection file
        open(unit=99, file=trim(options%connection_file))
        do
            read(99,*,iostat=errorcode) a1, a2
            if ( errorcode /= 0 ) then
                exit
            end if
            bondfileexists = .true.
            if ( a1 > bnat ) then
                write(*,'(a,i0,a)') "Atom index ", a1, " is greater than the number of atoms, check your connection file"
                stop
            end if
   
            if ( map(a1) == 0 .or. map(a2) == 0 ) then
                cycle
            end if 
            if (pah%neighbornumber(map(a1)) == 3) then
                write(*,'(a,i0,a)') "Atom # ", a1, " has too many neighbors, check your connection file"
                stop
            end if
            pah%neighbornumber(map(a1))=pah%neighbornumber(map(a1))+1
            pah%neighborlist(map(a1),pah%neighbornumber(map(a1)))=map(a2)        
            pah%neighbornumber(map(a2))=pah%neighbornumber(map(a2))+1
            pah%neighborlist(map(a2),pah%neighbornumber(map(a2)))=map(a1)        
        end do  
        close(99)
    end if


    if (.not. bondfileexists ) then
    ! #######################
        ! # find neighbor table #
        ! #######################
        do i=1, cnat
            do j=i+1, cnat
                if (dist(cnat,i,j,geom) < ccdist) then
                    pah%neighbornumber(i)=pah%neighbornumber(i)+1
                    pah%neighborlist(i,pah%neighbornumber(i))=j
                    pah%neighbornumber(j)=pah%neighbornumber(j)+1
                    pah%neighborlist(j,pah%neighbornumber(j))=i
                end if
            end do
        end do
    end if





    ! ####################################################
    ! # read (if provided) the preferred partition order #
    ! ####################################################
    pah%nbondlistentries = 0
    if (options%has_bondlistfile) then
        inquire(file=trim(options%bondlistfile),exist=bondfileexists)
        if (bondfileexists) then
            allocate(localbondlist(2,3*cnat/2))
            open(21,file=trim(options%bondlistfile))
            do
                read(21,*,iostat=errorcode) a1,a2
                if (errorcode == 0) then
                    pah%nbondlistentries =pah%nbondlistentries + 1
                    localbondlist(1,pah%nbondlistentries) = map(a1)
                    localbondlist(2,pah%nbondlistentries) = map(a2)
                    if (a1 > bnat .or. a2 > bnat) then
                        write(*,*)"Ooops, looks like your file: bondlist"
                        write(*,*)"does not correspond to your input file"
                        write(*,*)"cnat",cnat," a1,a2",a1,a2
                        stop
                    end if
                else if (errorcode == -1) then
                    exit
                else
                    write(*,*)"Ooops, reading error from the file: bondlist"
                    write(*,*)"Verify if the file is not corrupted"
                    stop
                end if
            end do
            allocate(pah%bondlist(2,pah%nbondlistentries))
            pah%bondlist = localbondlist(:,1:pah%nbondlistentries)
            deallocate(localbondlist)
            close(21)
            call clean_bond_list(pah)
        end if
    end if
!#### test

    if (pah%nbondlistentries == 0) then
        call bipartition(pah)
    end if


!    call print_bondlist(pah)
    ! #########################################################
    ! # find all substructures with removed one aromatic ring #
    ! #########################################################
!    allocate(lista(6,pah%nat))
!    call find_all_hexagons(pah%nat,pah,nhex,lista)
!    open(22,file='new.geometries')
!    do i=1,nhex
!        inertia=0.0d0
!        do l=1,pah%nat
!            inlist=.false.
!            do m=1,6
!                if (l == lista(m,i)) inlist=.true.
!            end do
!            if (inlist) cycle
!            do j=1,3
!                do k=1,3
!                    inertia(j,k)=inertia(j,k)+geom(j,l)*geom(k,l)
!                end do
!            end do
!        end do
!        write(22,'(i6)')pah%nat
!        write(22,'(1x,4i6,3f15.3)')pah%nat,0,0,0,inertia(1,1)+inertia(2,2)+inertia(3,3)
!        do j=1,pah%nat
!            inlist=.false.
!            do k=1,6
!                if (j == lista(k,i)) inlist=.true.
!            end do
!            if (inlist) then
!                write(22,'(a,3f20.10)')"B",(geom(k,j),k=1,3) 
!            else
!                write(22,'(a,3f20.10)')"C",(geom(k,j),k=1,3) 
!            end if
!        end do
!    end do
!    close(22)

    deallocate(map)
    deallocate(rmap)
    return

end subroutine read_input
!####################################################################################
!######################### end of subroutine read_input #############################


end module input_m
