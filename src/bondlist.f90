
module bondlist_m
    use accuracy_m
    use structure_m


    real(kreal), private, dimension(:,:), pointer :: mod_geom

contains

function dis(atom1, atom2) result(res)
    double precision :: res
    integer, intent(in) :: atom1, atom2
    
    double precision :: v(3)

    v = (mod_geom(:, atom1) + mod_geom(:, atom2))/2.0d0
    res = sqrt(dot_product(v,v))
   
end function

subroutine sort(pah, bondlist, len)
    type(structure), intent(in) :: pah
    integer, intent(inout), dimension(:,:) :: bondlist
    integer, intent(in) :: len
    integer :: i, j
    double precision :: D_vec(len)
    double precision :: d_tmp
    integer :: ia_tmp(2)
    

    do i = 1, len
        D_vec(i) = dis(pah%indexmapping(bondlist(1, i)), pah%indexmapping(bondlist(2,i)))
    end do

    do i = len, 2, -1
        do j = 1, i - 1
            if (D_vec(j) < D_vec(j+1)) then
                d_tmp = D_vec(j)
                D_vec(j) = D_vec(j+1)
                D_vec(j+1) = d_tmp
                ia_tmp = bondlist(:,j)
                bondlist(:,j) = bondlist(:,j+1)
                bondlist(:,j+1) = ia_tmp
            end if
        end do
    end do
end subroutine

function index_a(array, item) result(res)
    integer, intent(in), dimension(:) :: array
    integer, intent(in):: item
    integer :: res, i
    res = 0
    do i = 1, size(array)
        if ( array(i) == item )  then
            res = i
            exit
        end if
    end do
end function

subroutine bipartition(pah, bondlist, nlist, vex_g1, len_g1, vex_g2, len_g2, succ)
    use lapack
    implicit none
    type(structure), intent(in) :: pah

    integer, intent(inout), dimension(:,:) :: bondlist
    integer, intent(out) :: nlist
    logical, intent(out) :: succ

    integer, intent(inout), dimension(:) :: vex_g1, vex_g2
    integer, intent(out) :: len_g1, len_g2

    logical, dimension(pah%nat) :: g1, g2

    double precision , dimension(:,:), allocatable :: LP_mat
    double precision , dimension(:), allocatable :: ev  

    integer :: i, j, k, tmp
    integer :: nver

    nver = pah%nat

    allocate(LP_mat(nver, nver))
    allocate(ev(nver))
    LP_mat = 0.0d0

    do i = 1, nver
        LP_mat(i,i) = pah%neighbornumber(i)
        do k = 1, pah%neighbornumber(i)
            j = pah%neighborlist(i, k)
            LP_mat(i, j) = -1
        end do
    end do

    call syev('V', 'L', nver, LP_mat, nver, ev)
    if ( abs(ev(2)) < 1.0d-12 ) then
        succ = .false.
    else
        succ = .true.
        g1 = .false.
        g2 = .false.
        do i = 1, nver
            if (LP_mat(i,2)> 0) then
                g1(i) = .true.
            else
                g2(i) = .true.
            end if
        end do

        k = 0
        do i = 1, nver
            do j = 1, pah%neighbornumber(i)
                tmp = pah%neighborlist(i, j)
                if ( g1(i) .and. g2(tmp) ) then
                    k = k + 1
                    bondlist(1, k) = i
                    bondlist(2, k) = tmp
                end if
            end do      
        end do
        nlist = k

        !sort!
        call sort(pah, bondlist, nlist)

        k = 0
        j = 0
        do i = 1, nver
            if ( g1(i) ) then
                j = j + 1
                vex_g1(j) = i
            else if ( g2(i) ) then
                k = k + 1
                vex_g2(k) = i
            end if
        end do


        len_g1 = j
        len_g2 = k
    end if

    
    deallocate(LP_mat)
    deallocate(ev)

end subroutine


subroutine build_bondlist(pah, geom, map)
    use options_m
    use operator_m
    implicit none
    type(structure), intent(inout), target :: pah
    real(kreal), intent(in), dimension(:,:), target :: geom 
    integer, intent(in), dimension(:) :: map

    integer, dimension(2,pah%nat) :: bondlist

    integer :: i, j, k

    integer :: vertexes(pah%nat)
    
    integer :: nlist, nlist_total

    integer :: vex_g1(pah%nat), vex_g2(pah%nat), len_g1, len_g2

    type(structure_ptr), dimension(pah%nat) :: queue
    type(structure), pointer :: pah_ptr, cur_pah
    
    integer :: queue_index
    logical :: first, succ


    mod_geom => geom

    do i = 1, pah%nat
        vertexes(i) = i
    end do

    nlist_total = 0

    queue_index = 1
    queue(queue_index)%ptr => pah

    first = .true.
    do while (queue_index > 0)
        cur_pah => queue(queue_index)%ptr
        queue_index = queue_index -1

        call bipartition(cur_pah, bondlist(:,nlist_total+1:), nlist, vex_g1, len_g1, vex_g2, len_g2, succ)
        if ( .not. succ ) then
            exit
        end if
        do i = nlist_total+1, nlist_total + nlist
            bondlist(:, i) = map(cur_pah%indexmapping(bondlist(:, i)))
        end do

        nlist_total = nlist_total + nlist

        nullify(pah_ptr)
        allocate(pah_ptr)
        call create_noatoms_daughter(cur_pah, pah_ptr, len_g1, vex_g1, .false.) 
!        call cut_dangling_bonds(pah_ptr)
        if ( pah_ptr%nat > 50 ) then
            queue_index = queue_index + 1
            queue(queue_index)%ptr => pah_ptr
        else
            deallocate(pah_ptr)
        end if
        nullify(pah_ptr)
        allocate(pah_ptr)
        call create_noatoms_daughter(cur_pah, pah_ptr, len_g2, vex_g2, .false.) 
!        call cut_dangling_bonds(pah_ptr)
        if ( pah_ptr%nat > 30 ) then
            queue_index = queue_index + 1
            queue(queue_index)%ptr => pah_ptr
        else
            deallocate(pah_ptr)
        end if
        if ( first ) then
            first = .false.
        else
            deallocate(cur_pah)
            nullify(cur_pah)
        end if
    end do
       
    pah%nbondlistentries = nlist_total
    allocate(pah%bondlist(2, nlist_total))
    pah%bondlist = bondlist(:, :nlist_total)


    if ( options%verbose ) then
        print *, "------------------------------"
        print *, "generated bondlist"
    

        do i = 1, nlist_total
            print *, pah%bondlist(:, i)
        end do

        print *, "------------------------------"
    end if

end subroutine

end module
