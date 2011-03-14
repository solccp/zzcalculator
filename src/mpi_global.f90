
module mpi_global
    use mpi
    implicit none
    save
    logical, allocatable, dimension(:) :: image_busy
    integer :: image_count
    integer :: image_id 
    integer :: mpi_vli_type, mpi_row_type
 
contains
    subroutine mpi_global_init()
        use types_module
        integer :: error
    
        call MPI_Comm_rank ( MPI_COMM_WORLD, image_id, error )
        call MPI_Comm_size ( MPI_COMM_WORLD, image_count, error )
        allocate(image_busy(0:image_count-1))
        image_busy = .false.

        call MPI_TYPE_CONTIGUOUS(block_size+1, MPI_INTEGER8, mpi_vli_type, error)
        call MPI_TYPE_COMMIT(mpi_vli_type, error)

        call MPI_TYPE_VECTOR(3, 1, 3, MPI_INTEGER8, mpi_row_type, error)
        call MPI_TYPE_COMMIT(mpi_row_type, error)

    end subroutine

    subroutine mpi_global_finalize()
        integer :: error
    
        if (allocated(image_busy)) then
           deallocate(image_busy)
        end if
        call MPI_TYPE_FREE(mpi_vli_type, error)
        call MPI_TYPE_FREE(mpi_row_type, error)
    end subroutine

    subroutine mpi_update_states()
    end subroutine

    function get_free_node() result(res)
        integer :: res
        integer :: i
        res = -1
        do i = 1, image_count-1
            if ( .not. image_busy(i) ) then
                res = i
                exit
            end if
        end do
    end function

    subroutine send_structure(pah, level, dest_node)
        use types_module
        use structure_module
        integer, intent(in) :: dest_node
        integer(kint), intent(in) :: level
        type(structure), intent(in) :: pah
        integer :: error, i

        print*, image_id, 'sending structure' 
        call MPI_SEND(level, 1, MPI_INTEGER8, dest_node, 1, MPI_COMM_WORLD, ERROR)
        call MPI_SEND(pah%nat, 1, MPI_INTEGER8, dest_node, 2, MPI_COMM_WORLD, ERROR)
        call MPI_SEND(pah%neighbornumber, pah%nat, MPI_INTEGER8, dest_node, 3, MPI_COMM_WORLD, ERROR)
        do i = 1, 3
            call MPI_SEND(pah%neighborlist(1,i), pah%nat, MPI_INTEGER8, dest_node, 4, MPI_COMM_WORLD, ERROR)
        end do
        call MPI_SEND(pah%nbondlistentries, 1, MPI_INTEGER8, dest_node, 5, MPI_COMM_WORLD, ERROR)
        if ( pah%nbondlistentries > 0 ) then
            call MPI_SEND(pah%bondlist, pah%nbondlistentries, MPI_INTEGER8, dest_node, 6, MPI_COMM_WORLD, ERROR)
        end if

        image_busy(dest_node) = .true.
        print*, image_id, 'finished sending structure' 
    end subroutine

    subroutine recv_structure(pah, level)
        use types_module
        use structure_module
        type(structure), intent(inout) :: pah
        integer(kint), intent(out) :: level
        integer :: mpi_status(MPI_STATUS_SIZE) 
        integer :: error, i

        print*, image_id, 'getting structure' 

        call MPI_RECV(level, 1, MPI_INTEGER8, 0, 1, MPI_COMM_WORLD, mpi_status, ERROR)
        call MPI_RECV(pah%nat, 1, MPI_INTEGER8, 0, 2, MPI_COMM_WORLD, mpi_status, ERROR)
        call MPI_RECV(pah%neighbornumber, pah%nat, MPI_INTEGER8, 0, 3, MPI_COMM_WORLD, mpi_status, ERROR)
        do i = 1, 3
            call MPI_RECV(pah%neighborlist(1,i), pah%nat, MPI_INTEGER8, 0, 4, MPI_COMM_WORLD, mpi_status, ERROR)
        end do
        call MPI_RECV(pah%nbondlistentries, 1, MPI_INTEGER8, 0, 5, MPI_COMM_WORLD, mpi_status, ERROR)
        if ( pah%nbondlistentries > 0 ) then
            call MPI_RECV(pah%bondlist, pah%nbondlistentries, MPI_INTEGER8, 0, 6, MPI_COMM_WORLD, mpi_status, ERROR)
        end if
        print*, image_id, 'finished getting structure' 
    end subroutine

    subroutine send_polynomial(pah)
        use types_module
        use structure_module
        type(structure), intent(in) :: pah
        integer :: error

        print*, image_id, 'sending polynomial' 
        call MPI_SEND(pah%order, 1, MPI_INTEGER8, 0, 1, MPI_COMM_WORLD, ERROR)
        call MPI_SEND(pah%polynomial, pah%order+1, mpi_vli_type, 0, 1, MPI_COMM_WORLD, ERROR) 

        print*, image_id, 'finished sending polynomial' 
    end subroutine

    subroutine recv_polynomial(pah, from_node)
        use types_module
        use structure_module
        type(structure), intent(inout) :: pah
        integer, intent(in) :: from_node
        integer :: error
        integer :: mpi_status(MPI_STATUS_SIZE) 

        print*, image_id, 'getting polynomial' 

        call MPI_RECV(pah%order, 1, MPI_INTEGER8, from_node, 1, MPI_COMM_WORLD, mpi_status, ERROR)
        if ( .not. allocated(pah%polynomial) ) then
            allocate(pah%polynomial(pah%order+1))
        else if ( size(pah%polynomial) /= pah%order+1 ) then
            deallocate(pah%polynomial)
            allocate(pah%polynomial(pah%order+1))
        end if
        call MPI_RECV(pah%polynomial, pah%order+1, mpi_vli_type, from_node, 1, MPI_COMM_WORLD, mpi_status, ERROR)
        print*, image_id, 'finished getting polynomial' 
    end subroutine

end module
