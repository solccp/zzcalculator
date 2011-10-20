
module mpi_global
    use mpi
    use accuracy_m
    use structure_m
    use tree_m
    implicit none
    save
    integer :: image_count
    integer :: image_id 
    integer :: mpi_vli_type
 
contains
    subroutine mpi_global_init()
        integer :: error
    
        call MPI_Comm_rank ( MPI_COMM_WORLD, image_id, error )
        call MPI_Comm_size ( MPI_COMM_WORLD, image_count, error )

        call MPI_TYPE_CONTIGUOUS(block_size+1, MPI_INTEGER8, mpi_vli_type, error)
        call MPI_TYPE_COMMIT(mpi_vli_type, error)

    end subroutine

    subroutine mpi_global_finalize()
        integer :: error
    
        call MPI_TYPE_FREE(mpi_vli_type, error)
    end subroutine


    subroutine send_polynomial(pah_array, local_index, index_size)
        type(structure_ptr), dimension(:), intent(in) :: pah_array
        integer, intent(in) :: index_size
        integer, dimension(:), intent(in) :: local_index

        type(structure), pointer :: pah

        integer :: i

        integer :: error
        integer :: pos

        integer :: total_poly_size

        character, dimension(:), allocatable :: buff
        integer :: buff_size

!        print*, image_id, 'sending polynomial' 

        !start packing the buff

        pos = 0
        total_poly_size = 0          
        do i = 1, index_size
            pah => pah_array(local_index(i))%ptr
            total_poly_size = total_poly_size + 8*((pah%order+1)*(block_size+1)) + 4
        end do

        buff_size = total_poly_size + index_size*4 + 4 
        allocate(buff(buff_size))

        

        call MPI_PACK(index_size, 1, MPI_INTEGER4, buff, buff_size, pos, MPI_COMM_WORLD, ERROR)
        do i = 1, index_size
            pah => pah_array(local_index(i))%ptr
            call MPI_PACK(local_index(i), 1, MPI_INTEGER4, buff, buff_size, pos, MPI_COMM_WORLD, ERROR)
            call MPI_PACK(pah%order, 1, MPI_INTEGER4, buff, buff_size, pos, MPI_COMM_WORLD, ERROR)
            call MPI_PACK(pah%polynomial, pah%order+1, mpi_vli_type, buff, buff_size, pos, MPI_COMM_WORLD, ERROR) 
        end do

        call MPI_SEND(buff_size, 1, MPI_INTEGER4, 0, 200, MPI_COMM_WORLD, ERROR)
        call MPI_SEND(buff, pos, MPI_CHARACTER, 0, 201, MPI_COMM_WORLD, ERROR)

!        print*, image_id, 'finished sending polynomial' 
        deallocate(buff)
    end subroutine

    subroutine recv_polynomial(pah_array, from_node)
        type(structure_ptr), dimension(:), intent(inout) :: pah_array
        integer, intent(in) :: from_node
!        integer, intent(out) :: mpi_request
        
        type(structure), pointer :: pah

        integer :: error
        integer :: mpi_status(MPI_STATUS_SIZE) 

        integer :: buff_size
        character, dimension(:), allocatable :: buff
        integer :: pos
        integer :: local_index, index_size
        integer :: i

!        print*, image_id, 'getting polynomial from', from_node


        call MPI_RECV(buff_size, 1, MPI_INTEGER4, from_node, 200, MPI_COMM_WORLD, mpi_status, ERROR)
        allocate(buff(buff_size))
        call MPI_RECV(buff, buff_size, MPI_CHARACTER, from_node, 201, MPI_COMM_WORLD, mpi_status, ERROR)
        pos = 0
        call MPI_UNPACK(buff, buff_size, pos, index_size, 1, MPI_INTEGER4, MPI_COMM_WORLD, ERROR)

        do i = 1, index_size
            call MPI_UNPACK(buff, buff_size, pos, local_index, 1, MPI_INTEGER4, MPI_COMM_WORLD, ERROR)
            pah => pah_array(local_index)%ptr
            call MPI_UNPACK(buff, buff_size, pos, pah%order, 1, MPI_INTEGER4, MPI_COMM_WORLD, ERROR)
            if ( .not. allocated(pah%polynomial) .or. (size(pah%polynomial) /= (pah%order+1)) ) then
                if ( allocated(pah%polynomial) ) then
                    deallocate(pah%polynomial)
                end if
                allocate(pah%polynomial(pah%order+1))
            end if
                call MPI_UNPACK(buff, buff_size, pos, pah%polynomial, pah%order+1, mpi_vli_type, MPI_COMM_WORLD, ERROR) 
        end do

!        print*, image_id, 'finished getting polynomial', from_node
        deallocate(buff)
    end subroutine


end module
