module rbtree_key_m
    use iso_varying_string
    use structure_m
    type composite_key
        integer :: key1
        type(varying_string)  :: key2
    end type

    !===========================
    ! define your own datatype
    !===========================
    type nodedata
        type(composite_key) :: ckey
        type(pah_tree_node), pointer :: value => null()
    end type
end module
