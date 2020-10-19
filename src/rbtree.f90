module rbtree_m
    use iso_varying_string
    use structure_m
    use rbtree_key_m
    implicit none  

    type rb_treenode
        type(rb_treenode), pointer :: left => null()
        type(rb_treenode), pointer :: right => null()
        type(rb_treenode), pointer :: parent => null()
        logical :: red = .false.
        type(nodedata)  :: keyvalue
    end type rb_treenode
    
    type rb_treenode_ptr
        type(rb_treenode), pointer :: ptr => null()
    end type

    type RedBlackTree
        type(rb_treenode), pointer :: root => null()
        integer :: nitems = 0
    end type RedBlackTree

    type(rb_treenode), target, save :: nil
    public :: insert_node, RedBlackTree, nodedata, composite_key, getMin, getMax, isNil
    private :: nil

contains

function isNil(pos) result(res)
    type(rb_treenode_ptr) :: pos
    logical :: res
    res = isLeaf(pos%ptr)
end function

subroutine getMin_node(node, pos)
    type(rb_treenode), target :: node
    type(rb_treenode_ptr) :: pos
    type(rb_treenode), pointer :: x 

    x => node
    if ( isLeaf(x) ) then
        pos%ptr => null()
        return
    end if
    do while( .not. isLeaf(x%left) )
        x => x%left
    end do
    pos%ptr => x
end subroutine
subroutine getMin(tree, pos)
    type(RedBlackTree) :: tree
    type(rb_treenode_ptr) :: pos
    call getMin_node(tree%root, pos)
end subroutine

subroutine getMax_node(node, pos)
    type(rb_treenode), target :: node
    type(rb_treenode_ptr) :: pos
    type(rb_treenode), pointer :: x 

    x => node
    if ( isLeaf(x) ) then
        pos%ptr => nil
        return
    end if
    do while( .not. isLeaf(x%right) )
        x => x%right
    end do
    pos%ptr => x
end subroutine
subroutine getMax(tree, pos)
    type(RedBlackTree) :: tree
    type(rb_treenode_ptr) :: pos
    call getMax_node(tree%root, pos)
end subroutine

subroutine find(tree, ckey, pos)
    type(RedBlackTree) :: tree
    type(composite_key) :: ckey
    type(rb_treenode_ptr) :: pos
    type(rb_treenode), pointer :: x 

    integer :: comparator

    x => tree%root
    pos%ptr => nil

    do while ( .not. isLeaf(x) )
        if ( comparator(ckey, x%keyvalue%ckey) == 0 ) then
            pos%ptr => x
            exit
        else if ( comparator(ckey, x%keyvalue%ckey) == -1 ) then
            x => x%left
        else
            x => x%right
        end if
    end do
end subroutine
subroutine init(tree)
    type (RedBlackTree) :: tree
    tree%root => nil
end subroutine init

recursive subroutine print_tree(node, level)
    type(rb_treenode), pointer :: node
    integer :: level 
    if ( .not. isLeaf(node) ) then
        if ( node%red ) then
            write (*,'(i0,1x,a,1x,i0,1x,i0)') node%keyvalue%ckey%key1, 'RED', &
            node%left%keyvalue%ckey%key1, node%right%keyvalue%ckey%key1
        else
            write (*,'(i0,1x,a,1x,i0,1x,i0)') node%keyvalue%ckey%key1, 'BLK', &
            node%left%keyvalue%ckey%key1, node%right%keyvalue%ckey%key1
        end if
        if (.not. isLeaf(node%left)) then
            if ( .not. associated(node%left%parent, node) ) then
                print *, 'error in left child', node%keyvalue%ckey%key1
                stop
            end if
            call print_tree(node%left, level+1)
        end if
        if (.not. isLeaf(node%right)) then
            if ( .not. associated(node%right%parent, node) ) then
                print *, 'error in right child', node%keyvalue%ckey%key1
                stop
            end if
            call print_tree(node%right, level+1)
        end if
    end if
end subroutine
subroutine erase_node(tree, node)
    type (RedBlackTree), target :: tree
    type(rb_treenode), pointer :: node, y, x 
    type(rb_treenode_ptr) :: tmp

    if ( isLeaf(node) ) then
        return
    end if

    y => node
    if ( .not. isLeaf(y%left) .and. .not. isLeaf(y%right) ) then
        call getMin_node(node%right, tmp)
        y => tmp%ptr
    end if

    if ( .not. isLeaf(y%left) ) then
        x => y%left
    else
        x => y%right
    end if

    x%parent => y%parent
    if ( isLeaf(x%parent) ) then
        tree%root => x
    else if ( isLeft(y%parent, y) ) then
        y%parent%left => x
    else
        y%parent%right => x
    end if

    if ( .not. associated(y, node) ) then
        node%keyvalue = y%keyvalue
    end if

    if ( .not. y%red ) then
        call deleteFixup(tree, x)
    end if

    y%left => nil
    y%right => nil
    y%parent => nil
    y%keyvalue%value => null()
    nil%parent => tree%root
    
    call finalize_node(y)
    tree%nitems = tree%nitems - 1

end subroutine
subroutine deleteFixup(tree, node)
    type (RedBlackTree), target :: tree
    type(rb_treenode), pointer :: node, x, w, xp
    type(rb_treenode_ptr) :: tmp

    x => node
    do while ( .not. associated(x, tree%root) .and. .not. x%red ) 
        xp => x%parent
        if ( isLeft(xp, x) ) then
            w => xp%right
            if ( w%red ) then
                w%red = .false.
                xp%red = .true.
                call rotate_left(tree, xp)
                w => xp%right
            end if
            if ( isLeaf(w) ) then
                exit
            end if
            if ( .not. w%left%red .and. .not. w%right%red ) then
                w%red = .true.
                x => xp
            else
                if ( .not. w%right%red ) then
                    w%left%red = .false.
                    w%red = .true.
                    call rotate_right(tree, w)
                    w => w%parent
                end if
                w%red = w%parent%red
                w%parent%red = .false.
                w%right%red = .false.
                call rotate_left(tree, w%parent)
                x => tree%root
            end if
        else
            w => xp%left 
            if ( w%red ) then
                w%red = .false.
                xp%red = .true.
                call rotate_right(tree, xp)
                w => xp%left
            end if
            if ( isLeaf(w) ) then
                exit
            end if
            if ( .not. w%left%red .and. .not. w%right%red ) then
                w%red = .true.
                x => xp
            else
                if ( .not. w%left%red ) then
                    w%right%red = .false.
                    w%red = .true.
                    call rotate_left(tree, w)
                    w => w%parent
                end if
                w%red = w%parent%red
                w%parent%red = .false.
                w%right%red = .false.
                call rotate_right(tree, w%parent)
                x => tree%root
            end if
        end if
    end do
    x%red = .false.
end subroutine


subroutine insert_node(tree, keyvalue, inserted)
    type (RedBlackTree), target :: tree
    type(nodedata), intent(inout) :: keyvalue
    logical, intent(out) :: inserted
    type(rb_treenode), pointer :: node, p, leaf
    integer :: comparator

    if (.not. associated(tree%root)) return

    p => tree%root
    leaf => nil
    do while ( .not. isLeaf(p) )
        leaf => p       
        if ( comparator(keyvalue%ckey, p%keyvalue%ckey) == 0 ) then
            keyvalue = p%keyvalue
            inserted = .false.
            return
        else if ( comparator(keyvalue%ckey, p%keyvalue%ckey) < 0 ) then
            p => p%left
        else
            p => p%right
        end if
    end do

    allocate(node)
    node%keyvalue = keyvalue
    node%left => nil
    node%right => nil
    node%red = .true.

    inserted = .true.
    tree%nitems = tree%nitems + 1

    node%parent => leaf
    if ( isLeaf(leaf) ) then
        tree%root => node
    else if ( comparator(node%keyvalue%ckey, leaf%keyvalue%ckey) < 0 ) then
        leaf%left => node
    else
        leaf%right => node
    end if
    call insertFixup(tree, node)
end subroutine


subroutine insertFixup(tree, node)
    type (RedBlackTree), target :: tree
    type(rb_treenode), pointer :: node, p, u
    p => node
    do while ( p%parent%red ) 
        if ( isLeft(p%parent%parent, p%parent) ) then
            u => p%parent%parent%right
            if ( u%red ) then
                p%parent%red = .false.
                u%red = .false.
                p%parent%parent%red = .true.
                p => p%parent%parent
            else
                if ( isRight(p%parent, p) ) then
                    p => p%parent
                    call rotate_left(tree, p)
                end if
                p%parent%red = .false.
                p%parent%parent%red = .true.
                call rotate_right(tree, p%parent%parent)
            end if
        else
            u => p%parent%parent%left
            if ( u%red ) then
                p%parent%red = .false.
                u%red = .false.
                p%parent%parent%red = .true.
                p => p%parent%parent
            else
                if ( isLeft(p%parent, p) ) then
                    p => p%parent
                    call rotate_right(tree, p)
                end if
                p%parent%red = .false.
                p%parent%parent%red = .true.
                call rotate_left(tree, p%parent%parent)
            end if
        end if
    end do
    tree%root%red = .false.
end subroutine

subroutine rotate_left(tree, node)
    type(RedBlackTree), intent(inout) :: tree
    type(rb_treenode), target, intent(inout) :: node
    type(rb_treenode), pointer :: x 
    type(rb_treenode), pointer :: y 

    x => node
    y => x % right
    x % right => y % left

    if ( .not. isLeaf(y % left) ) then
       y % left % parent => x
    end if

    y % parent => x % parent

    if ( isLeaf(x % parent) ) then
       tree % root => y
    else
        if ( isLeft( x%parent, x) ) then 
            x % parent % left => y
        else
            x % parent % right => y
        end if
    end if
    y % left => x
    x % parent => y
end subroutine

subroutine rotate_right(tree, node)
    type(RedBlackTree), intent(inout) :: tree
    type(rb_treenode), target, intent(inout) :: node
    type(rb_treenode), pointer :: x 
    type(rb_treenode), pointer :: y 

    x => node
    y => x % left
    x % left => y % right

    if (.not. isLeaf(y % right)) then
        y % right % parent => x
    end if

    y % parent => x % parent

    if ( isLeaf(x % parent) ) then
       tree % root => y
    else
        if ( isRight( x%parent, x) ) then 
            x % parent % right => y
        else
            x % parent % left => y
        end if
    end if
    y % right => x
    x % parent => y
end subroutine

!========================================
! Check if node y is the left child of x
!========================================
pure function isLeft(x, y) result(res)
    type (rb_treenode), pointer :: x, y
    logical :: res
    res = associated(x % left, y)
end function

pure function isRight(x, y) result(res)
    type (rb_treenode), pointer :: x, y
    logical :: res
    res = associated(x % right, y)
end function

pure function isLeaf(x) result(res)
    type (rb_treenode), pointer :: x
    logical :: res
    res = associated(x, nil)
end function

recursive subroutine finalize_node(node)
    type (rb_treenode), pointer :: node
    if ( .not. isLeaf(node) ) then
        call finalize_node(node%left)
        call finalize_node(node%right)
        deallocate(node)
    end if
end subroutine

subroutine finalize(tree)
    type(RedBlackTree), intent(inout) :: tree
    call finalize_node(tree%root) 
end subroutine


end module

function comparator(lhs, rhs) result(res)
    use rbtree_m
    type(composite_key) :: lhs, rhs
    integer :: res
    if ( lhs%key1 < rhs%key1 ) then
        res = -1
    else if ( lhs%key1 == rhs%key1 ) then
        if ( lhs%key2 < rhs%key2 ) then
            res = -1
        else if ( lhs%key2 == rhs%key2 ) then
            res = 0
        else 
            res = 1
        end if
    else
        res = 1
    end if
end function

