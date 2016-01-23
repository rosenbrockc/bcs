! Module laplace defined in file laplace.f90

subroutine f90wrap_laplace_iterator_initialise(this)
    use laplace, only: laplace_iterator
    implicit none
    
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_laplace_iterator_initialise

subroutine f90wrap_laplace_iterator_finalise(this)
    use laplace, only: laplace_iterator
    implicit none
    
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_laplace_iterator_finalise

subroutine f90wrap_norm(ret_norm, vector, n0)
    use laplace, only: norm
    implicit none
    
    real(8), intent(out) :: ret_norm
    real(8), dimension(n0) :: vector
    integer :: n0
    !f2py intent(hide), depend(vector) :: n0 = shape(vector,0)
    ret_norm = norm(vector=vector)
end subroutine f90wrap_norm

! End of module laplace defined in file laplace.f90

