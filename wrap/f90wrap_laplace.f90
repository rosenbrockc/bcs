! Module laplace defined in file ../bcs/laplace.f90

subroutine f90wrap_laplace_iterator__array__s(this, nd, dtype, dshape, dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%s)) then
        dshape(1:1) = shape(this_ptr%p%s)
        dloc = loc(this_ptr%p%s)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__s

subroutine f90wrap_laplace_iterator__array__q(this, nd, dtype, dshape, dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%q)) then
        dshape(1:1) = shape(this_ptr%p%q)
        dloc = loc(this_ptr%p%q)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__q

subroutine f90wrap_laplace_iterator__array__ls(this, nd, dtype, dshape, dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ls)) then
        dshape(1:1) = shape(this_ptr%p%ls)
        dloc = loc(this_ptr%p%ls)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__ls

subroutine f90wrap_laplace_iterator__array__lq(this, nd, dtype, dshape, dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%lq)) then
        dshape(1:1) = shape(this_ptr%p%lq)
        dloc = loc(this_ptr%p%lq)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__lq

subroutine f90wrap_laplace_iterator__get__ninputs(this, ninputs)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: ninputs
    
    this_ptr = transfer(this, this_ptr)
    ninputs = this_ptr%p%ninputs
end subroutine f90wrap_laplace_iterator__get__ninputs

subroutine f90wrap_laplace_iterator__set__ninputs(this, ninputs)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(in) :: ninputs
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ninputs = ninputs
end subroutine f90wrap_laplace_iterator__set__ninputs

subroutine f90wrap_laplace_iterator__get__ncorrs(this, ncorrs)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: ncorrs
    
    this_ptr = transfer(this, this_ptr)
    ncorrs = this_ptr%p%ncorrs
end subroutine f90wrap_laplace_iterator__get__ncorrs

subroutine f90wrap_laplace_iterator__set__ncorrs(this, ncorrs)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(in) :: ncorrs
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ncorrs = ncorrs
end subroutine f90wrap_laplace_iterator__set__ncorrs

subroutine f90wrap_laplace_iterator__get__nused(this, nused)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nused
    
    this_ptr = transfer(this, this_ptr)
    nused = this_ptr%p%nused
end subroutine f90wrap_laplace_iterator__get__nused

subroutine f90wrap_laplace_iterator__set__nused(this, nused)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(in) :: nused
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nused = nused
end subroutine f90wrap_laplace_iterator__set__nused

subroutine f90wrap_laplace_iterator__get__sigma2(this, sigma2)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    real(8), intent(out) :: sigma2
    
    this_ptr = transfer(this, this_ptr)
    sigma2 = this_ptr%p%sigma2
end subroutine f90wrap_laplace_iterator__get__sigma2

subroutine f90wrap_laplace_iterator__set__sigma2(this, sigma2)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    real(8), intent(in) :: sigma2
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%sigma2 = sigma2
end subroutine f90wrap_laplace_iterator__set__sigma2

subroutine f90wrap_laplace_iterator__get__eta(this, eta)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    real(8), intent(out) :: eta
    
    this_ptr = transfer(this, this_ptr)
    eta = this_ptr%p%eta
end subroutine f90wrap_laplace_iterator__get__eta

subroutine f90wrap_laplace_iterator__set__eta(this, eta)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    real(8), intent(in) :: eta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eta = eta
end subroutine f90wrap_laplace_iterator__set__eta

subroutine f90wrap_laplace_iterator__get__lambda(this, lambda)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    real(8), intent(out) :: lambda
    
    this_ptr = transfer(this, this_ptr)
    lambda = this_ptr%p%lambda
end subroutine f90wrap_laplace_iterator__get__lambda

subroutine f90wrap_laplace_iterator__set__lambda(this, lambda)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    real(8), intent(in) :: lambda
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lambda = lambda
end subroutine f90wrap_laplace_iterator__set__lambda

subroutine f90wrap_laplace_iterator__array__alpha(this, nd, dtype, dshape, dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%alpha)) then
        dshape(1:1) = shape(this_ptr%p%alpha)
        dloc = loc(this_ptr%p%alpha)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__alpha

subroutine f90wrap_laplace_iterator__array__subphi(this, nd, dtype, dshape, &
    dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%subphi)) then
        dshape(1:2) = shape(this_ptr%p%subphi)
        dloc = loc(this_ptr%p%subphi)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__subphi

subroutine f90wrap_laplace_iterator__array__indices(this, nd, dtype, dshape, &
    dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%indices)) then
        dshape(1:1) = shape(this_ptr%p%indices)
        dloc = loc(this_ptr%p%indices)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__indices

subroutine f90wrap_laplace_iterator__array__selected(this, nd, dtype, dshape, &
    dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%selected)) then
        dshape(1:1) = shape(this_ptr%p%selected)
        dloc = loc(this_ptr%p%selected)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__selected

subroutine f90wrap_laplace_iterator__array__sigma(this, nd, dtype, dshape, dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%sigma)) then
        dshape(1:2) = shape(this_ptr%p%sigma)
        dloc = loc(this_ptr%p%sigma)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__sigma

subroutine f90wrap_laplace_iterator__array__mu(this, nd, dtype, dshape, dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%mu)) then
        dshape(1:1) = shape(this_ptr%p%mu)
        dloc = loc(this_ptr%p%mu)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__mu

subroutine f90wrap_laplace_iterator__array__theta(this, nd, dtype, dshape, dloc)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in) :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%theta)) then
        dshape(1:1) = shape(this_ptr%p%theta)
        dloc = loc(this_ptr%p%theta)
    else
        dloc = 0
    end if
end subroutine f90wrap_laplace_iterator__array__theta

subroutine f90wrap_laplace_iterator__get__initialized(this, initialized)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    logical, intent(out) :: initialized
    
    this_ptr = transfer(this, this_ptr)
    initialized = this_ptr%p%initialized
end subroutine f90wrap_laplace_iterator__get__initialized

subroutine f90wrap_laplace_iterator__set__initialized(this, initialized)
    use laplace, only: laplace_iterator
    implicit none
    type laplace_iterator_ptr_type
        type(laplace_iterator), pointer :: p => NULL()
    end type laplace_iterator_ptr_type
    integer, intent(in)   :: this(2)
    type(laplace_iterator_ptr_type) :: this_ptr
    logical, intent(in) :: initialized
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%initialized = initialized
end subroutine f90wrap_laplace_iterator__set__initialized

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
    real(8), intent(in), dimension(n0) :: vector
    integer :: n0
    !f2py intent(hide), depend(vector) :: n0 = shape(vector,0)
    ret_norm = norm(vector=vector)
end subroutine f90wrap_norm

! End of module laplace defined in file ../bcs/laplace.f90

