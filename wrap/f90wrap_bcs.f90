! Module bcs defined in file bcs.f90

subroutine f90wrap_gram_schmidt(a, q, n0, n1, n2, n3)
    use bcs, only: gram_schmidt
    implicit none
    
    real(8), dimension(n0,n1) :: a
    real(8), dimension(n2,n3) :: q
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(a) :: n1 = shape(a,1)
    integer :: n2
    !f2py intent(hide), depend(q) :: n2 = shape(q,0)
    integer :: n3
    !f2py intent(hide), depend(q) :: n3 = shape(q,1)
    call gram_schmidt(A=a, Q=q)
end subroutine f90wrap_gram_schmidt

subroutine f90wrap_orthogonalize_to_set_list(randvec, veclist, n0, n1, n2)
    use bcs, only: orthogonalize_to_set_list
    implicit none
    
    real(8), dimension(n0) :: randvec
    real(8), dimension(n1,n2) :: veclist
    integer :: n0
    !f2py intent(hide), depend(randvec) :: n0 = shape(randvec,0)
    integer :: n1
    !f2py intent(hide), depend(veclist) :: n1 = shape(veclist,0)
    integer :: n2
    !f2py intent(hide), depend(veclist) :: n2 = shape(veclist,1)
    call orthogonalize_to_set_list(randvec=randvec, veclist=veclist)
end subroutine f90wrap_orthogonalize_to_set_list

subroutine f90wrap_reweight_penalty(j0, epsilon, ret_reweight_penalty, fxn)
    use bcs, only: reweight_penalty
    implicit none
    
    real(8) :: j0
    real(8) :: epsilon
    real(8), intent(out) :: ret_reweight_penalty
    character(6), optional :: fxn
    ret_reweight_penalty = reweight_penalty(j0=j0, epsilon=epsilon, fxn=fxn)
end subroutine f90wrap_reweight_penalty

subroutine f90wrap_do_wrapped(full_pi, y, sigma2, eta, js, error_bars, n0, n1, &
    n2, n3, n4)
    use bcs, only: do_wrapped
    implicit none
    
    real(8), dimension(n0,n1) :: full_pi
    real(8), dimension(n2) :: y
    real(8) :: sigma2
    real(8) :: eta
    real(8), dimension(n3) :: js
    real(8), dimension(n4) :: error_bars
    integer :: n0
    !f2py intent(hide), depend(full_pi) :: n0 = shape(full_pi,0)
    integer :: n1
    !f2py intent(hide), depend(full_pi) :: n1 = shape(full_pi,1)
    integer :: n2
    !f2py intent(hide), depend(y) :: n2 = shape(y,0)
    integer :: n3
    !f2py intent(hide), depend(js) :: n3 = shape(js,0)
    integer :: n4
    !f2py intent(hide), depend(error_bars) :: n4 = shape(error_bars,0)
    call do_wrapped(full_pi=full_pi, y=y, sigma2=sigma2, eta=eta, js=js, &
        error_bars=error_bars)
end subroutine f90wrap_do_wrapped

subroutine f90wrap_partition_holdout_set(nfits, nsets, nholdout, fitlist, &
    holdlist, n0, n1, n2, n3)
    use bcs, only: partition_holdout_set
    implicit none
    
    integer :: nfits
    integer :: nsets
    integer :: nholdout
    integer, dimension(n0,n1) :: fitlist
    integer, dimension(n2,n3) :: holdlist
    integer :: n0
    !f2py intent(hide), depend(fitlist) :: n0 = shape(fitlist,0)
    integer :: n1
    !f2py intent(hide), depend(fitlist) :: n1 = shape(fitlist,1)
    integer :: n2
    !f2py intent(hide), depend(holdlist) :: n2 = shape(holdlist,0)
    integer :: n3
    !f2py intent(hide), depend(holdlist) :: n3 = shape(holdlist,1)
    call partition_holdout_set(nfits=nfits, nsets=nsets, nholdout=nholdout, &
        fitlist=fitlist, holdlist=holdlist)
end subroutine f90wrap_partition_holdout_set

subroutine f90wrap_newunit(ret_newunit, unit)
    use bcs, only: newunit
    implicit none
    
    integer, intent(out) :: ret_newunit
    integer, optional :: unit
    ret_newunit = newunit(unit=unit)
end subroutine f90wrap_newunit

subroutine f90wrap_file_exists(ret_file_exists, filename)
    use bcs, only: file_exists
    implicit none
    
    logical, intent(out) :: ret_file_exists
    character(*) :: filename
    ret_file_exists = file_exists(filename=filename)
end subroutine f90wrap_file_exists

subroutine f90wrap_value_count(line, ret_value_count, length)
    use bcs, only: value_count
    implicit none
    
    character(1024) :: line
    integer, intent(out) :: ret_value_count
    integer :: length
    ret_value_count = value_count(line=line, length=length)
end subroutine f90wrap_value_count

subroutine f90wrap_linevalue_count(filename, n, commentchar, nlines, nvalues)
    use bcs, only: linevalue_count
    implicit none
    
    character(1024) :: filename
    integer :: n
    character(1) :: commentchar
    integer :: nlines
    integer :: nvalues
    call linevalue_count(filename=filename, n=n, commentchar=commentchar, &
        nlines=nlines, nvalues=nvalues)
end subroutine f90wrap_linevalue_count

subroutine f90wrap_write_results(js, fit_rms, fit_err, hold_rms_, hold_err_, &
    sigma2_, n0, n1, n2, n3, n4, n5, n6)
    use bcs, only: write_results
    implicit none
    
    real(8), dimension(n0,n1) :: js
    real(8), dimension(n2) :: fit_rms
    real(8), dimension(n3) :: fit_err
    real(8), optional, dimension(n4) :: hold_rms_
    real(8), optional, dimension(n5) :: hold_err_
    real(8), optional, dimension(n6) :: sigma2_
    integer :: n0
    !f2py intent(hide), depend(js) :: n0 = shape(js,0)
    integer :: n1
    !f2py intent(hide), depend(js) :: n1 = shape(js,1)
    integer :: n2
    !f2py intent(hide), depend(fit_rms) :: n2 = shape(fit_rms,0)
    integer :: n3
    !f2py intent(hide), depend(fit_err) :: n3 = shape(fit_err,0)
    integer :: n4
    !f2py intent(hide), depend(hold_rms_) :: n4 = shape(hold_rms_,0)
    integer :: n5
    !f2py intent(hide), depend(hold_err_) :: n5 = shape(hold_err_,0)
    integer :: n6
    !f2py intent(hide), depend(sigma2_) :: n6 = shape(sigma2_,0)
    call write_results(js=js, fit_rms=fit_rms, fit_err=fit_err, hold_rms_=hold_rms_, &
        hold_err_=hold_err_, sigma2_=sigma2_)
end subroutine f90wrap_write_results

! End of module bcs defined in file bcs.f90

