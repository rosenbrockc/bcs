  !!<summary>Wrapper module for calculating finding a solution vector to an
  !!underdetermined linear system by minimizing the ell_1 norm with Bayesian
  !!constraints. It is assumed that a collection of 'data sets' are available
  !!with each one having an 1) experimental measurement and 2) a vector of
  !!coefficients from evaluating the set within some basis that matches the
  !!constraints of compressive sensing. The collection of data sets' function
  !!evaluation vectors is called 'full_pi' and the measurements vector is called
  !!'y'.</summary>
module bcs
  use laplace

contains
  !!<summary>The purpose of this routine is to find a set of structures that are
  !!uncorrelated (in the sense discussed by Candes in his Compressive
  !!Sampling review).</summary>
  !!<comments author="GLWH" date="March 2013">Original routine by LJN March 2013;
  !!Updated to Vidvud's algorithm by GLWH March 2013.</comments>
  !!<comments>The basic algorithm is this:
  !! 1. Start with an empty set of selected data sets (a data set is a vector of
  !!    basis function evaluations for a single CS measurement).
  !! 2. Generate a random vector "pi" on the unit hypersphere.
  !! 3. Orthogonalize "pi" to the subspace spanned by the current set of data sets.
  !! 4. Normalize "pi".
  !! 5. Find the nearest data set to the orthonormalized "pi" using Gus' ArcCos distance.
  !! 6. Add this data set to the set of selected sets.
  !! 7. Construct an orthonormal basis in the space of selected data sets.
  !! 8. Go back to 2 until the desired size of selected set is reached.
  !!</comments>
  !!<parameter name="full_pi">The entire matrix of data sets (as defined above).</parameter>
  !!<parameter name="nrandsets">The number of random sets to select.</parameter>
  !!<parameter name="selected">The indices of the randomly selected data sets (rows in the
  !!@CREF[param.full_pi] matrix.</parameter>
  !!<parameter name="seed">The seed for the random number generator; set using the
  !!system clock if unspecified.</parameter>
  SUBROUTINE choose_n_random_sets(full_pi, nrandsets, selected, seed)
    real(dp), allocatable :: full_pi(:,:)
    integer, intent(in) :: nrandsets
    integer, pointer :: selected(:)
    integer, optional, allocatable :: seed(:)

    !!<local name="n, clock">The length of the random seed and system clock time.</local>
    !!<local name="nsets, nbasis">The number of data sets and basis function evaluations
    !!(rows and columns) in @CREF[param.full_pi].</local>
    !!<local name="randset">A randomly data set of basis function evaluations.</local>
    !!<local name="tempset">Variable used for normalizing the data set rows in
    !!@CREF[param.full_pi] for comparison with @CREF[randset].</local>
    !!<local name="closest">*Minimum* distance to the random vector for the current data set so
    !!far when sorting through data sets to find the one closest to the random vector.</local>
    !!<local name="distance">Distance to the random vector for the current data set
    !!when sorting through data sets.</local>
    !!<local name="keepidx">Index of the current data set in the @CREF[param.full_pi]
    !!that is closest.</local>
    !!<local name="inList">Used to hold the result from checking whether a data set is
    !!already in the solution vector of indices @CREF[param.selected].</local>
    integer n, clock
    integer nsets, nbasis
    real(dp), allocatable :: randset(:), tempset(:)
    real(dp) closest, distance ! For sorting the candidate structures
    integer i, j, keepIdx
    integer inList(1)

    if (.not. present(seed)) then
       call random_seed(size = n)
       allocate(seed(n))
       call system_clock(count = clock)
       seed = clock + 37 * (/ (i -1, i = 1, n)/)
    end if
    
    ! Seed the random number generator using the provided or calculated seed.
    call random_seed(put = seed)
    deallocate(seed)

    nsets = size(full_pi, 1)
    nbasis = size(full_pi, 2)

    allocate(randset(nbasis), tempset(nbasis))
    allocate(selected(nrandsets))
    nrandsets = 0

    ! Loop until the desired number of uncorrelated data sets are found.
    do i = 1, nrandsets
       ! Get a random vector of basis function evaluations.
       call random_number(randset)
       ! Rescale the vector so that each number is between -1 and 1
       randset = randset *2 - 1
       ! Before sorting through the list of data sets and determining
       ! which set this vector is "closest" to, orthogonalize it
       ! with respect to all the structures already chosen
       randset = randset/norm(randset)
       
       if (i > 1) then 
          call orthogonalize_to_set_list(randset, full_pi(selected(1:i-1),:))
       end if

       !We are hoping that the minimum distance to the random vector will be less than this
       !it isn't really a risky assumption.
       closest = 1000 
       keepIdx = -1
       do j = 1, nrandsets
          ! Normalize this data set's basis function evaluations' vector.
          tempset(:) = full_pi(j,:)/norm(full_pi(j,:))
          ! Compute the distance (arc length) between this data set's normalized
          ! evaluations' vector to the random vector drawn above. 
          distance = acos(dot_product(tempset,randset))
          ! Is data set j in the list already?
          inList = maxloc(selected, mask=(selected == j))
          ! If not and this distance is smaller than all previous, save this data set.
          if (distance < closest .and. inList(1) == 0) then
             closest = distance
             keepIdx = j
          end if
       end do
       ! Add the data set that was the closest to the random vector in our final list.
       selected(i) = keepIdx
    end do
  end subroutine choose_n_random_sets

  !!<summary>This is a "modified" Gram-Schmidt orthogonalization routine.</summary>
  !!<comments>Modified Gram-Schmidt routine retrieved January 2015 from
  !!http://terminus.sdsu.edu/SDSU/Math543_s2008/Lectures/07/lecture.pdf; while it
  !!isn't as good as householder, it is really simple.</comments>
  !!<parameter name="A">The matrix whose columns represent the vectors to orthogonalize.</parameter>
  !!<parameter name="Q">The orthogonal matrix from the QR decomposition by
  !!classic Gram-Schmidt.</parameter>
  subroutine gram_schmidt(A, Q)
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(out) :: Q(size(A, 1), size(A, 2))

    !!<local name="qii">The ith unit vector for Q.</local>
    !!<local name="ncol, nrow">The number of columns and rows in A.</local>
    real(dp) :: qi(size(A, 1))
    integer ncol
    integer i, j

    ncol = size(A, 2)
    Q = A
    
    do i = 1, ncol
       qi = Q(:, i)/norm(Q(:, i))
       do j = (i+1), ncol
          Q(:,j) = Q(:,j) - dot_product(qi, Q(:, j))*qi
       end do
    end do
  end subroutine gram_schmidt
  
  !!<summary>This routine takes a random vector on a unit hypersphere and
  !!orthogonalizes it with respect to an existing list of vectors
  !!(themselves not necessarily orthogonal). The output vector is normalized.</summary>
  !!<parameter name="randvec">The random vector to orthogonalize.</parameter>
  !!<parameter name="veclist">The list of existing vectors to orthogonalize against.</parameter>
  subroutine orthogonalize_to_set_list(randvec, veclist)
    real(dp), intent(inout) :: randvec(:)
    real(dp), intent(in)    :: veclist(:,:)

    !!<local name="ortholist">The list of orthogonalized data sets in *row* format.</local>
    !!<local name="tr_ortholist">The *transposed* list of orthogonalized data sets
    !!returned by Gram-Schmidt in *column* format.</local>
    real(dp) :: ortholist(size(veclist, 1), size(veclist, 2))
    real(dp) :: tr_ortholist(size(veclist,2), size(veclist, 1))
    integer :: i

    call gram_schmidt(transpose(veclist), tr_ortholist)
    ortholist = transpose(tr_ortholist)

    ! You want to do this projection in a loop (rather than "vectorized"),
    ! updating the vector as you go because it is numerically stable that
    ! way. (See discussions about stability of Gram-Schmidt...same issues
    ! here are important here.)
    do i = 1, size(veclist, 1)
       randvec = randvec - dot_product(randvec, ortholist(i,:))*ortholist(i,:)
    enddo
    randvec = randvec/norm(randvec)
  end subroutine orthogonalize_to_set_list

  !!<summary>Returns the value of the weight that should be used for the
  !!weighting matrix in the bcs reweighting scheme.</summary>
  !!<comments author="Conrad W. Rosenbrock" date="11 Jun 2013">
  !!This method is based on Candes "Enhancing sparsity by reweighted ell 1 minimization",
  !!eq 9 and discussion following additional penalty functions are based on the
  !!transfer functions used in machine learning -> neural networks because they
  !!share similarities in the requirements of being bounded, differentiable, etc.    
  !!</comments>
  !!<parameter name="j0">The current value of the j-coefficient for the
  !!corresponding basis</parameter>
  !!<parameter name="epsilon">!A fudge factor for numerical stability.</parameter>
  !!<parameter name="fxn">The name of the penalty function to use. Possible values:
  !!logsum, logsig, arctan, quarti, hexics, octics.</parameter>
  real(dp) function reweight_penalty(j0, epsilon, fxn)
    real(dp), intent(in) :: j0 
    real(dp), intent(in) :: epsilon 
    character(len=6), intent(in), optional :: fxn 
 
    select case (fxn)
    case ('logsum')
       reweight_penalty = j0 + epsilon
    case ('logsig')
       reweight_penalty = 4*exp(-epsilon + j0)*(1 + exp(epsilon - j0))**2
    case ('arctan')
       reweight_penalty = j0**2 + epsilon**2
    case ('quarti')
       reweight_penalty = 1+(j0 + 0.7*epsilon)**4
    case ('hexics')
       reweight_penalty = 1+(j0 + 0.57*epsilon)**6
    case ('octics')
       reweight_penalty = 1+(j0 + 0.44*epsilon)**8
    case default
       reweight_penalty = j0 + epsilon
    end select
  end function reweight_penalty
  
  !!<summary>Finds the optimal value for sigma^2 and checks its validity to detect
  !!data that isn't well suited to the BCS routine.</summary>
  !!<parameter name="y">The experimental measurements for each data set.</parameter>
  !!<parameter name="full_pi">The full matrix of basis function evaluations for
  !!each data set for which we have measurements.</parameter>
  !!<parameter name="nfit">The number of data sets allowed to be used in the fit.</parameter>
  real(dp) function get_sigma2(full_pi, y, nfit)
    real(dp), pointer, intent(in) :: full_pi(:,:), y(:)
    integer, intent(in) :: nfit

    !This estimate for sigma^2 is just (<y_ii^2> - <y_ii>^2)/100.
    get_sigma2 = (sum( (/( y(i)**2, i=1,nfit )/) )/nfit - (sum( (/( y(i), i=1,nfit )/) )/nfit)**2)/100
    if (get_sigma2 > 100) then
       print *, "There is a problem with your data (in function get_sigma2). "
       print *, "sigma^2 =  ", get_sigma2
       print *, "This value is way too big, there may be some very poor data points in the fitting set."
       stop
    end if
  end function get_sigma2
  
  !!<summary>This routine wraps the Bayesian CS routine with an outside loop. Only
  !!regular ell_1 minimization is done (as opposed to p<1 reweighted minimization).</summary>
  !!<parameter name="y">The experimental measurements for each data set.</parameter>
  !!<parameter name="full_pi">The full matrix of basis function evaluations for
  !!each data set for which we have measurements.</parameter>
  !!<parameter name="sigma2">Initial noise variance (default : std(t)^2/1e2).</parameter>
  !!<parameter name="eta">Threshold for stopping the iteration algorithm.</parameter>
  !!<parameter name="js">The solution vector for the reweighted BCS.</parameter>
  !!<parameter name="error_bars">Error bars for the solution vector @CREF[param.js]; taking
  !!from the diagonal entries of the covariance matrix.</parameter>
  subroutine do_normal(full_pi, y, sigma2, eta, js, error_bars)
    real(dp), pointer, intent(in) :: full_pi(:,:), y(:)
    real(dp), intent(in) :: sigma2, eta, jcutoff
    character(len=6), intent(in) :: penaltyfxn
    real(dp), intent(out) :: js(size(full_pi, 2)), error_bars(size(full_pi, 2))

    !!<local name="iterator">The laplace iterator to find the solution vector.</local>
    !!<local name="returnsigma2">The value of sigma2 calculated by the laplace
    !!iterator if one was not provided. We don't actually output this anywhere.</local>
    real(dp) returnsigma2
    type(laplace_iterator), pointer :: iterator
    
    allocate(iterator)
    call iterator%initialize(full_pi, y, sigma2, eta)
    call iterator%iterate(js, error_bars, returnsigma2)
    deallocate(this%iterator)
  end subroutine do_normal
  
  !!<summary>This routine wraps the Bayesian CS routine with an outside loop that
  ! applies the p<1 norm reweighting scheme.</summary>
  !!<comments>See: E. J. Candès, M. Wakin and S. Boyd. Enhancing sparsity by
  !! reweighted \ell_1 minimization. J. Fourier Anal. Appl., 14 877-905.</comments>
  !!<parameter name="y">The experimental measurements for each data set.</parameter>
  !!<parameter name="full_pi">The full matrix of basis function evaluations for
  !!each data set for which we have measurements.</parameter>
  !!<parameter name="sigma2">Initial noise variance (default : std(t)^2/1e2).</parameter>
  !!<parameter name="eta">Threshold for stopping the iteration algorithm.</parameter>
  !!<parameter name="jcutoff">The minimum value a J coefficient has to have
  !!in order to be kept between reweighting iterations.</parameter>
  !!<parameter name="penaltyfxn">The name of the penalty function to use. Possible values:
  !!logsum, logsig, arctan, quarti, hexics, octics.</parameter>
  !!<parameter name="js">The solution vector for the reweighted BCS.</parameter>
  !!<parameter name="error_bars">Error bars for the solution vector @CREF[param.js]; taking
  !!from the diagonal entries of the covariance matrix.</parameter>
  subroutine do_reweighted(full_pi, y, sigma2, eta, jcutoff, penaltyfxn, js, error_bars)
    real(dp), pointer, intent(in) :: full_pi(:,:), y(:)
    real(dp), intent(in) :: sigma2, eta, jcutoff
    character(len=6), intent(in) :: penaltyfxn
    real(dp), intent(out) :: js(size(full_pi, 2)), error_bars(size(full_pi, 2))

    !!<local name="returnsigma2">The value of sigma2 calculated by the laplace
    !!iterator if one was not provided. We don't actually output this anywhere.</local>
    !!<local name="w_pi">The *weighted* @CREF[param.full_pi] matrix.</local>
    !!<local name="weight_matrix">The matrix of adjusted/re-weighted J values.</local>
    !!<local name="nsets, nbasis">The number of data sets (rows) and basis functions
    !!(columns) in the @CREF[param.full_pi] matrix.</local>
    !!<local name="copyJs">A copy of the current js solution vector for manipulation.</local>
    !!<local name="ell0, prevell0">The iteration continues until the ell_0 norm of
    !!the solution vector stops changing significantly. These variables keep track of
    !!the previous and current ell_0 norm.</local>
    !!<local name="iterator">The laplace iterator to find the solution vector.</local>
    !!<local name="i_0">The maximum number of non-negligible coefficient values
    !!in order for the compressive sensing mathematics to be applicable.</local>
    !!<local name="maxidx">The index of the maximum value in the copy of the J
    !!solution vector. Helps determine the value of the reweighting epsilon parameter.</local>
    !!<local name="jmax">The value of the largest j value outside of the range
    !!specifed by 'i_0'.</local>
    !!<local name="epsilon">The reweighting parameter that controls how aggresively the
    !!weight matrix is altered.</local>
    real(dp) returnsigma2
    real(dp), allocatable :: w_pi(:,:)
    real(dp), allocatable :: weight_matrix(:,:)
    integer :: nsets, nbasis
    real(dp), allocatable :: copyJs(size(full_pi,2))
    integer :: ell0, prevell0
    type(laplace_iterator), pointer :: iterator
    integer  :: i_0, maxidx(1)
    real(dp) :: jmax(1), epsilon
    integer :: i

    !We start off with the identity matrix for the weight matrix and then
    !update it after each iteration.
    allocate(weight_matrix(nbasis,nbasis), w_pi(nsets,nbasis))
    allocate(iterator)
    weight_matrix = 0
    do i = 1, this%nCorrs
       weight_matrix(i,i) = 1
    end do

    !Here we make a rough estimate of the number of J coefficients that will
    !have non-negligible values (i.e. > 1e-3). We hope since the solution is
    !supposed to be sparse that this condition will hold.
    i_0 = nsets/(4*log(real(nbasis)/real(nsets)))
    prevell0 = 3000
    
    do while (.true.)
       w_pi = matmul(full_pi, weight_matrix)
       call iterator%initialize(w_pi, y, sigma2, eta)
       call iterator%iterate(js, error_bars, returnsigma2)

       js(:) = matmul(weight_matrix, js)
       ell0 = count(abs(js) > jcutoff)
       if (abs(ell0 - prevell0) < 5) then
          exit
       else
          prevell0 = ell0
       end if

       !Set the largest i_0 j-values to zero so we can get an estimate of
       !the size of the smallest coefficients. They affect the reweighting
       !matrix and penalty functions.
       copyJs = abs(js)
       do i=1, i_0-1
          maxidx = maxloc(copyJs)
          copyJs(maxidx(1)) = 0
       end do
       jmax = maxval(copyJs)
       if (jmax(1) > 10e-3) then
          epsilon = jmax(1)
       else
          epsilon = 10e-3
       end if

       !Adjust the values of the diagonals on the weighting matrix using the
       !specified penalty function.
       do i=1, nbasis
          !This is the inverse W matrix.
          weight_matrix(i,i) = reweight_penalty(js(i), epsilon, penaltyfxn)
       end do
     
       call iterator%reset()
    end do

    !Perform cleanup on the iterator object. This will call its finalizer once
    !gfortran implements them. For now we rely on the built-in deallocation.
    deallocate(iterator)    
  end subroutine do_reweighted

  !!<summary>Performs BCS for the given sensing matrix and measurement vector.</summary>
  !!<parameter name="y">The experimental measurements for each data set.</parameter>
  !!<parameter name="full_pi">The full matrix of basis function evaluations for
  !!each data set for which we have measurements.</parameter>
  !!<parameter name="nfits">The number of random fits to perform with these args.</parameter>
  !!<parameter name="js">The solution vector for the reweighted BCS. There is one row
  !!in this matrix for each fit performed.</parameter>
  !!<parameter name="error_bars">Error bars for the solution vector @CREF[param.js]; taking
  !!from the diagonal entries of the covariance matrix. There is one row in this matrix
  !!for each fit performed.</parameter>
  !!<parameter name="sigma2">Initial noise variance (default : std(t)^2/1e2).</parameter>
  !!<parameter name="eta">Threshold for stopping the iteration algorithm.</parameter>
  !!<parameter name="jcutoff">The minimum value a J coefficient has to have
  !!in order to be kept between reweighting iterations.</parameter>
  !!<parameter name="penaltyfxn">The name of the penalty function to use. Possible values:
  !!logsum, logsig, arctan, quarti, hexics, octics.</parameter>
  subroutine do_bcs(full_pi, y, nfits, js, error_bars, sigma2, eta, jcutoff, penaltyfxn, reweight)
    real(dp), pointer, intent(in) :: full_pi(:,:), y(:)
    integer :: nfits
    real(dp), intent(out) :: js(nfits, size(full_pi, 2)), error_bars(nfits, size(full_pi, 2))
    real(dp), intent(in) :: sigma2, eta, jcutoff
    character(len=6), intent(in) :: penaltyfxn

  end subroutine do_bcs
end module bcs
