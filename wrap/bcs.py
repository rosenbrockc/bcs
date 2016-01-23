import _bcs
import f90wrap.runtime
import logging

class Bcs(f90wrap.runtime.FortranModule):
    """
    Module bcs
    
    
    Defined at ../bcs/bcs.f90 lines 9-1049
    
    """
    @staticmethod
    def gram_schmidt(a, q):
        """
        gram_schmidt(a, q)
        
        
        Defined at ../bcs/bcs.f90 lines 145-164
        
        Parameters
        ----------
        a : float array
        q : float array
        
        """
        _bcs.f90wrap_gram_schmidt(a=a, q=q)
    
    @staticmethod
    def orthogonalize_to_set_list(randvec, veclist):
        """
        orthogonalize_to_set_list(randvec, veclist)
        
        
        Defined at ../bcs/bcs.f90 lines 172-194
        
        Parameters
        ----------
        randvec : float array
        veclist : float array
        
        """
        _bcs.f90wrap_orthogonalize_to_set_list(randvec=randvec, veclist=veclist)
    
    @staticmethod
    def reweight_penalty(j0, epsilon, fxn=None):
        """
        reweight_penalty = reweight_penalty(j0, epsilon[, fxn])
        
        
        Defined at ../bcs/bcs.f90 lines 211-232
        
        Parameters
        ----------
        j0 : float
        epsilon : float
        fxn : str
        
        Returns
        -------
        reweight_penalty : float
        
        """
        reweight_penalty = _bcs.f90wrap_reweight_penalty(j0=j0, epsilon=epsilon, \
            fxn=fxn)
        return reweight_penalty
    
    @staticmethod
    def do_wrapped(full_pi, y, sigma2, eta, js, error_bars):
        """
        do_wrapped(full_pi, y, sigma2, eta, js, error_bars)
        
        
        Defined at ../bcs/bcs.f90 lines 292-307
        
        Parameters
        ----------
        full_pi : float array
        y : float array
        sigma2 : float
        eta : float
        js : float array
        error_bars : float array
        
        """
        _bcs.f90wrap_do_wrapped(full_pi=full_pi, y=y, sigma2=sigma2, eta=eta, js=js, \
            error_bars=error_bars)
    
    @staticmethod
    def bestsolution(hold_pi, hold_y, trackedjs):
        """
        bestsolution = bestsolution(hold_pi, hold_y, trackedjs)
        
        
        Defined at ../bcs/bcs.f90 lines 317-339
        
        Parameters
        ----------
        hold_pi : float array
        hold_y : float array
        trackedjs : float array
        
        Returns
        -------
        bestsolution : int
        
        """
        bestsolution = _bcs.f90wrap_bestsolution(hold_pi=hold_pi, hold_y=hold_y, \
            trackedjs=trackedjs)
        return bestsolution
    
    @staticmethod
    def isringing(idx, trackedell0s):
        """
        isringing = isringing(idx, trackedell0s)
        
        
        Defined at ../bcs/bcs.f90 lines 350-398
        
        Parameters
        ----------
        idx : int
        trackedell0s : int array
        
        Returns
        -------
        isringing : bool
        
        """
        isringing = _bcs.f90wrap_isringing(idx=idx, trackedell0s=trackedell0s)
        return isringing
    
    @staticmethod
    def partition_holdout_set(nfits, nsets, nholdout, fitlist, holdlist):
        """
        partition_holdout_set(nfits, nsets, nholdout, fitlist, holdlist)
        
        
        Defined at ../bcs/bcs.f90 lines 569-610
        
        Parameters
        ----------
        nfits : int
        nsets : int
        nholdout : int
        fitlist : int array
        holdlist : int array
        
        """
        _bcs.f90wrap_partition_holdout_set(nfits=nfits, nsets=nsets, nholdout=nholdout, \
            fitlist=fitlist, holdlist=holdlist)
    
    @staticmethod
    def newunit():
        """
        newunit, unit = newunit()
        
        
        Defined at ../bcs/bcs.f90 lines 615-628
        
        
        Returns
        -------
        newunit : int
        unit : int
        
        """
        newunit, unit = _bcs.f90wrap_newunit()
        return newunit, unit
    
    @staticmethod
    def file_exists(filename):
        """
        file_exists = file_exists(filename)
        
        
        Defined at ../bcs/bcs.f90 lines 633-636
        
        Parameters
        ----------
        filename : str
        
        Returns
        -------
        file_exists : bool
        
        """
        file_exists = _bcs.f90wrap_file_exists(filename=filename)
        return file_exists
    
    @staticmethod
    def value_count(line, length):
        """
        value_count = value_count(line, length)
        
        
        Defined at ../bcs/bcs.f90 lines 644-673
        
        Parameters
        ----------
        line : str
        length : int
        
        Returns
        -------
        value_count : int
        
        """
        value_count = _bcs.f90wrap_value_count(line=line, length=length)
        return value_count
    
    @staticmethod
    def linevalue_count(filename, n, commentchar):
        """
        nlines, nvalues = linevalue_count(filename, n, commentchar)
        
        
        Defined at ../bcs/bcs.f90 lines 683-719
        
        Parameters
        ----------
        filename : str
        n : int
        commentchar : str
        
        Returns
        -------
        nlines : int
        nvalues : int
        
        """
        nlines, nvalues = _bcs.f90wrap_linevalue_count(filename=filename, n=n, \
            commentchar=commentchar)
        return nlines, nvalues
    
    @staticmethod
    def write_results(js, fit_rms, fit_err, hold_rms_=None, hold_err_=None, \
        sigma2_=None):
        """
        write_results(js, fit_rms, fit_err[, hold_rms_, hold_err_, sigma2_])
        
        
        Defined at ../bcs/bcs.f90 lines 730-772
        
        Parameters
        ----------
        js : float array
        fit_rms : float array
        fit_err : float array
        hold_rms_ : float array
        hold_err_ : float array
        sigma2_ : float array
        
        """
        _bcs.f90wrap_write_results(js=js, fit_rms=fit_rms, fit_err=fit_err, \
            hold_rms_=hold_rms_, hold_err_=hold_err_, sigma2_=sigma2_)
    
    _dt_array_initialisers = []
    

bcs = Bcs()

class Laplace(f90wrap.runtime.FortranModule):
    """
    Module laplace
    
    
    Defined at ../bcs/laplace.f90 lines 2-773
    
    """
    class Laplace_Iterator(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=laplace_iterator)
        
        
        Defined at ../bcs/laplace.f90 lines 11-55
        
        """
        def __init__(self, handle=None):
            """
            self = Laplace_Iterator()
            
            
            Defined at ../bcs/laplace.f90 lines 11-55
            
            
            Returns
            -------
            this : Laplace_Iterator
            	Object to be constructed
            
            
            Automatically generated constructor for laplace_iterator
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            self._handle = _bcs.f90wrap_laplace_iterator_initialise()
        
        def __del__(self):
            """
            Destructor for class Laplace_Iterator
            
            
            Defined at ../bcs/laplace.f90 lines 11-55
            
            Parameters
            ----------
            this : Laplace_Iterator
            	Object to be destructed
            
            
            Automatically generated destructor for laplace_iterator
            """
            if self._alloc:
                _bcs.f90wrap_laplace_iterator_finalise(this=self._handle)
        
        @property
        def s(self):
            """
            Element s ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__s(self._handle)
            if array_handle in self._arrays:
                s = self._arrays[array_handle]
            else:
                s = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__s)
                self._arrays[array_handle] = s
            return s
        
        @s.setter
        def s(self, s):
            self.s[...] = s
        
        @property
        def q(self):
            """
            Element q ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__q(self._handle)
            if array_handle in self._arrays:
                q = self._arrays[array_handle]
            else:
                q = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__q)
                self._arrays[array_handle] = q
            return q
        
        @q.setter
        def q(self, q):
            self.q[...] = q
        
        @property
        def ls(self):
            """
            Element ls ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__ls(self._handle)
            if array_handle in self._arrays:
                ls = self._arrays[array_handle]
            else:
                ls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__ls)
                self._arrays[array_handle] = ls
            return ls
        
        @ls.setter
        def ls(self, ls):
            self.ls[...] = ls
        
        @property
        def lq(self):
            """
            Element lq ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__lq(self._handle)
            if array_handle in self._arrays:
                lq = self._arrays[array_handle]
            else:
                lq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__lq)
                self._arrays[array_handle] = lq
            return lq
        
        @lq.setter
        def lq(self, lq):
            self.lq[...] = lq
        
        @property
        def ninputs(self):
            """
            Element ninputs ftype=integer  pytype=int
            
            
            Defined at ../bcs/laplace.f90 line 20
            
            """
            return _bcs.f90wrap_laplace_iterator__get__ninputs(self._handle)
        
        @ninputs.setter
        def ninputs(self, ninputs):
            _bcs.f90wrap_laplace_iterator__set__ninputs(self._handle, ninputs)
        
        @property
        def ncorrs(self):
            """
            Element ncorrs ftype=integer  pytype=int
            
            
            Defined at ../bcs/laplace.f90 line 20
            
            """
            return _bcs.f90wrap_laplace_iterator__get__ncorrs(self._handle)
        
        @ncorrs.setter
        def ncorrs(self, ncorrs):
            _bcs.f90wrap_laplace_iterator__set__ncorrs(self._handle, ncorrs)
        
        @property
        def nused(self):
            """
            Element nused ftype=integer  pytype=int
            
            
            Defined at ../bcs/laplace.f90 line 20
            
            """
            return _bcs.f90wrap_laplace_iterator__get__nused(self._handle)
        
        @nused.setter
        def nused(self, nused):
            _bcs.f90wrap_laplace_iterator__set__nused(self._handle, nused)
        
        @property
        def sigma2(self):
            """
            Element sigma2 ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 24
            
            """
            return _bcs.f90wrap_laplace_iterator__get__sigma2(self._handle)
        
        @sigma2.setter
        def sigma2(self, sigma2):
            _bcs.f90wrap_laplace_iterator__set__sigma2(self._handle, sigma2)
        
        @property
        def eta(self):
            """
            Element eta ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 24
            
            """
            return _bcs.f90wrap_laplace_iterator__get__eta(self._handle)
        
        @eta.setter
        def eta(self, eta):
            _bcs.f90wrap_laplace_iterator__set__eta(self._handle, eta)
        
        @property
        def lambda(self):
            """
            Element lambda ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 24
            
            """
            return _bcs.f90wrap_laplace_iterator__get__lambda(self._handle)
        
        @lambda.setter
        def lambda(self, lambda):
            _bcs.f90wrap_laplace_iterator__set__lambda(self._handle, lambda)
        
        @property
        def alpha(self):
            """
            Element alpha ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 28
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__alpha(self._handle)
            if array_handle in self._arrays:
                alpha = self._arrays[array_handle]
            else:
                alpha = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__alpha)
                self._arrays[array_handle] = alpha
            return alpha
        
        @alpha.setter
        def alpha(self, alpha):
            self.alpha[...] = alpha
        
        @property
        def subphi(self):
            """
            Element subphi ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 28
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__subphi(self._handle)
            if array_handle in self._arrays:
                subphi = self._arrays[array_handle]
            else:
                subphi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__subphi)
                self._arrays[array_handle] = subphi
            return subphi
        
        @subphi.setter
        def subphi(self, subphi):
            self.subphi[...] = subphi
        
        @property
        def indices(self):
            """
            Element indices ftype=integer pytype=int
            
            
            Defined at ../bcs/laplace.f90 line 33
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__indices(self._handle)
            if array_handle in self._arrays:
                indices = self._arrays[array_handle]
            else:
                indices = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__indices)
                self._arrays[array_handle] = indices
            return indices
        
        @indices.setter
        def indices(self, indices):
            self.indices[...] = indices
        
        @property
        def selected(self):
            """
            Element selected ftype=integer pytype=int
            
            
            Defined at ../bcs/laplace.f90 line 33
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__selected(self._handle)
            if array_handle in self._arrays:
                selected = self._arrays[array_handle]
            else:
                selected = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__selected)
                self._arrays[array_handle] = selected
            return selected
        
        @selected.setter
        def selected(self, selected):
            self.selected[...] = selected
        
        @property
        def sigma(self):
            """
            Element sigma ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 39
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__sigma(self._handle)
            if array_handle in self._arrays:
                sigma = self._arrays[array_handle]
            else:
                sigma = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__sigma)
                self._arrays[array_handle] = sigma
            return sigma
        
        @sigma.setter
        def sigma(self, sigma):
            self.sigma[...] = sigma
        
        @property
        def mu(self):
            """
            Element mu ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 39
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__mu(self._handle)
            if array_handle in self._arrays:
                mu = self._arrays[array_handle]
            else:
                mu = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__mu)
                self._arrays[array_handle] = mu
            return mu
        
        @mu.setter
        def mu(self, mu):
            self.mu[...] = mu
        
        @property
        def theta(self):
            """
            Element theta ftype=real(dp) pytype=float
            
            
            Defined at ../bcs/laplace.f90 line 39
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bcs.f90wrap_laplace_iterator__array__theta(self._handle)
            if array_handle in self._arrays:
                theta = self._arrays[array_handle]
            else:
                theta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bcs.f90wrap_laplace_iterator__array__theta)
                self._arrays[array_handle] = theta
            return theta
        
        @theta.setter
        def theta(self, theta):
            self.theta[...] = theta
        
        @property
        def initialized(self):
            """
            Element initialized ftype=logical pytype=bool
            
            
            Defined at ../bcs/laplace.f90 line 41
            
            """
            return _bcs.f90wrap_laplace_iterator__get__initialized(self._handle)
        
        @initialized.setter
        def initialized(self, initialized):
            _bcs.f90wrap_laplace_iterator__set__initialized(self._handle, initialized)
        
        def __str__(self):
            ret = ['<laplace_iterator>{\n']
            ret.append('    s : ')
            ret.append(repr(self.s))
            ret.append(',\n    q : ')
            ret.append(repr(self.q))
            ret.append(',\n    ls : ')
            ret.append(repr(self.ls))
            ret.append(',\n    lq : ')
            ret.append(repr(self.lq))
            ret.append(',\n    ninputs : ')
            ret.append(repr(self.ninputs))
            ret.append(',\n    ncorrs : ')
            ret.append(repr(self.ncorrs))
            ret.append(',\n    nused : ')
            ret.append(repr(self.nused))
            ret.append(',\n    sigma2 : ')
            ret.append(repr(self.sigma2))
            ret.append(',\n    eta : ')
            ret.append(repr(self.eta))
            ret.append(',\n    lambda : ')
            ret.append(repr(self.lambda))
            ret.append(',\n    alpha : ')
            ret.append(repr(self.alpha))
            ret.append(',\n    subphi : ')
            ret.append(repr(self.subphi))
            ret.append(',\n    indices : ')
            ret.append(repr(self.indices))
            ret.append(',\n    selected : ')
            ret.append(repr(self.selected))
            ret.append(',\n    sigma : ')
            ret.append(repr(self.sigma))
            ret.append(',\n    mu : ')
            ret.append(repr(self.mu))
            ret.append(',\n    theta : ')
            ret.append(repr(self.theta))
            ret.append(',\n    initialized : ')
            ret.append(repr(self.initialized))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    @staticmethod
    def norm(vector):
        """
        norm = norm(vector)
        
        
        Defined at ../bcs/laplace.f90 lines 60-64
        
        Parameters
        ----------
        vector : float array
        
        Returns
        -------
        norm : float
        
        """
        norm = _bcs.f90wrap_norm(vector=vector)
        return norm
    
    _dt_array_initialisers = []
    

laplace = Laplace()

class Matrix_Sets(f90wrap.runtime.FortranModule):
    """
    Module matrix_sets
    
    
    Defined at ../bcs/matrix_sets.f90 lines 6-171
    
    """
    pass
    _dt_array_initialisers = []
    

matrix_sets = Matrix_Sets()

class Num_Types(f90wrap.runtime.FortranModule):
    """
    Module num_types
    
    
    Defined at ../bcs/num_types.f90 lines 2-16
    
    """
    pass
    _dt_array_initialisers = []
    

num_types = Num_Types()

