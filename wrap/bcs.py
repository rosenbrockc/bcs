import _bcs
import f90wrap.runtime

class Num_Types(f90wrap.runtime.FortranModule):
    """
    Module num_types
    
    
    Defined at num_types.f90 lines 2-16
    
    """
    pass
    _dt_array_initialisers = []
    

num_types = Num_Types()

class Laplace(f90wrap.runtime.FortranModule):
    """
    Module laplace
    
    
    Defined at laplace.f90 lines 2-768
    
    """
    class Laplace_Iterator(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=laplace_iterator)
        
        
        Defined at laplace.f90 lines 11-56
        
        """
        def __init__(self, handle=None):
            """
            self = Laplace_Iterator()
            
            
            Defined at laplace.f90 lines 11-56
            
            
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
            
            
            Defined at laplace.f90 lines 11-56
            
            Parameters
            ----------
            this : Laplace_Iterator
            	Object to be destructed
            
            
            Automatically generated destructor for laplace_iterator
            """
            if self._alloc:
                _bcs.f90wrap_laplace_iterator_finalise(this=self._handle)
        
        _dt_array_initialisers = []
        
    @staticmethod
    def norm(vector):
        """
        norm = norm(vector)
        
        
        Defined at laplace.f90 lines 61-65
        
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
    
    
    Defined at matrix_sets.f90 lines 6-171
    
    """
    pass
    _dt_array_initialisers = []
    

matrix_sets = Matrix_Sets()

class Bcs(f90wrap.runtime.FortranModule):
    """
    Module bcs
    
    
    Defined at bcs.f90 lines 9-923
    
    """
    @staticmethod
    def gram_schmidt(a, q):
        """
        gram_schmidt(a, q)
        
        
        Defined at bcs.f90 lines 145-164
        
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
        
        
        Defined at bcs.f90 lines 172-194
        
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
        
        
        Defined at bcs.f90 lines 211-232
        
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
        
        
        Defined at bcs.f90 lines 292-307
        
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
    def partition_holdout_set(nfits, nsets, nholdout, fitlist, holdlist):
        """
        partition_holdout_set(nfits, nsets, nholdout, fitlist, holdlist)
        
        
        Defined at bcs.f90 lines 443-484
        
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
    def newunit(unit=None):
        """
        newunit = newunit([unit])
        
        
        Defined at bcs.f90 lines 489-502
        
        Parameters
        ----------
        unit : int
        
        Returns
        -------
        newunit : int
        
        """
        newunit = _bcs.f90wrap_newunit(unit=unit)
        return newunit
    
    @staticmethod
    def file_exists(filename):
        """
        file_exists = file_exists(filename)
        
        
        Defined at bcs.f90 lines 507-510
        
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
        
        
        Defined at bcs.f90 lines 518-547
        
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
    def linevalue_count(filename, n, commentchar, nlines, nvalues):
        """
        linevalue_count(filename, n, commentchar, nlines, nvalues)
        
        
        Defined at bcs.f90 lines 557-593
        
        Parameters
        ----------
        filename : str
        n : int
        commentchar : str
        nlines : int
        nvalues : int
        
        """
        _bcs.f90wrap_linevalue_count(filename=filename, n=n, commentchar=commentchar, \
            nlines=nlines, nvalues=nvalues)
    
    @staticmethod
    def write_results(js, fit_rms, fit_err, hold_rms_=None, hold_err_=None, \
        sigma2_=None):
        """
        write_results(js, fit_rms, fit_err[, hold_rms_, hold_err_, sigma2_])
        
        
        Defined at bcs.f90 lines 604-646
        
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

