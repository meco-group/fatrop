from fatropy cimport OCPBuilder
from fatropy cimport FatropAlg
from fatropy cimport FatropApplication
from fatropy cimport OCPSolutionSampler
from fatropy cimport ParameterSetter
from fatropy cimport FatropStats
from libcpp.memory cimport shared_ptr 
from libcpp.vector cimport vector
from cpython cimport array
import json
import numpy as np
from casadi import Callback
# class FatropCasFun(Callback):
#     def __init__(self, myocp, name, params, res):
#         Callback.__init__(self)
#         self.name = name
#         self.n_in = len(params)
#         self.n_out = len(res)
#         self.paramstetters
#         self.construct(name, self.n_in, self.n_out)
#         self.myocp = myocp
#     def get_n_in(self): return self.n_in
#     def get_n_out(self): return self.n_out
#     def eval(self, arg):
#         return self.fun(arg)

# interface the FatropStats struct to python
cdef class PyFatropStats:
    cdef FatropStats stats
    # double compute_sd_time
    @property
    def compute_sd_time(self):
        return self.stats.compute_sd_time
    # double duinf_time
    @property
    def duinf_time(self):
        return self.stats.duinf_time
    # double eval_hess_time
    @property
    def eval_hess_time(self):
        return self.stats.eval_hess_time
    # double eval_jac_time
    @property
    def eval_jac_time(self):
        return self.stats.eval_jac_time
    # double eval_cv_time
    @property
    def eval_cv_time(self):
        return self.stats.eval_cv_time
    # double eval_grad_time
    @property
    def eval_grad_time(self):
        return self.stats.eval_grad_time
    # double eval_obj_time
    @property
    def eval_obj_time(self):
        return self.stats.eval_obj_time
    # double initialization_time
    @property
    def initialization_time(self):
        return self.stats.initialization_time
    # double time_total
    @property
    def time_total(self):
        return self.stats.time_total
    # int eval_hess_count
    @property
    def eval_hess_count(self):
        return self.stats.eval_hess_count
    # int eval_jac_count
    @property
    def eval_jac_count(self):
        return self.stats.eval_jac_count 
    # double eval_hess_time
    @property
    def eval_hess_time(self):
        return self.stats.eval_hess_time
    # double eval_jac_time
    @property
    def eval_jac_time(self):
        return self.stats.eval_jac_time
    # double eval_cv_time
    @property
    def eval_cv_time(self):
        return self.stats.eval_cv_time
    # double eval_grad_time
    @property
    def eval_grad_time(self):
        return self.stats.eval_grad_time
    # double eval_obj_time
    @property
    def eval_obj_time(self):
        return self.stats.eval_obj_time
    # double initialization_time
    @property
    def initialization_time(self):
        return self.stats.initialization_time
    # double time_total
    @property
    def time_total(self):
        return self.stats.time_total
    # int eval_hess_count
    @property
    def eval_hess_count(self):
        return self.stats.eval_hess_count
    # int eval_jac_count
    @property
    def eval_jac_count(self):
        return self.stats.eval_jac_count
    # int eval_cv_count
    @property
    def eval_cv_count(self):
        return self.stats.eval_cv_count
    # int eval_grad_count
    @property
    def eval_grad_count(self):
        return self.stats.eval_grad_count
    # int eval_obj_count
    @property
    def eval_obj_count(self):
        return self.stats.eval_obj_count
    # int iterations_count
    @property
    def iterations_count(self):
        return self.stats.iterations_count
    def Print(self):
        self.stats.Print()
    


    
cdef class OCP:
    cdef OCPBuilder *myOCPBuilder  # Hold a C++ instance which we are wrapping
    cdef public dict OCPspecs # Public dict attribute to contain specs as defined in json file
    cdef shared_ptr[FatropApplication] myFatropApplication
    
    def __cinit__(self, functions, specfile):
        self.myOCPBuilder = new OCPBuilder(functions.encode('utf-8'),specfile.encode('utf-8'))
        specfile_object = open(specfile.encode('utf-8'),"r")
        self.OCPspecs = json.load(specfile_object)
        specfile_object.close()
        self.myFatropApplication = self.myOCPBuilder.Build()

    def Optimize(self):
        return self.myFatropApplication.get().Optimize()

    def SampleMaxEnt(self, alpha):
        return self.myOCPBuilder.SampleMaxEnt(alpha)
    def WarmStart(self):
        return self.myFatropApplication.get().WarmStart()
    def Sample(self, name):
        # retrieve sampler
        cdef shared_ptr[OCPSolutionSampler] sampler = self.myOCPBuilder.GetSampler(name.encode('utf-8'))
        # allocate buffer
        cdef vector[double] buffer = vector[double](sampler.get().Size())
        # use sampler
        sampler.get().Sample(buffer)
        n_rows = sampler.get().n_rows()
        n_cols = sampler.get().n_cols()
        K = sampler.get().K()
        # deallocate sampler
        if n_cols == 1:
            return np.asarray(buffer).reshape((K, n_rows))
        else:
            res = np.asarray(buffer).reshape((n_rows, n_cols, K), order = 'F')
            return np.moveaxis(res, [0,1,2], [1, 2, 0])
    def SetValue(self, name, double[::1] value):
        # retrieve parameter setter
        cdef shared_ptr[ParameterSetter] paramsetter = self.myOCPBuilder.GetParameterSetter(name.encode('utf-8'))
        paramsetter.get().SetValue(&value[0])
        return
    def to_function(self, params, samplers):
        pass
    def GetStats(self):
        res = PyFatropStats()
        res.stats = self.myOCPBuilder.fatropalg.get().GetStats()
        return res 
    # @property
    # def sd_time(self):
    #     return self.myOCPBuilder.fatropalg.get().sd_time

    # @property
    # def hess_time(self):
    #     return self.myOCPBuilder.fatropalg.get().hess_time

    # @property
    # def jac_time(self):
    #     return self.myOCPBuilder.fatropalg.get().jac_time

    # @property
    # def cv_time(self):
    #     return self.myOCPBuilder.fatropalg.get().cv_time

    # @property
    # def grad_time(self):
    #     return self.myOCPBuilder.fatropalg.get().grad_time

    # @property
    # def obj_time(self):
    #     return self.myOCPBuilder.fatropalg.get().obj_time

    # @property
    # def init_time(self):
    #     return self.myOCPBuilder.fatropalg.get().init_time

    # @property
    # def total_time(self):
    #     return self.myOCPBuilder.fatropalg.get().total_time

    # def SetBounds(self):
    #     self.myOCPBuilder.SetBounds()

    # def SetInitial(self):
    #     self.myOCPBuilder.SetInitial()
    def SetParams(self, stage_params_in, global_params_in):
        self.myOCPBuilder.ocptempladapter.get().SetParams(stage_params_in, global_params_in)
    def SetInitial(self, initial_u, initial_x):
        self.myOCPBuilder.ocptempladapter.get().SetInitial(self.myOCPBuilder.fatropdata, initial_u, initial_x)

    # Attribute access
    @property
    def initial_u(self):
        return self.myOCPBuilder.initial_u
    @initial_u.setter
    def initial_u(self, initial_u):
        self.myOCPBuilder.initial_u = initial_u

    # Attribute access
    @property
    def initial_x(self):
        return self.myOCPBuilder.initial_x
    @initial_x.setter
    def initial_x(self, initial_x):
        self.myOCPBuilder.initial_x = initial_x
    
    # Attribute access
    @property
    def lower(self):
        return self.myOCPBuilder.lower
    @lower.setter
    def lower(self, lower):
        self.myOCPBuilder.lower = lower

    # Attribute access
    @property
    def upper(self):
        return self.myOCPBuilder.upper
    @upper.setter
    def upper(self, upper):
        self.myOCPBuilder.upper = upper

    # Attribute access
    @property
    def lowerF(self):
        return self.myOCPBuilder.lowerF
    @lowerF.setter
    def lowerF(self, lowerF):
        self.myOCPBuilder.lowerF = lowerF

    # Attribute access
    @property
    def upperF(self):
        return self.myOCPBuilder.upperF
    @upperF.setter
    def upperF(self, upperF):
        self.myOCPBuilder.upperF = upperF

    # Attribute access
    @property
    def x_curr(self):
        nels = self.myOCPBuilder.fatropdata.get().x_curr.nels()
        retval = np.empty(nels)
        for ii in range(nels):
           retval[ii] = self.myOCPBuilder.fatropdata.get().x_curr.get_el(ii)
        return retval
    
    # Attribute access
    @property
    def x_next(self):
        nels = self.myOCPBuilder.fatropdata.get().x_next.nels()
        retval = np.empty(nels)
        for ii in range(nels):
           retval[ii] = self.myOCPBuilder.fatropdata.get().x_next.get_el(ii)
        return retval

    # Attribute access
    @property
    # TODO make this more efficient
    def u0_sol(self):
        nu = self.OCPspecs["nu"]
        retval = np.empty(nu)
        for ii in range(nu):
           retval[ii] = self.myOCPBuilder.fatropdata.get().x_curr.get_el(ii)
        return retval

    # Attribute access
    @property
    # TODO make this more efficient
    def u_sol(self):
        nu = self.OCPspecs["nu"]
        nx_plus_nu = self.OCPspecs["nx"]+nu
        K = self.OCPspecs["K"]
        retval = np.empty((nu,K-1))
        for ii in range(K-1):
            for jj in range(nu):               
                retval[jj,ii] = self.myOCPBuilder.fatropdata.get().x_curr.get_el(jj+ii*(nx_plus_nu))
        return retval

    @property
    # TODO make this more efficient
    def x_sol(self):
        nx = self.OCPspecs["nx"]
        nu = self.OCPspecs["nu"]
        nx_plus_nu = nx+nu
        K = self.OCPspecs["K"]
        retval = np.ones((nx,K))
        for ii in range(K-1):
            for jj in range(nx):               
                retval[jj,ii] = self.myOCPBuilder.fatropdata.get().x_curr.get_el(nu+jj+ii*(nx_plus_nu))
        for jj in range(nx):
            retval[jj,K-1] = self.myOCPBuilder.fatropdata.get().x_curr.get_el(jj+(K-1)*(nx_plus_nu))
        return retval
    

    # Attribute access
    @property
    def n_eqs(self):
        return self.myOCPBuilder.fatropdata.get().n_eqs

    # Attribute access
    @property
    def n_ineqs(self):
        return self.myOCPBuilder.fatropdata.get().n_ineqs

    def __dealloc__(self):
        del self.myOCPBuilder
