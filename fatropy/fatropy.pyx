from fatropy cimport OCPBuilder
from fatropy cimport FatropAlg
import json
import numpy as np

cdef class OCP:
    cdef OCPBuilder *myOCPBuilder  # Hold a C++ instance which we are wrapping
    cdef public dict OCPspecs # Public dict attribute to contain specs as defined in json file
    
    def __cinit__(self, functions, specfile):
        self.myOCPBuilder = new OCPBuilder(functions.encode('utf-8'),specfile.encode('utf-8'))
        specfile_object = open(specfile.encode('utf-8'),"r")
        self.OCPspecs = json.load(specfile_object)
        specfile_object.close()

    def Optimize(self):
        return self.myOCPBuilder.fatropalg.GetRawPtr().Optimize()

    @property
    def sd_time(self):
        return self.myOCPBuilder.fatropalg.GetRawPtr().sd_time

    # @property
    # def hess_time(self):
    #     return self.myOCPBuilder.fatropalg.GetRawPtr().hess_time

    # @property
    # def jac_time(self):
    #     return self.myOCPBuilder.fatropalg.GetRawPtr().jac_time

    # @property
    # def cv_time(self):
    #     return self.myOCPBuilder.fatropalg.GetRawPtr().cv_time

    # @property
    # def grad_time(self):
    #     return self.myOCPBuilder.fatropalg.GetRawPtr().grad_time

    # @property
    # def obj_time(self):
    #     return self.myOCPBuilder.fatropalg.GetRawPtr().obj_time

    @property
    def init_time(self):
        return self.myOCPBuilder.fatropalg.GetRawPtr().init_time

    @property
    def total_time(self):
        return self.myOCPBuilder.fatropalg.GetRawPtr().total_time

    def SetBounds(self):
        self.myOCPBuilder.SetBounds()

    def SetInitial(self):
        self.myOCPBuilder.SetInitial()

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
        nels = self.myOCPBuilder.fatropdata.GetRawPtr().x_curr.nels()
        retval = np.empty(nels)
        for ii in range(nels):
           retval[ii] = self.myOCPBuilder.fatropdata.GetRawPtr().x_curr.get_el(ii)
        return retval
    
    # Attribute access
    @property
    def x_next(self):
        nels = self.myOCPBuilder.fatropdata.GetRawPtr().x_next.nels()
        retval = np.empty(nels)
        for ii in range(nels):
           retval[ii] = self.myOCPBuilder.fatropdata.GetRawPtr().x_next.get_el(ii)
        return retval

    # Attribute access
    @property
    # TODO make this more efficient
    def u0_sol(self):
        nu = self.OCPspecs["nu"]
        retval = np.empty(nu)
        for ii in range(nu):
           retval[ii] = self.myOCPBuilder.fatropdata.GetRawPtr().x_next.get_el(ii)
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
                retval[jj,ii] = self.myOCPBuilder.fatropdata.GetRawPtr().x_next.get_el(jj+ii*(nx_plus_nu))
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
                retval[jj,ii] = self.myOCPBuilder.fatropdata.GetRawPtr().x_next.get_el(nu+jj+ii*(nx_plus_nu))
        for jj in range(nx):
            retval[jj,K-1] = self.myOCPBuilder.fatropdata.GetRawPtr().x_next.get_el(jj+(K-1)*(nx_plus_nu))
        return retval

    # Attribute access
    @property
    def n_eqs(self):
        return self.myOCPBuilder.fatropdata.GetRawPtr().n_eqs

    # Attribute access
    @property
    def n_ineqs(self):
        return self.myOCPBuilder.fatropdata.GetRawPtr().n_ineqs

    def __dealloc__(self):
        del self.myOCPBuilder
