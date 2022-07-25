from fatropy cimport OCPBuilder

cdef class PyOCP:
    cdef OCPBuilder *myOCPBuilder  # Hold a C++ instance which we are wrapping

    def __cinit__(self, functions, specfile):
        self.myOCPBuilder = new OCPBuilder(functions.encode('utf-8'),specfile.encode('utf-8'))

    def Optimize(self):
        return self.myOCPBuilder.fatropalg.GetRawPtr().Optimize()

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
        retval = vector[double](nels)
        for ii in range(nels):
           retval[ii] = self.myOCPBuilder.fatropdata.GetRawPtr().x_curr.get_el(ii)
        return retval
    
    # Attribute access
    @property
    def x_next(self):
        nels = self.myOCPBuilder.fatropdata.GetRawPtr().x_next.nels()
        retval = vector[double](nels)
        for ii in range(nels):
           retval[ii] = self.myOCPBuilder.fatropdata.GetRawPtr().x_next.get_el(ii)
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
