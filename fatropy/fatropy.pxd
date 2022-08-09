from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "SmartPtr.hpp" namespace "fatrop":
    cdef cppclass RefCountPtr[T]:
        RefCountPtr() except +
        T* GetRawPtr()

cdef extern from "FatropAlg.hpp" namespace "fatrop":
    cdef cppclass FatropAlg:
        int Optimize()
        double sd_time
        # double hess_time
        # double jac_time
        # double cv_time
        # double grad_time
        # double obj_time
        double init_time
        double total_time

cdef extern from "LinearAlgebraBlasfeo.hpp" namespace "fatrop":
    cdef cppclass FatropVecBF:
        int offset()
        int nels()
        double get_el(const int ai)

cdef extern from "FatropData.hpp" namespace "fatrop":
    cdef cppclass FatropData:
        FatropVecBF x_curr
        FatropVecBF x_next
        int n_eqs
        int n_ineqs      

cdef extern from "FatropParams.hpp" namespace "fatrop":
    cdef cppclass FatropParams:
        int max_iter
        int tol

cdef extern from "OCPBuilder.hpp" namespace "fatrop":
    cdef cppclass OCPBuilder:
        OCPBuilder(const string &functions, const string &json_spec_file) except +
        RefCountPtr[FatropAlg] fatropalg
        RefCountPtr[FatropData] fatropdata
        RefCountPtr[FatropParams] fatropparams
        vector[double] initial_u
        vector[double] initial_x
        vector[double] lower
        vector[double] upper
        vector[double] lowerF
        vector[double] upperF
        void SetBounds()
        void SetInitial()
