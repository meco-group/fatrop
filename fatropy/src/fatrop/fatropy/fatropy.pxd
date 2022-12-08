from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr 

# cdef extern from "SmartPtr.hpp" namespace "fatrop":
#     cdef cppclass shared_ptr[T]:
#         shared_ptr() except +
#         T* GetRawPtr()

cdef extern from "FatropAlg.hpp" namespace "fatrop":
    cdef cppclass FatropAlg:
        int Optimize()
        # double sd_time
        # double hess_time
        # double jac_time
        # double cv_time
        # double grad_time
        # double obj_time
        # double init_time
        # double total_time

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

cdef extern from "FatropApplication.hpp" namespace "fatrop":
    cdef cppclass FatropApplication:
        void Initialize()
        void Reset()
        void SetBounds(const vector[double]& lower, const vector[double]& upper)
        void SetInitial(const vector[double]& initial)
        void GetSolution(vector[double]& sol)
        void WarmStart()
        int Optimize()

cdef extern from "OCPBuilder.hpp" namespace "fatrop":
    cdef cppclass OCP:
        void SetParams(const vector[double] &stage_params_in, const vector[double] &global_params_in)
        void SetInitial(const shared_ptr[FatropData] &fatropdata, vector[double] &initial_u, vector[double] &initial_x)
    
cdef extern from "OCPBuilder.hpp" namespace "fatrop":
    cdef cppclass OCPBuilder:
        OCPBuilder(const string &functions, const string &json_spec_file) except +
        shared_ptr[FatropAlg] fatropalg
        shared_ptr[FatropData] fatropdata
        shared_ptr[FatropParams] fatropparams
        shared_ptr[OCP] ocptempladapter
        vector[double] initial_u
        vector[double] initial_x
        vector[double] lower
        vector[double] upper
        vector[double] lowerF
        vector[double] upperF
        void SetBounds()
        void SetInitial()
        shared_ptr[FatropApplication] Build()