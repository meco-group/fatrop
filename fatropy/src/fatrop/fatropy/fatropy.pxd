from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr 

# cdef extern from "SmartPtr.hpp" namespace "fatrop":
#     cdef cppclass shared_ptr[T]:
#         shared_ptr() except +
#         T* GetRawPtr()
cdef extern from "FatropStats.hpp" namespace "fatrop":
    cdef cppclass FatropStats:
        double compute_sd_time
        double duinf_time
        double eval_hess_time
        double eval_jac_time
        double eval_cv_time
        double eval_grad_time
        double eval_obj_time
        double initialization_time
        double time_total
        int eval_hess_count
        int eval_jac_count
        int eval_cv_count
        int eval_grad_count
        int eval_obj_count
        int iterations_count
        void Print()
cdef extern from "FatropAlg.hpp" namespace "fatrop":
    cdef cppclass FatropAlg:
        int Optimize()
        FatropStats GetStats()
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
    cdef cppclass OCPSolutionSampler:
        OCPSolutionSampler(const OCPSolutionSampler& cp)
        int Sample(vector[double]& sample)
        int Size()
        int n_rows()
        int n_cols()
        int K()
cdef extern from "OCPBuilder.hpp" namespace "fatrop":
    cdef cppclass ParameterSetter:
        void SetValue(const double* value)

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
        shared_ptr[OCPSolutionSampler] GetSampler(const string &sampler_name)
        shared_ptr[ParameterSetter] GetParameterSetter(const string &parameter_setter_name)