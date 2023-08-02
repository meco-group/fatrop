from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr 

# cdef extern from "SmartPtr.hpp" namespace "fatrop":
#     cdef cppclass shared_ptr[T]:
#         shared_ptr() except +
#         T* GetRawPtr()
cdef extern from "assign_ptr.hpp":
    cdef void assign_shared_ptr[T1, T2](shared_ptr[T1]& lhs, shared_ptr[T2]& rhs)

cdef extern from "LinearAlgebraBlasfeo.hpp" namespace "fatrop":
    cdef cppclass FatropVecBF:
        int offset()
        int nels()
        double get_el(const int ai)
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
        void print()
    

cdef extern from "StageOCPExpressions.hpp" namespace "fatrop":
    cdef cppclass StageExpressionEvaluatorBase:
        StageExpressionEvaluatorBase(const StageExpressionEvaluatorBase& cp)
        int n_rows() const
        int n_cols() const
        int K() const
        int size() const
        # vector[double] Eval(const FatropVecBF& solution, const vector[double]& global_params, const vector[double]& stage_params)
cdef extern from "StageOCPApplication.hpp" namespace "fatrop":
    cdef cppclass FatropSolution:
        pass
cdef extern from "StageOCPApplication.hpp" namespace "fatrop":
    cdef cppclass StageOCPSolution(FatropSolution):
        vector[double] evaluate(const StageExpressionEvaluatorBase &evaluator) const

cdef extern from "StageOCPExpressions.hpp" namespace "fatrop":
    cdef cppclass StageControlGridSampler(StageExpressionEvaluatorBase):
        StageControlGridSampler(const StageControlGridSampler& cp)
        # int Evaluate(const FatropVecBF& solution, const vector[double]& global_params, const vector[double]& stage_params, vector[double] &sample)
        # int n_rows()
        # int n_cols()
        # int K() const
        # int Size()
cdef extern from "StageOCPExpressions.hpp" namespace "fatrop":
    cdef cppclass StageExpressionEvaluatorFactory:
        StageExpressionEvaluatorFactory(const StageExpressionEvaluatorFactory& cp)
        StageControlGridSampler at_control() const
# cdef extern from "BasicOCPSamplers.hpp" namespace "fatrop":
#     cdef cppclass ParameterSetter:
#         void SetValue(vector[double]& global_params, vector[double]& stage_params, const double* value)
cdef extern from "StageOCPApplication.hpp" namespace "fatrop::StageOCPApplication":
    cdef cppclass AppParameterSetter:
        void set_value(const double* value)

cdef extern from "StageOCPApplication.hpp" namespace "fatrop":
    cdef cppclass StageOCPApplication:
        StageOCPApplication(const StageOCPApplication& cp)
        int optimize()
        void build()
        const FatropVecBF& last_solution_primal() const
        vector[double] &global_parameters()
        vector[double] &stage_parameters()
        # vector[double] &InitialGuessPrimal()
        void set_initial(vector[double] &initial_u, vector[double]& initial_x)
        void set_initial(const FatropSolution &initial_guess)
        # shared_ptr[OCPSolutionSampler] GetSampler(const string &sampler_name)
        AppParameterSetter get_parameter_setter(const string &sampler_name)
        vector[double] sample(const string &sampler_name) const
        FatropStats get_stats() const
        const int nx_
        const int nu_
        const int K_
        StageExpressionEvaluatorFactory get_expression(const string &sampler_name)
        const StageOCPSolution & last_stageocp_solution() const
        void set_option[T](const string &option_name, T value)
        void set_params(const vector[double] &global_params, const vector[double] &stage_params)

cdef extern from "StageOCPApplication.hpp" namespace "fatrop":
    cdef cppclass StageOCPApplicationFactory:
        @staticmethod
        StageOCPApplication from_rockit_interface(const string &functions, const string &json_spec_file) # except +


# cdef extern from "FatropAlg.hpp" namespace "fatrop":
#     cdef cppclass FatropAlg:
#         int Optimize()
#         FatropStats GetStats()
#         # double sd_time
#         # double hess_time
#         # double jac_time
#         # double cv_time
#         # double grad_time
#         # double obj_time
#         # double init_time
#         # double total_time


# cdef extern from "FatropData.hpp" namespace "fatrop":
#     cdef cppclass FatropData:
#         FatropVecBF x_curr
#         FatropVecBF x_next
#         int n_eqs
#         int n_ineqs      

# cdef extern from "FatropOptions.hpp" namespace "fatrop":
#     cdef cppclass FatropOptions:
#         int max_iter
#         int tol

# cdef extern from "FatropApplication.hpp" namespace "fatrop":
#     cdef cppclass FatropApplication:
#         void Initialize()
#         void Reset()
#         void SetBounds(const vector[double]& lower, const vector[double]& upper)
#         void SetInitial(const vector[double]& initial)
#         void GetSolution(vector[double]& sol)
#         void WarmStart()
#         int Optimize()

# cdef extern from "OCPBuilder.hpp" namespace "fatrop":
#     cdef cppclass OCP:
#         void SetParams(const vector[double] &stage_params_in, const vector[double] &global_params_in)
#         void SetInitial(const shared_ptr[FatropData] &fatropdata, vector[double] &initial_u, vector[double] &initial_x)

# cdef extern from "OCPBuilder.hpp" namespace "fatrop":
#     cdef cppclass OCPSolutionSampler:
#         OCPSolutionSampler(const OCPSolutionSampler& cp)
#         int Sample(vector[double]& sample)
#         int Size()
#         int n_rows()
#         int n_cols()
#         int K()
# cdef extern from "OCPBuilder.hpp" namespace "fatrop":
#     cdef cppclass ParameterSetter:
#         void SetValue(const double* value)

# cdef extern from "OCPBuilder.hpp" namespace "fatrop":
#     cdef cppclass OCPBuilder:
#         OCPBuilder(const string &functions, const string &json_spec_file) except +
#         shared_ptr[FatropAlg] fatropalg
#         shared_ptr[FatropData] fatropdata
#         shared_ptr[FatropOptions] fatropparams
#         shared_ptr[OCP] ocptempladapter
#         vector[double] initial_u
#         vector[double] initial_x
#         vector[double] lower
#         vector[double] upper
#         vector[double] lowerF
#         vector[double] upperF
#         void SetBounds()
#         void SetInitial()
#         shared_ptr[FatropApplication] Build()
#         shared_ptr[OCPSolutionSampler] GetSampler(const string &sampler_name)
#         shared_ptr[ParameterSetter] GetParameterSetter(const string &parameter_setter_name)
#         # int SampleMaxEnt(double alpha)