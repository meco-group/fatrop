#ifndef BASICOCPAPPLICATIONINCLUDED
#define BASICOCPAPPLICATIONINCLUDED
#include "ocp/StageOCP.hpp"
#include "solver/AlgBuilder.hpp"
#include "ocp/OCPAdapter.hpp"
#include "solver/FatropAlg.hpp"
#include "ocp/FatropOCP.hpp"
#include "ocp/FatropOCPBuilder.hpp"
#include "ocp/OCPAbstract.hpp"
#include "ocp/StageOCPExpressions.hpp"
#include "solver/FatropStats.hpp"
#include <map>
#include "json/json.h"
#include <fstream>
#include <sstream>
#include <string>

namespace fatrop
{
    class FatropSolution
    {
    public:
        // void GetPrimalSolution(vector<double> &result);
        // defautl copy constructor
        FatropSolution(const FatropSolution &other) = default;
        const std::vector<double> &primal_solution() { return sol_primal_; };

    protected:
        FatropSolution();
        void set_dims(const NLPDims &dims);
        void set_primal_solution(const FatropVecBF &sol);
        void set_solution(const FatropVecBF &sol_primal, const FatropVecBF &sol_dual, const FatropVecBF &sol_zL, const FatropVecBF &sol_zU);

    protected:
        std::vector<double> sol_primal_;
        std::vector<double> sol_dual_;
        std::vector<double> sol_zL_;
        std::vector<double> sol_zU_;
        friend class NLPApplication;
    };
    // TODO move this class to a separate file
    class NLPApplication
    {
    public:
        NLPApplication();

    protected:
        void build(const std::shared_ptr<FatropNLP> &nlp);

    public:
        int optimize();
        // TODO: make this protected and use last_solution instead and choose other name
        const FatropVecBF &last_solution_primal() const;
        const FatropVecBF &last_solution_dual() const;
        const FatropVecBF &last_solution_zL() const;
        const FatropVecBF &last_solution_zU() const;
        FatropVecBF &initial_guess_primal() const;
        FatropVecBF &initial_guess_dual() const;
        FatropVecBF &initial_guess_zL() const;
        FatropVecBF &initial_guess_zU() const;
        FatropStats get_stats() const;
        NLPDims get_nlp_dims();
        template <typename T>
        void set_option(const std::string &option_name, T value);

    public:
        void set_initial(const FatropSolution &initial_guess);
        const FatropOptions &get_options() const;

    protected:
        std::shared_ptr<FatropOptions> fatropoptions_;
        std::shared_ptr<FatropData> fatropdata_;
        std::shared_ptr<FatropNLP> nlp_;
        bool dirty = true;
        std::shared_ptr<FatropPrinter> printer_;

    private:
        const std::shared_ptr<OCPAbstract> ocp_;
        std::shared_ptr<Journaller> journaller_;
        std::shared_ptr<FatropAlg> fatropalg_;
    };

    // TODO move this class to a separate file
    class OCPApplication : public NLPApplication
    {
    public:
        OCPApplication(const std::shared_ptr<OCPAbstract> &ocp);

    public:
        void build();

    public:
        using NLPApplication::set_initial;
        void set_initial(std::vector<double> &initial_u, std::vector<double> &initial_x);
        OCPDims get_ocp_dims();
        void set_params(const std::vector<double> &global_params, const std::vector<double> &stage_params)
        {
            adapter->set_parameters(global_params, stage_params);
        }

    private:
        const std::shared_ptr<OCPAbstract> ocp_;

    protected:
        std::vector<double> &global_parameters();
        std::vector<double> &stage_parameters();
        std::shared_ptr<OCPAdapter> adapter;
    };

    struct StageOCPSolution : public FatropSolution
    {
    public:
        StageOCPSolution(const std::shared_ptr<OCP> &app);
        void get_u(std::vector<double> &result) const;
        void get_x(std::vector<double> &result) const;

    protected:
        StageOCPSolution();
        void set_dims(const OCPDims &dims);
        void set_parameters(const std::vector<double> &global_params, const std::vector<double> &stage_params);
        int nx;
        int nu;
        int n_stage_params;
        int n_global_params;
        int K;

    public:
        // todo make this deprecated, only use Eval
        void sample(const StageControlGridSampler &sampler, std::vector<double> &result) const;
        std::vector<double> evaluate(const StageExpressionEvaluatorBase &evaluator) const;
        void evaluate(const StageExpressionEvaluatorBase &evaluator, std::vector<double> &result) const;

    protected:
        std::vector<double> global_params;
        std::vector<double> stage_params;
        friend class StageOCPApplication;
    };

    // struct StageOCPApplicationAbstract : public OCPApplication
    // {
    //     StageOCPApplicationAbstract(const shared_ptr<StageOCP> &ocp);
    // };
    class StageOCPApplication : public OCPApplication
    {
    public:
        StageOCPApplication(const std::shared_ptr<StageOCP> &ocp);

        void set_initial_u(const std::vector<double> &initial_guess_u) const;
        void set_initial_x(const std::vector<double> &initial_guess_x) const;

    public:
        class AppParameterSetter : public ParameterSetter
        {
        public:
            AppParameterSetter(const std::shared_ptr<OCPAdapter> &adapter, const std::shared_ptr<ParameterSetter> &ps);

        public:
            void set_value(const double value[]);
            void set_value(const std::initializer_list<double> il_);

        private:
            const std::shared_ptr<OCPAdapter> adapter_;
        };

    public:
        StageExpressionEvaluatorFactory get_expression(const std::string &sampler_name);
        StageExpressionEvaluatorFactory get_expression_evaluator(const std::shared_ptr<StageExpression> &expr);
        AppParameterSetter get_parameter_setter(const std::string &setter_name);
        void build();
        int optimize();
        const StageOCPSolution &last_stageocp_solution() const;
        const FatropSolution &last_solution() const { return last_stageocp_solution(); };
        const std::vector<std::string> parameter_names() const;
        const std::vector<std::string> stage_expression_names() const;

    public:
        const int nx_;
        const int nu_;
        const int n_stage_params_;
        const int K_;

    private:
        StageOCPSolution last_solution_;

    protected:
        std::map<std::string, std::shared_ptr<StageExpression>> stage_expressions;
        std::map<std::string, std::shared_ptr<ParameterSetter>> param_setters;
        friend class StageOCPApplicationFactory;
    };

    class StageOCPApplicationFactory
    {
    public:
        static StageOCPApplication from_rockit_interface(const std::string &functions, const std::string &json_spec_file);
    };
}
#endif // BASICOCPAPPLICATIONINCLUDED