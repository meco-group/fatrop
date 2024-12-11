/*
 * Fatrop - A fast trajectory optimization solver
 *  Copyright (C) 2022 - 2024 Lander Vanroye, KU Leuven. All rights reserved.
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
#ifndef __fatrop_ocp_StageOCPApplication_hpp__
#define __fatrop_ocp_StageOCPApplication_hpp__
#include "OCPApplication.hpp"
namespace fatrop
{
    struct StageOCPSolution : public FatropSolution
    {
    public:
        StageOCPSolution(const std::shared_ptr<OCP> &app);
        void get_u(std::vector<double> &result) const;
        void get_x(std::vector<double> &result) const;
        int get_nx() const
        {
            return nx;
        };
        int get_nu() const
        {
            return nu;
        };
        int get_K() const
        {
            return K;
        };

    protected:
        StageOCPSolution();
        void set_dims(const OCPDims &dims);
        void set_parameters(const std::vector<double> &global_params, const std::vector<double> &stage_params);
        fatrop_int nx;
        fatrop_int nu;
        fatrop_int n_stage_params;
        fatrop_int n_global_params;
        fatrop_int K;

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
    class StageOCPApplication : public OCPAbstractApplication
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
            void set_value(const std::vector<double> il_)
            {
                set_value(il_.data());
            };

        private:
            const std::shared_ptr<OCPAdapter> adapter_;
        };

    public:
        StageExpressionEvaluatorFactory get_expression(const std::string &sampler_name);
        StageExpressionEvaluatorFactory get_expression_evaluator(const std::shared_ptr<StageExpression> &expr);
        AppParameterSetter get_parameter_setter(const std::string &setter_name);
        void set_value(const std::string &setter_name, const std::vector<double>& value);
        void build();
        fatrop_int optimize();
        const StageOCPSolution &last_solution() const;
        // const FatropSolution &last_solution() const { return last_solution(); };
        const std::vector<std::string> parameter_names() const;
        const std::vector<std::string> stage_expression_names() const;

    public:
        const fatrop_int nx_;
        const fatrop_int nu_;
        const fatrop_int n_stage_params_;
        const fatrop_int K_;

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
#endif // __fatrop_ocp_StageOCPApplication_hpp__