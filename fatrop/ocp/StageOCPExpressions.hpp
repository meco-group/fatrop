/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye <lander.vanroye@kuleuven.be>
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
#ifndef BASICOCPSAMPLERSINCLUDED
#define BASICOCPSAMPLERSINCLUDED
#include "function_evaluation/FunctionEvaluation.hpp"
#include "ocp/OCPAdapter.hpp"
#include <memory>
// todo change the name of this file to evaluators
namespace fatrop
{
    class StageExpression
    {
    public:
        virtual void evaluate(const double *u, const double *x, const double *global_params, const double *stage_params, double *res) const = 0;
        virtual int n_rows() const = 0;
        virtual int n_cols() const= 0;
        int size() const;
    };
    class IndexEpression : public StageExpression
    {
    public:
        IndexEpression(const bool control, const std::vector<int> offsets_in, const std::vector<int> offsets_out);
        void evaluate(const double *u, const double *x, const double *global_params, const double *stage_params, double *res) const override;
        int n_rows() const override;
        int n_cols() const override;

    private:
        const int _no_var;
        const std::vector<int> _offsets_in;
        const std::vector<int> _offsets_out;
        const bool _control;
    };
    class EvalBaseSE : public StageExpression
    {
    public:
        EvalBaseSE(const std::shared_ptr<EvalBase> &evalbase);
        void evaluate(const double *u, const double *x, const double *global_params, const double *stage_params, double *res) const override;
        int n_rows() const override;
        int n_cols() const override;

    private:
        std::shared_ptr<EvalBase> evalbase_;
        const int n_rows_;
        const int n_cols_;
    };
    class StageExpressionEvaluatorBase
    {
    public:
        virtual int size() const = 0;
        virtual int n_rows() const= 0;
        virtual int n_cols() const = 0;
        virtual int evaluate(const std::vector<double> &solution, const std::vector<double> &global_params, const std::vector<double> &stage_params, std::vector<double> &result) const = 0;
        std::vector<double> evaluate(const std::vector<double> &solution, const std::vector<double> &global_params, const std::vector<double> &stage_params) const
        {
            std::vector<double> res(size());
            evaluate(solution, global_params, stage_params, res);
            return res;
        };
        virtual ~StageExpressionEvaluatorBase() = default;

    protected:
        // virtual int Evaluate()
    };
    class StageControlGridSampler : public StageExpressionEvaluatorBase
    {
    public:
        StageControlGridSampler(int nu, int nx, int no_stage_params, int K, const std::shared_ptr<StageExpression> &eval);
        int evaluate(const std::vector<double> &solution, const std::vector<double> &global_params, const std::vector<double> &stage_params, std::vector<double> &result) const;
        int size() const;
        int n_rows() const;
        int n_cols() const;
        int K() const;
        const int nu;
        const int nx;
        const int no_stage_params;
        const int K_;

    private:
        std::shared_ptr<StageExpression> eval_;
    };

    class OCPTimeStepSampler : public StageExpressionEvaluatorBase
    {
    public:
        OCPTimeStepSampler(int nu, int nx, int no_stage_params, int K, int k, const std::shared_ptr<StageExpression> &eval) : nu(nu), nx(nx), no_stage_params(no_stage_params), K_(K), k_(k), eval_(eval)
        {
        }
        int evaluate(const std::vector<double> &solution, const std::vector<double> &global_params, const std::vector<double> &stage_params, std::vector<double> &result) const
        {
            eval_->evaluate(solution.data() + (k_ < K_ - 1 ? k_ * (nu + nx) : (K_ - 2) * (nu + nx)), solution.data() + (k_ < K_ - 1 ? k_ * (nu + nx) + nu : (K_ - 1) * (nu + nx)), global_params.data(), stage_params.data() + k_ * no_stage_params, result.data());
            return 0;
        }
        int size() const
        {
            return eval_->size();
        }
        int n_rows() const
        {
            return eval_->n_rows();
        }
        int n_cols() const
        {
            return eval_->n_cols();
        }
        const int nu;
        const int nx;
        const int no_stage_params;
        const int K_;
        const int k_;

    private:
        std::shared_ptr<StageExpression> eval_;
    };
    class StageExpressionEvaluatorFactory
    {
    public:
        StageExpressionEvaluatorFactory(const std::shared_ptr<StageExpression> &eval, int nu, int nx, int no_stage_params, int K) : nu(nu), nx(nx), no_stage_params(no_stage_params), K(K), eval_(eval){};
        // at_t0()
        OCPTimeStepSampler at_t0() const
        {
            return OCPTimeStepSampler(nu, nx, no_stage_params, K, 0, eval_);
        }
        // at_tf()
        OCPTimeStepSampler at_tf() const
        {
            return OCPTimeStepSampler(nu, nx, no_stage_params, K, K - 1, eval_);
        }
        // at_tk(int k)
        OCPTimeStepSampler at_tk(const int k) const
        {
            return OCPTimeStepSampler(nu, nx, no_stage_params, K, k, eval_);
        }
        // evaluate at control grid
        StageControlGridSampler at_control() const
        {
            return StageControlGridSampler(nu, nx, no_stage_params, K, eval_);
        }
        const int nu;
        const int nx;
        const int no_stage_params;
        const int K;
        std::shared_ptr<StageExpression> eval_;
    };

    class ParameterSetter
    {
    public:
        ParameterSetter(const std::vector<int> &offsets_in, const std::vector<int> &offsets_out, const int no_stage_params, const int no_var, const int K, const bool global);
        void set_value(std::vector<double> &global_params, std::vector<double> &stage_params, const double value[]);
        void set_value(std::vector<double> &global_params, std::vector<double> &stage_params, const std::initializer_list<double> il_);

    private:
        const std::vector<int> _offsets_in;
        const std::vector<int> _offsets_out;
        const int no_stage_params;

    protected:
        const int _no_var;
        const int K;
        const bool _global;
    };

} // namespace fatrop

#endif //  BASICOCPSAMPLERSINCLUDED