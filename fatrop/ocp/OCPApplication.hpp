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
#ifndef __fatrop_ocp_OCPApplication_hpp__
#define __fatrop_ocp_OCPApplication_hpp__
#include "fatrop/solver/FatropStats.hpp"
#include "fatrop/ocp/StageOCPExpressions.hpp"
#include <map>
#include <fstream>
#include <sstream>
#include <string>

namespace fatrop
{
    // forward declarations to hide the implementation details
    class Journaller;
    class FatropAlg;
    class StageOCP;
    class NLPDims;
    class FatropVecBF;
    class FatropNLP;
    class FatropOptions;
    class FatropData;
    class OCPAdapter;
    class FatropPrinter;
    class OCPAbstract;
    class OCPDims;
    class OCP;
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
        void build(const std::shared_ptr<FatropNLP> &nlp, const std::shared_ptr<FatropNLP> &nlp_resto);

    public:
        fatrop_int optimize() const;
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
        void set_initial(const FatropSolution &initial_guess) const;
        void set_initial(const std::vector<double> &initial_guess_primal_) const;
        const FatropOptions &get_options() const;

    protected:
        std::shared_ptr<FatropOptions> fatropoptions_;
        std::shared_ptr<FatropData> fatropdata_;
        std::shared_ptr<FatropNLP> nlp_;
        bool dirty = true;
        std::shared_ptr<FatropPrinter> printer_;
        friend class OcpSolverDriver;
        
    private:
        const std::shared_ptr<OCPAbstract> ocp_;
        std::shared_ptr<Journaller> journaller_;
        std::shared_ptr<FatropAlg> fatropalg_;
    };

    class OCPApplication : public NLPApplication
    {
    public:
        OCPApplication(const std::shared_ptr<OCP> &ocp);
        OCPApplication();
        void build();
    public:
        using NLPApplication::set_initial;
        void set_initial(std::vector<double> &initial_u, std::vector<double> &initial_x);
        OCPDims get_ocp_dims();
        std::shared_ptr<OCP> ocp_;
    };

    class OCPAbstractApplication : public OCPApplication
    {
    public:
        OCPAbstractApplication(const std::shared_ptr<OCPAbstract> &ocp);
        void set_params(const std::vector<double> &global_params, const std::vector<double> &stage_params);

    protected:
        std::shared_ptr<OCPAdapter> adapter;

    protected:
        std::vector<double> &global_parameters();
        std::vector<double> &stage_parameters();
    };
}
#endif // __fatrop_ocp_OCPApplication_hpp__