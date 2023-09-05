#include "fatrop-casadi-problem.hpp"
#include "fatrop-casadi-solver.hpp"
#include "single-stage-problem.hpp"
#include "single-stage-fatrop.hpp"
#include <casadi/casadi.hpp>
#include "ocp/StageOCPApplication.hpp"
#include <limits>
#include "utilities.hpp"
#define INF std::numeric_limits<double>::infinity()
using namespace fatrop;
using namespace fatrop::fatrop_casadi;
int main()
{
    // auto x = casadi::MX::sym("x", 2);
    // auto u = casadi::MX::sym("u", 2);
    // std::cout << ConstraintHelper(x<u).g_ineq << std::endl;
    StageProblem ocp;
    auto x = ocp.state(2);
    auto u = ocp.control();
    auto e = 1. - x(0) * x(0);
    double dt = .5;
    ocp.set_next(x, x + casadi::MX::vertcat({x(1), e * x(1) - x(0) + u}) * dt);
    ocp.subject_to(ocp.at_t0(x(0)) == 1.0);
    ocp.subject_to(ocp.at_t0(x(1)) == 0.0);
    ocp.subject_to((-0.25 <x(1))< INF, false, true);
    ocp.subject_to((-1.0 <x(1))< 1, false, true);

    ocp.add_objective(ocp.sum(u(0) * u(0) + x(0) * x(0) + x(1) * x(1), true, false));
    ocp.add_objective(ocp.at_tf(x(1) * x(1)));

    auto fatrop_method = std::make_shared<StageProblemFatropMethod>(&ocp);
    ocp.method = fatrop_method;
    fatrop_method -> transcribe(20);

    auto solver = std::make_shared<FatropCasadiSolver>(*fatrop_method);

    auto fatrop_solver = OCPApplication(solver);
    fatrop_solver.build();
    fatrop_solver.optimize();







    // typedef std::vector<double> vd;
    // // define variables
    // auto x = casadi::MX::sym("x", 2);
    // auto u = casadi::MX::sym("u", 1);
    // auto p_stage = casadi::MX::sym("p_stage", 0);
    // auto p_global = casadi::MX::sym("p_global", 0);
    // auto e = 1. - x(0) * x(0);
    // double dt = .5;
    // auto obj = u * u + x(0) * x(0) + x(1) * x(1);
    // auto obj_term = x(1) * x(1);
    // auto x_next = x + casadi::MX::vertcat({x(1), e * x(1) - x(0) + u}) * dt;
    // auto eq_initial = cs::MX::vertcat({x(0) -1.0, x(1)});
    // auto stage_initial = MicroStage::create(MicroStageSyms{x, u, p_stage, p_global}, obj, x_next, eq_initial, cs::MX(), std::vector<double>(), std::vector<double>());
    // auto stage_middle = MicroStage::create(MicroStageSyms{x, u, p_stage, p_global}, obj, x_next, cs::MX(), cs::MX::vertcat({x(1),u}), vd{-0.25, -1} , vd{INF, 1});
    // auto stage_terminal = MicroStage::create(MicroStageSyms{x, cs::MX(), p_stage, p_global}, obj_term, cs::MX(), cs::MX(), x(1), vd{-0.25} , vd{INF});

    // // set up the problem
    // auto problem = FatropCasadiProblem();
    // // add the microstages
    // problem.push_back(stage_initial);
    // for (int i = 0; i < 19; i++)
    // {
    //     problem.push_back(stage_middle);
    // }
    // problem.push_back(stage_terminal);

    // // // set up the solver
    // auto solver = std::make_shared<FatropCasadiSolver>(problem);
    // auto fatrop_solver = OCPApplication(solver);
    // fatrop_solver.build();
    // fatrop_solver.optimize();
    return 0;
}