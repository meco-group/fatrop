#include "fatrop-casadi-problem.hpp"
#include "fatrop-casadi-solver.hpp"
#include <casadi/casadi.hpp>
#include "ocp/StageOCPApplication.hpp"
using namespace fatrop;
using namespace fatrop::specification;
int main()
{
    // define variables
    auto x = casadi::MX::sym("x", 2);
    auto u = casadi::MX::sym("u", 1);
    auto p_stage = casadi::MX::sym("p_stage", 0);
    auto p_global = casadi::MX::sym("p_global", 0);
    auto e = 1. - x(0)*x(0) - x(1)*x(1);
    double dt = 0.2;
    auto lb_middle = std::vector<double>{-100, -100};
    auto ub_middle = std::vector<double>{100, 100};
    // define functions
    auto x_next = casadi::MX::vertcat({e*x(0) -x(1) + u, x(0)})*dt;
    auto obj = u*u;
    auto L = casadi::Function("L", {x, u, p_stage, p_global}, {obj});
    auto L_term = casadi::Function("L", {x, casadi::MX::sym("dummy", 0), p_stage, p_global}, {0});
    auto dynamics = casadi::Function("dynamics", {x, u, p_stage, p_global}, {x_next});
    auto g_eq0 = casadi::Function("g_eq0", {x, u, p_stage, p_global}, {cs::MX::vertcat({x(0), x(1)})});
    auto g_term = casadi::Function("g_eq0", {x, casadi::MX::sym("dummy",0), p_stage, p_global}, {cs::MX::vertcat({x(0) -1.0, x(1) -1.0})});
    auto g_ineq_middle = casadi::Function("g_ineq_middle", {x, u, p_stage, p_global}, {cs::MX::vertcat({x(0), x(1)})});
    auto g_ineq_term = casadi::Function("g_ineq_term", {x, casadi::MX::sym("dummy", 0), p_stage, p_global}, {cs::MX::vertcat({x(0), x(1)})});

    // define StageQuantities
    auto stage_initial = StageQuantities::create(L, dynamics, g_eq0, casadi::Function(), std::vector<double>(), std::vector<double>()).expand();
    auto stage_middle = StageQuantities::create(L, dynamics, casadi::Function(), g_ineq_middle, lb_middle, ub_middle).expand();
    auto stage_terminal = StageQuantities::create(L_term, casadi::Function(), g_term, g_ineq_term, lb_middle, ub_middle).expand();



    // auto stage_initial = StageQuantities::create(L, dynamics, g_eq0, casadi::Function(), std::vector<double>(), std::vector<double>()).expand();
    // auto stage_middle = StageQuantities::create(L, dynamics, casadi::Function(), casadi::Function(), std::vector<double>(), std::vector<double>()).expand();
    // auto stage_terminal = StageQuantities::create(L_term, casadi::Function(), casadi::Function(), casadi::Function(), std::vector<double>(), std::vector<double>()).expand();

    // set up the problem
    auto problem = FatropCasadiProblem();
    // add the stage quantities
    problem.push_back(stage_initial);
    for(int i =0; i< 10; i++)
    {
        problem.push_back(stage_middle);
    }
    problem.push_back(stage_terminal);

    // // set up the solver
    auto solver = std::make_shared<FatropCasadiSolver>(problem);
    auto fatrop_solver = OCPApplication(solver);
    fatrop_solver.build();
    fatrop_solver.optimize();
    return 0;
}