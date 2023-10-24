#include "spectrop.hpp"
#include "ocp/StageOCPApplication.hpp"
using namespace std;
using namespace fatrop;
using namespace fatrop::spectrop;

int main()
{
  auto ocp = Ocp();
  auto x = ocp.state(2);
  auto u = ocp.control();
  auto p = ocp.parameter();

  auto e = 1. - x(0) * x(0);
  double dt = .5;
  auto x_next = x + casadi::MX::vertcat({x(1), e * x(1) - x(0) + u}) * dt;

  /*
    =----  initial stage ----=
  */
  auto initial_stage = ocp.new_stage(); // states and controls are derived automatically
  /* constraint   */ initial_stage.subject_to(x(0) == 1.);
  /* constraint   */ initial_stage.subject_to(x(1) == 0.);
  /* dynamics     */ initial_stage.set_next(x, x_next);
  /* objective    */ initial_stage.add_objective(u * u + cs::MX::sumsqr(x));

  /*
    =----  middle stage ----=
  */
  auto middle_stage = ocp.new_stage(19); // 18 is the number of nodes, states, controls and parameters are derived automatically
  /* constraints  */
  middle_stage.subject_to(-0.25 < x(1));
  middle_stage.subject_to((-1.0 < u) < 1);
  /* dynamics     */ middle_stage.set_next(x, x_next);
  /* objective    */ middle_stage.add_objective(u * u + cs::MX::sumsqr(x));

  /*
    =----  terminal stage ----=
  */
  auto terminal_stage = ocp.new_stage(); // the last stage also has x2 as state
  /* constraints  */ terminal_stage.subject_to(-0.25 < x(1));
  /* objective    */ terminal_stage.add_objective(x(1)*x(1));

  auto fatrop_impl = std::make_shared<FatropOcpImpl>(ocp);
  auto fatrop_solver = OCPApplication(fatrop_impl);
  fatrop_solver.build();
  fatrop_solver.optimize();

  // ocp.sample(x);           // check which stages have x as state, evaluate and concatenate in a matrix
  // ocp.set_initial(x, 0.5); // set initial value of x for all stages that have x as state
  return 0;
}