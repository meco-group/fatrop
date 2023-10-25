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
  auto x_test = ocp.state();
  auto u_test = ocp.control(10);

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

  auto solver = SolverFatrop();
  solver.transcribe(ocp);
  auto dummy = std::vector<cs::MX>();
  auto dummy1 = std::vector<cs::MX>();

  cs::Function func = solver.to_function(ocp, dummy, dummy1);
  ocp.to_function(dummy, dummy1);
  // print shape of elements of dummy
  for (auto &el : dummy)
    std::cout << el.size1() << " " << el.size2() << std::endl;
  auto dummyin0_MX = cs::MX::sym("dummy", func.sparsity_in(0));
  auto dummyin1_MX = cs::MX::sym("dummy", func.sparsity_in(1));
  auto dummyin2_MX = cs::MX::sym("dummy", func.sparsity_in(2));
  auto dummyin0 = cs::DM::zeros(func.sparsity_in(0));
  auto dummyin1 = cs::DM::zeros(func.sparsity_in(1));
  auto dummyin2 = cs::DM::zeros(func.sparsity_in(2));
  auto res = func(std::vector<cs::MX>{dummyin0_MX, dummyin1_MX, dummyin2_MX});
  cs::Function func2 = cs::Function("func", {dummyin0_MX, dummyin1_MX, dummyin2_MX}, {res[0]});
  func(std::vector<cs::DM>{dummyin0, dummyin1, dummyin2});

  // ocp.sample(x);           // check which stages have x as state, evaluate and concatenate in a matrix
  // ocp.set_initial(x, 0.5); // set initial value of x for all stages that have x as state
  return 0;
}