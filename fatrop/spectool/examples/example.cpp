#include "fatrop/spectool/spectool.hpp"
#include "fatrop/ocp/StageOCPApplication.hpp"
using namespace std;
using namespace fatrop;
using namespace fatrop::spectool;

int main()
{
  auto ocp = Ocp();
  auto x = ocp.state(2);
  auto u = ocp.control();
  auto p = ocp.parameter();
  auto dt = ocp.hybrid();

  auto e = 1. - x(0) * x(0);
  // double dt = .5;
  auto x_next = x + casadi::MX::vertcat({x(1), e * x(1) - x(0) + u}) * dt;


  // first stage
  auto stage = ocp.new_stage(20);
  stage.at_t0().subject_to(x(0) == 1.);
  stage.at_t0().subject_to(x(1) == 0.);
  stage.add_objective(u * u + casadi::MX::sumsqr(x), at::t0, at::mid);
  stage.add_objective(casadi::MX::sumsqr(x), at::tf);
  stage.add_objective(dt*dt, at::t0, at::mid, at::tf);
  stage.set_next(x, x_next);
  stage.set_next(dt, dt);


  // connect stage -> stage2
  stage.at_tf().set_next(x, x_next);


  // second stage
  auto stage2 = ocp.new_stage(20);
  stage2.add_objective(u * u + casadi::MX::sumsqr(x-1), at::t0, at::mid);
  stage2.add_objective(dt*dt, at::t0, at::mid, at::tf);
  stage2.add_objective(casadi::MX::sumsqr(x), at::tf);
  stage2.set_next(x, x_next);
  stage2.subject_to(-0.25 < x(1), at::tf);
  stage2.set_next(dt, dt);


  ocp.set_initial(dt, 0.5);

  // auto my_ustage = uStage(20);
  // my_ustage.register_state({x});
  // my_ustage.register_control({u});
  // my_ustage.add_objective(dt*dt);


  // /*
  //   =----  initial ustage ----=
  // */
  // auto initial_ustage = ocp.new_ustage(); // states and controls are derived automatically
  // /* constraint   */ initial_ustage.subject_to(x(0) == 1.);
  // /* constraint   */ initial_ustage.subject_to(x(1) == 0.);
  // /* dynamics     */ initial_ustage.set_next(x, x_next);
  // /* objective    */ initial_ustage.add_objective(u * u + cs::MX::sumsqr(x));

  // /*
  //   =----  middle stage ----=
  // */
  // auto middle_ustage = ocp.new_ustage(18); // 18 is the number of control intervals, states, controls and parameters are derived automatically
  // /* constraints  */ middle_ustage.subject_to(-0.25 < x(1));
  //                    middle_ustage.subject_to((-1.0 < u) < 1);
  // /* dynamics     */ middle_ustage.set_next(initial_ustage.dynamics());
  // /* objective    */ middle_ustage.add_objective(u * u + cs::MX::sumsqr(x));

  // /*
  //   =----  terminal stage ----=
  // */
  // auto terminal_ustage = ocp.new_ustage(); // the last stage also has x2 as state
  // /* constraints  */ terminal_ustage.subject_to(-0.25 < x(1));
  //                    terminal_ustage.subject_to(x(1) == p);
  // /* objective    */ terminal_ustage.add_objective(x(1)*x(1)+p);


  cs::Function ocp_func = ocp.to_function({p}, {ocp.at_t0(u), ocp.sample(x).second,p, ocp.at_t0(dt), ocp.at_tf(dt)}, {{"jit",false}});
  auto ret = ocp_func({cs::DM(1.23)});
  std::cout << ret << std::endl;
  return 0;
}