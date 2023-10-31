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
    =----  initial ustage ----=
  */
  auto initial_ustage = ocp.new_ustage(); // states and controls are derived automatically
  /* constraint   */ initial_ustage.subject_to(x(0) == 1.);
  /* constraint   */ initial_ustage.subject_to(x(1) == 0.);
  /* dynamics     */ initial_ustage.set_next(x, x_next);
  /* objective    */ initial_ustage.add_objective(u * u + cs::MX::sumsqr(x));

  /*
    =----  middle stage ----=
  */
  auto middle_ustage = ocp.new_ustage(18); // 18 is the number of control intervals, states, controls and parameters are derived automatically
  /* constraints  */ middle_ustage.subject_to(-0.25 < x(1));
                     middle_ustage.subject_to((-1.0 < u) < 1);
  /* dynamics     */ middle_ustage.set_next(initial_ustage.dynamics());
  /* objective    */ middle_ustage.add_objective(u * u + cs::MX::sumsqr(x));

  /*
    =----  terminal stage ----=
  */
  auto terminal_ustage = ocp.new_ustage(); // the last stage also has x2 as state
  /* constraints  */ terminal_ustage.subject_to(-0.25 < x(1));
                     terminal_ustage.subject_to(x(1) == p);
  /* objective    */ terminal_ustage.add_objective(x(1)*x(1)+p);


  cs::Function ocp_func = ocp.to_function({p}, {ocp.at_t0(u), ocp.sample(x),p});
  auto ret = ocp_func({cs::DM(1.23)});
  std::cout << ret << std::endl;
  return 0;
}