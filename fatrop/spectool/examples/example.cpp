#include "fatrop/spectool/spectool.hpp"
#include "fatrop/fatrop.hpp"
// #include "fatropy/src/python_func_import.hpp"
using namespace std;
using namespace fatrop;
using namespace fatrop::spectool;

int main()
{
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
    stage.add_objective(u * u + casadi::MX::sumsqr(x), at::t0, at::mid);
    stage.add_objective(casadi::MX::sumsqr(x), at::tf);
    stage.add_objective(dt * dt, at::t0, at::mid, at::tf);
    stage.set_next(x, x_next);
    stage.set_next(dt, dt);

    // connect stage -> stage2
    stage.at_tf().set_next(x, x_next);

    // second stage
    auto stage2 = ocp.new_stage(20);
    stage2.add_objective(u * u + casadi::MX::sumsqr(x - 1), at::t0, at::mid);
    stage2.add_objective(dt * dt, at::t0, at::mid, at::tf);
    stage2.add_objective(casadi::MX::sumsqr(x), at::tf);
    stage2.set_next(x, x_next);
    stage2.subject_to(-0.25 < x(1), at::tf);
    stage2.set_next(dt, dt);
    auto dum = ocp.control();
    ocp.at_tf().subject_to(dum == 0.5);

    ocp.set_initial(dt, 0.5);

    ocp.solver("fatrop");
    cs::Function ocp_func = ocp.to_function("example_ocp", {p}, {ocp.at_t0(u), ocp.sample(x).second, p, ocp.at_t0(dt), ocp.at_tf(dt), ocp.at_tf(dum)});
    auto ret = ocp_func({cs::DM(1.23)});
    std::cout << ret << std::endl;
  }

  // it's somewhat of a hack, but fatrop can also be used to solve small-scale dense problems, by using only one ustage.
  // {
  //   auto ss = Ocp();
  //   auto x = ss.state();
  //   auto y = ss.state();
  //   auto z = ss.state();
  //   auto stage = ss.new_ustage(1);
  //   stage.add_objective(y * sin(x) + x * cos(y) + z * z);
  //   stage.subject_to(x + y + z == 0.);
  //   ss.solver("fatrop");
  //   auto funcc = ss.to_function("example_dense_smallscale", {}, {stage.at_t0(x)});
  //   std::cout << funcc(std::vector<cs::DM>{}) << std::endl;
  //   std::cout << funcc(std::vector<cs::DM>{}) << std::endl;
  // }
  // {

  //   auto ss = Ocp();
  //   auto x = ss.state(10);
  //   auto stage = ss.new_ustage(1);
  //   stage.add_objective(sum1(sin(x)));
  //   ss.solver("fatrop");
  //   auto funcc = ss.to_function("example_dense", {}, {stage.at_t0(x)});
  //   std::cout << funcc(std::vector<cs::DM>{}) << std::endl;
  //   std::cout << funcc(std::vector<cs::DM>{}) << std::endl;
  // }
  {
    auto ss = Ocp();
    auto x = ss.control(10);
    auto ustage = uStage();
    ustage.register_control({x});
    ustage.add_objective(sum1(sin(x)));
    ustage.subject_to(sumsqr(x) == 1.);
    auto ustage_dup = ustage.clone();
    ss.add_ustage(ustage);
    ss.add_ustage(ustage_dup);
    ss.solver("fatrop", {{"jit", true}});
    std::cout <<cs::MX::evalf(cs::MX::substitute(sumsqr(x), x, cs::MX::ones(10))) << std::endl;
    auto test_func = cs::Function("test", {x}, {x});
    std::cout << test_func(cs::DM::ones(10)) << std::endl;
    auto funcc = ss.to_function("example_dense", {}, {ustage.at_t0(x)});
    std::cout << funcc(std::vector<cs::DM>{}) << std::endl;
    // std::cout << funcc(std::vector<cs::DM>{}) << std::endl;
  }
  // {
  //   auto func = fatropy::get_func_from_py("example", "get_func");
  //   func(std::vector<cs::DM>{2.5});
  // }

  return 0;
}