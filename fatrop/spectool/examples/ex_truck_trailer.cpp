#include "fatrop/spectool/spectool.hpp"
#include "fatrop/fatrop.hpp"
#include <casadi/casadi.hpp>
using namespace std;
using namespace fatrop;
using namespace fatrop::spectool;

int main()
{
    auto ocp = Ocp();
    // # truck:
    double L0 = 0.4;
    double M0 = 0.1;
    double W0 = 0.2;

    // # trailer1:
    double L1 = 1.1;
    double M1 = 0.2;
    double W1 = 0.2;

    // # trailer2:
    double L2 = 0.8;
    double M2 = 0.1;
    double W2 = 0.2;
    double pi = 3.14159265358979323846;

    // # Trailer model
    auto theta2 = ocp.state();
    auto x2     = ocp.state();
    auto y2     = ocp.state();

    auto theta1 = ocp.state();
    auto x1     = x2 + L2*cos(theta2) + M1*cos(theta1);
    auto y1     = y2 + L2*sin(theta2) + M1*sin(theta1);

    auto theta0 = ocp.state();
    auto x0     = x1 + L1*cos(theta1) + M0*cos(theta0);
    auto y0     = y1 + L1*sin(theta1) + M0*sin(theta0);

    auto delta0 = ocp.state();
    auto ddelta0 = ocp.control();
    auto v0     = ocp.state();
    auto dv0     = ocp.control();

    auto beta01 = theta0 - theta1;
    auto beta12 = theta1 - theta2;

    auto dtheta0 = v0/L0*tan(delta0);
    auto dtheta1 = v0/L1*sin(beta01) - M0/L1*cos(beta01)*dtheta0;
    auto v1 = v0*cos(beta01) + M0*sin(beta01)*dtheta0;

    auto dtheta2 = v1/L2*sin(beta12) - M1/L2*cos(beta12)*dtheta1;
    auto v2 = v1*cos(beta12) + M1*sin(beta12)*dtheta1;
    auto T = ocp.state();
    auto d_theta2 =  dtheta2;
    auto d_x2 =      v2*cos(theta2);
    auto d_y2 =      v2*sin(theta2);
    auto d_theta1 =  dtheta1;
    auto d_theta0 =  dtheta0;
    auto d_delta0 =  ddelta0;
    auto d_v0 =      dv0;
    auto intg = IntegratorRk4({
    {theta2 ,  d_theta2},
    {x2 ,      d_x2},
    {y2 ,      d_y2},
    {theta1 ,  d_theta1},
    {theta0 ,  d_theta0},
    {delta0 ,  d_delta0},
    {v0 ,      d_v0}
    }, T/50.);
    auto stages = std::vector<Stage>(); 
    for(int i =0; i<1; i++)
    {
        stages.push_back(ocp.new_stage(50));
        stages.back().set_next(theta2, intg(theta2));
        stages.back().set_next(x2, intg(x2));
        stages.back().set_next(y2, intg(y2));
        stages.back().set_next(theta1, intg(theta1));
        stages.back().set_next(theta0, intg(theta0));
        stages.back().set_next(delta0, intg(delta0));
        stages.back().set_next(v0, intg(v0));
        stages.back().set_next(T, T);
        stages.back().set_next(T, T);
        double speedf = 1;
        stages.back().subject_to(-.2 * speedf <= (v0 <= .2* speedf), at::t0, at::mid, at::tf);
        stages.back().subject_to(-1 <= (dv0 <= 1), at::t0, at::mid);
        stages.back().subject_to(-pi/6 <= (delta0 <= pi/6), at::t0, at::mid, at::tf);
        stages.back().subject_to(-pi/10<= (ddelta0 <= pi/10), at::t0, at::mid);
        stages.back().subject_to(-pi/2 <= (beta01 <= pi/2), at::t0, at::mid, at::tf);
        stages.back().subject_to(-pi/2 <= (beta12 <= pi/2), at::t0, at::mid, at::tf);
        stages.back().add_objective(beta01*beta01, at::t0, at::mid, at::tf);
        stages.back().add_objective(beta12*beta12, at::t0, at::mid, at::tf);
    }
    
    // # parameters
    double x2_t0 = 0.;
    double y2_t0 = 0.;
    double theta2_t0 = 0.;
    double theta1_t0 = 0.;
    double theta0_t0 = 0.;
    double x2_tf = 0.;
    double y2_tf = -2.;
    double theta2_tf = 2*pi/4;
    double theta1_tf = 2*pi/4;
    double theta0_tf = 2*pi/4;
    // # Initial constraints
    ocp.at_t0().subject_to(x2 == x2_t0);
    ocp.at_t0().subject_to(y2 == y2_t0);
    ocp.at_t0().subject_to(theta2 == theta2_t0);
    ocp.at_t0().subject_to(theta1 == theta1_t0);
    ocp.at_t0().subject_to(theta0 == theta0_t0);
    // # Final constraint
    ocp.at_tf().subject_to(x2 == x2_tf);
    ocp.at_tf().subject_to(y2 == y2_tf);
    ocp.at_tf().subject_to(theta2 == 2*pi/4);
    ocp.at_tf().subject_to(beta12 == theta1_tf - theta2_tf);
    ocp.at_tf().subject_to(beta01 == theta0_tf - theta1_tf);

    // # objective
    ocp.at_tf().add_objective(T);
    
    ocp.set_initial(T, cs::MX(20.));
    ocp.set_initial(theta0, cs::MX(.1));
    ocp.set_initial(theta1, cs::MX(0));
    ocp.set_initial(theta2, cs::MX(0));
    ocp.set_initial(v0,     cs::MX(-.2));

    auto a = ocp.parameter();

    ocp.solver("fatrop", casadi::Dict({{"post_expand", true}});
    auto func = ocp.to_function({a}, {ocp.at_tf(T)}));
    func(cs::DM(0.0));
    return 0;
}