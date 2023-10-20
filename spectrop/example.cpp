#include "spectrop.hpp"
using namespace std;
using namespace fatrop;
using namespace fatrop::spectrop;
int main()
{
    auto ocp = Ocp();
    auto x = ocp.state();
    auto u = ocp.control();
    auto p = ocp.parameter();

    // /*
    //   =----  initial stage ----=
    // */
    auto initial_stage = ocp.new_stage(); // states and controls are derived automatically
    /* constraint   */ initial_stage.subject_to(x == 0);
    /* dynamics     */ initial_stage.set_next(x, x + u + 1);
    /* objective    */ initial_stage.add_objective(u*u);

    // /*
    //   =----  middle stage ----=
    // */
    auto middle_stage = ocp.new_stage(18); // 18 is the number of nodes // the states, controls and parameters are derived automatically
    /* constraints  */ middle_stage.subject_to(x <= 1);
    /* objective    */ middle_stage.add_objective(u * u);
    /* dynamics     */ middle_stage.set_next(initial_stage.dynamics());

    // /*
    //   =----  terminal stage ----=
    // */
    auto terminal_stage = ocp.new_stage(); // the last stage also has x2 as state
    auto x2 = ocp.state();
    /* objective    */ terminal_stage.add_objective(x2 * x2);
    /* dynamics     */ middle_stage.set_next(x2, x + u + 1);

    // ocp.sample(x);           // check which stages have x as state, evaluate and concatenate in a matrix
    // ocp.set_initial(x, 0.5); // set initial value of x for all stages that have x as state
    return 0;
}