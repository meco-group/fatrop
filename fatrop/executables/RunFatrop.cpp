#include <iostream>
#include <fstream>
#include <string>
#include <ocp/StageOCPApplication.hpp>
using namespace fatrop;
using namespace std;
int main(int argc, char **argv)
{
    if (argc == 3)
    {
        //// dynamic memory allocation  
        StageOCPApplication app = StageOCPApplicationFactory::from_rockit_interface(argv[1], argv[2]);
        auto eval_expression = app.get_expression("control_u").at_t0();
        vector<double> u0_result(eval_expression.size());
        ///  no dynamic memory allocation
        app.optimize();
        // ///  retrieve solution
        app.last_stageocp_solution().evaluate(eval_expression, u0_result);
        ///  initialize next solver run 
        app.set_initial(app.last_solution());
        ///  change solver options
        app.set_option("tol", 1e-6);
        app.set_option("mu_init", 1e-5);
        app.set_option("bound_push", 1e-7); // all *_bound_push variables
        app.set_option("warm_start_mult_bound_push", 1e-7); 
        app.set_option("accept_every_trial_step", false);
        app.set_option("iterative_refinement", false); // fast_step_computation
        app.set_option("warm_start_init_point", true);
        app.set_option("print_level", 0);
        // cout << app -> GetOptions();
        app.optimize();
    }
    else
    {
        cout << "run me as RunFatrop f.so f.json" << endl;
    }
}