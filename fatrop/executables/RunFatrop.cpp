#include <iostream>
#include <fstream>
#include <string>
#include <ocp/StageOCPApplication.hpp>
using namespace fatrop;
int main(int argc, char **argv)
{
    if (argc == 3)
    {
        //// dynamic memory allocation  
        shared_ptr<StageOCPApplication> app = StageOCPApplicationBuilder::FromRockitInterface(argv[1], argv[2]);
        auto eval_expression = app->GetExpression("control_u")->at_t0();
        vector<double> u0_result(eval_expression -> Size());
        ///  no dynamic memory allocation
        app->Optimize();
        // ///  retrieve solution
        app->LastStageOCPSolution().Eval(eval_expression, u0_result);
        ///  initialize next solver run 
        app->SetInitial(app->LastStageOCPSolution());
        ///  change solver options
        app->SetOption("tol", 1e-6);
        app->SetOption("mu_init", 1e-5);
        app->SetOption("bound_push", 1e-7); // all *_bound_push variables
        app->SetOption("warm_start_mult_bound_push", 1e-7); 
        app->SetOption("accept_every_trial_step", false);
        app->SetOption("iterative_refinement", false); // fast_step_computation
        app->SetOption("warm_start_init_point", true);
        // cout << app -> GetOptions();
        app->Optimize();
    }
    else
    {
        cout << "run me as RunFatrop f.so f.json" << endl;
    }
}