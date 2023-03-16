#include <iostream>
#include <fstream>
#include <string>
#include <ocp/StageOCPApplication.hpp>
using namespace fatrop;
int main(int argc, char **argv)
{
    if (argc == 3)
    {
        shared_ptr<StageOCPApplication> app = StageOCPApplicationBuilder::FromRockitInterface(argv[1], argv[2]);
        app->SetOption("tol", 1e-6);
        app->Optimize();
        vector<double> result = (app->LastStageOCPSolution()).Eval(app->GetExprEvaluator("control_u")->at_t0());
        app->SetInitial(app->LastStageOCPSolution());
        // app -> SetOption("warm_start_init_point", true);
        app -> Optimize();
    }
    else
    {
        cout << "run me as RunFatrop f.so f.json" << endl;
    }
}