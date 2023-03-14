#include <iostream>
#include <fstream>
#include <string>
#include <ocp/BasicOCPApplication.hpp>
using namespace fatrop;
int main(int argc, char **argv)
{

    if (argc == 3)
    {
        shared_ptr<BasicOCPApplication> app = BasicOCPApplicationBuilder::FromRockitInterface(argv[1], argv[2]);
        app -> SetNumericOption("tol", 1e-6);
        app -> Optimize();
        vector<double> result = (app->LastBasicOCPSolution()).Eval(app->GetExprEvaluator("control_u")->at_t0());
    }
    else
    {
        cout << "run me as RunFatrop f.so f.json" << endl;
    }
}