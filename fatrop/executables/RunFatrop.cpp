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
        app->Optimize();
        vector<double> result = (app->LastBasicOCPSolution()).Eval(app->GetEvaluator("control_u")->at_t0());
    }
    else
    {
        cout << "run me as RunFatrop f.so f.json" << endl;
    }
}