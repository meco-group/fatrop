#include <iostream>
#include <fstream>
#include "ocp/OCPBuilder.hpp"
#include <string>
#include <ocp/BasicOCPApplication.hpp>
using namespace fatrop;
int main(int argc, char **argv)
{
    // if (argc == 3)
    // {
        // OCPBuilder ocpbuilder(argv[1], argv[2]);
        shared_ptr<BasicOCPApplication> app = BasicOCPApplicationBuilder::FromRockitInterface("./casadi_codegen.so", "./casadi_codegen.json");
        // shared_ptr<OCPBuilder> app =make_shared<OCPBuilder>("./casadi_codegen.so", "./casadi_codegen.json");
        cout << "Building the solver " << endl;
        app -> Build();
        // shared_ptr<FatropApplication> solver = ocpbuilder.Build();
        // usage of parameter setter
        // ocpbuilder.GetParameterSetter("target_pos")->SetValue({1., 1., 1.});
        // cout << "Calling the solver " << endl;
        app-> Optimize();
        // usage of parameter sampler
        // auto res = app-> GetSampler("Frame") -> Sample(app-> LastSolution(), app->GlobalParameters(), app->StageParameters());
    // }
    // else
    // {
    //     cout << "run me as RunFatrop f.so f.json" << endl;
    // }
}