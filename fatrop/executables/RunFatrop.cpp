#include <iostream>
#include <fstream>
#include "ocp/OCPBuilder.hpp"
#include <string>
using namespace fatrop;
int main(int argc, char **argv)
{
    // if (argc == 3)
    // {
        // OCPBuilder ocpbuilder(argv[1], argv[2]);
        OCPBuilder ocpbuilder("./casadi_codegen.so", "./casadi_codegen.json");
        shared_ptr<FatropApplication> solver = ocpbuilder.Build();
        // usage of parameter setter
        // ocpbuilder.GetParameterSetter("target_pos")->SetValue({1., 1., 1.});
        solver->Optimize();
        ocpbuilder.SampleMaxEnt(100.0);
        // usage of parameter sampler
        // auto res = ocpbuilder.GetSampler("state_pos") -> Sample();
    // }
    // else
    // {
    //     cout << "run me as RunFatrop f.so f.json" << endl;
    // }
}