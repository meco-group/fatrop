#include <iostream>
#include <fstream>
#include "ocp/OCPBuilder.hpp"
#include <string>
using namespace fatrop;
int main(int argc, char **argv)
{
    // if (argc == 3)
    // {
    // string funcs = "/home/lander/projects/fatrop_motivation_examples/foobar/casadi_codegen.so";
    // string jsons = "/home/lander/projects/fatrop_motivation_examples/foobar/casadi_codegen.json";
    string funcs = "/home/lander/projects/fatrop-invariants-tutorial/foobar/casadi_codegen.so";
    string jsons = "/home/lander/projects/fatrop-invariants-tutorial/foobar/casadi_codegen.json";
    // OCPBuilder ocpbuilder(argv[1], argv[2]);
    OCPBuilder ocpbuilder(funcs, jsons);
    shared_ptr<FatropApplication> solver = ocpbuilder.Build();
    solver->Optimize();
    double target[] = {-1., 1., 1.};
    ocpbuilder.GetParameterSetter("target_pos")->SetValue(target);
    solver->Optimize();
    auto sampler = ocpbuilder.GetSampler("state_x2");
    vector<double> res(sampler -> Size());
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    sampler -> Sample(res);
    double el = blasfeo_toc(&timer);
    for (auto resi : res)
    {
        cout << resi << endl;
    }
    cout << "sample time " << el << endl;
    // }
    // else
    // {
    //     cout << "run me as RunFatrop f.so f.json" << endl;
    // }
}
