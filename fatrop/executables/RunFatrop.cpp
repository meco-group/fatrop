#include <iostream>
#include <fstream>
#include "ocp/OCPBuilderAL.hpp"
#include "ocp/OCPBuilder.hpp"
#include <string>
using namespace fatrop;

int main(int argc, char **argv)
{
    // if (argc == 3)
    if (true)
    {

#ifdef ENABLE_MULTITHREADING
        cout << "Multithreading enabled" << endl;
#else
        cout << "Multithreading disabled" << endl;
#endif

        // OCPBuilder ocpbuilder(argv[1], argv[2]);
        // OCPBuilder ocpbuilder("../robot.so","../robot.json");
        // OCPBuilderAL ocpbuilderal("../Rocket_example.so", "../Rocket_example.json");
        OCPBuilder ocpbuilder("../casadi_codegen.so", "../casadi_codegen.json");
        //OCPBuilder ocpbuilder("../Rocket_example.so", "../Rocket_example.json");
        // OCPBuilderAL ocpbuilderal("../invariants.so", "../invariants.json");
        // OCPBuilder ocpbuilder("../invariants.so", "../invariants.json");
        // ocpbuilderal.GN = false;
        // ocpbuilderal.DDP = false;
        // shared_ptr<FatropApplication> solveral = ocpbuilderal.Build();
        // ocpbuilder.GN = true;
        // ocpbuilder.DDP = true;
        shared_ptr<FatropApplication> solver = ocpbuilder.Build();
        // vector<double> result = vector<double>(ocpbuilder.fatropocp->GetNLPDims().nvars);
        // solveral->Optimize();
        // solveral->GetSolution(result);
        // solver->SetInitial(result);
        solver->Optimize();
    }
    else
    {
        cout << "run me as RunFatrop f.so f.json" << endl;
    }
}
