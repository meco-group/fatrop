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
        OCPBuilderAL ocpbuilderal("../robot.so", "../robot.json");
        ocpbuilderal.GN = false;
        ocpbuilderal.DDP = false;
        shared_ptr<FatropApplication> solveral = ocpbuilderal.Build();
        OCPBuilder ocpbuilder("../robot.so", "../robot.json");
        ocpbuilder.GN = false;
        ocpbuilder.DDP = false;
        shared_ptr<FatropApplication> solver = ocpbuilder.Build();
        vector<double> result = vector<double>(ocpbuilder.fatropocp->GetNLPDims().nvars);
        solveral->Optimize();
        solveral->GetSolution(result);
        solver->SetInitial(result);
        solver->Optimize();
    }
    else
    {
        cout << "run me as RunFatrop f.so f.json" << endl;
    }
}
