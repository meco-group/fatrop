#include <iostream>
#include <fstream>
#include "ocp/OCPBuilder.hpp"
#include <string>
using namespace fatrop;

int main(int argc, char **argv)
{
    if (argc == 3)
    {
        OCPBuilder ocpbuilder(argv[1], argv[2]);
        ocpbuilder.fatropalg->Optimize();
    }
    else
    {
        cout << "run me as ./RunFatrop f.so f.json" << endl;
    }
}