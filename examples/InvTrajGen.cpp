#include <iostream>
#include <fstream>
#include "ocp/OCPBuilder.hpp"
#include <string>
using namespace fatrop;
int main()
{
    OCPBuilder ocpbuilder("./../../specification_examples/Rocket/f.so", "./../../specification_examples/Rocket/test.json");
    ocpbuilder.fatropalg->Optimize();
    // vector<double> initial_x = result["initial_x"].get_number_array<double>("%lf");
    // cout << "end" << endl;
}