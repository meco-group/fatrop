#include <iostream>
#include <fstream>
#include "ocp/OCPBuilder.hpp"
#include <string>
using namespace fatrop;
int main()
{
    OCPBuilder ocpbuilder("./../../specification_examples/Robot2/f.so", "./../../specification_examples/Robot2/test.json");
    // vector<double> initial_x = result["initial_x"].get_number_array<double>("%lf");
    // cout << "end" << endl;
}