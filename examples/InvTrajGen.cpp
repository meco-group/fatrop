#include <iostream>
#include <fstream>
#include "ocp/OCPBuilder.hpp"
#include <string>
using namespace fatrop;
int main()
{
    OCPBuilder ocpbuilder("./../../ocp_specification/f.so", "./../../ocp_specification/test.json");
    // vector<double> initial_x = result["initial_x"].get_number_array<double>("%lf");
    // cout << "end" << endl;
}