#include <iostream>
#include <fstream>
// #include "ocp/OCPBuilder.hpp"
#include <string>
#include <ocp/BasicOCPApplication.hpp>
using namespace fatrop;
int main(int argc, char **argv)
{
    if (argc == 3)
    {
        shared_ptr<BasicOCPApplication> app = BasicOCPApplicationBuilder::FromRockitInterface(argv[1], argv[2]);
        app -> Optimize();
    }
    else
    {
        cout << "run me as RunFatrop f.so f.json" << endl;
    }
}