#include "OCP/OCPTemplateBasic.hpp"
#include "OCP/OCPTemplate.hpp"
using namespace fatrop;
int main()
{
    OCPTemplateBasic ocptemplatebasic = OCPTemplateBasic::from_shared_lib("./f.so", 10);
}