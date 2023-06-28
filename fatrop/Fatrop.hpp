#ifndef FATROP_INCLUDED
#define FATROP_INCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "ocp/OCPKKT.hpp"
#include "ocp/OCPAdapter.hpp"
#include "ocp/OCPAbstact.hpp"
#include "auxiliary/FatropVector.hpp"
#include "solver/FatropAlg.hpp"
#include "solver/FatropData.hpp"
#include "ocp/OCPScalingMethod.hpp"
#include "ocp/OCPNoScaling.hpp"
#include "solver/AlgStrategy.hpp"
#include "solver/FatropOptions.hpp"
#include "function_evaluation/CasadiCodegen.hpp"
#include "solver/AlgBuilder.hpp"
#include "ocp/StageOCPApplication.hpp"
#include "ocp/FatropOCPBuilder.hpp"

// #include "SparseSolvers/InterfaceMUMPS.hpp"
#endif //FATROP_INCLUDED