#include "ocp/BFOCPBasic.hpp"
#include "ocp/BFOCP.hpp"
#include "ocp/BFOCPAdapter.hpp"
#include "aux/SmartPtr.hpp"
#include "ocp/FatropOCP.hpp"
#include "ocp/OCPLinearSolver.hpp"
#include "ocp/OCPLSRiccati.hpp"
#include "aux/FatropVector.hpp"
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "ocp/OCPNoScaling.hpp"
#include "solver/FatropParams.hpp"
#include "solver/FatropAlg.hpp"
#include "sparse/SparseOCP.hpp"

using namespace fatrop;
int main()
{
    // RefCountPtr<BFOCP> ocptemplatebasic =
    //     new BFOCPBasic(BFOCPBasic::from_shared_lib("./../../OCP_specification/f.so", 3));
    // RefCountPtr<OCP> ocptempladapter = new BFOCPAdapter(ocptemplatebasic);
    // RefCountPtr<OCPLinearSolver> ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims());
    // RefCountPtr<FatropParams> params = new FatropParams();
    // RefCountPtr<OCPScalingMethod> ocpscaler = new OCPNoScaling(params);
    // FatropOCP ocpalg(ocptempladapter, ocplsriccati, ocpscaler);
    // NLPDims nlpdims = ocpalg.GetNLPDims();
    // // blasfeo_timer timer;
    // FatropMemoryVecBF ux(nlpdims.nvars, 5);
    // FatropMemoryVecBF lags(nlpdims.neqs, 5);
    // // int N = 1000;
    // FatropVecBF scalesx = ux[1];
    // FatropVecBF dux = ux[3];
    // FatropVecBF grad = ux[2];
    // FatropVecBF cv = lags[2];
    // FatropVecBF scaleslam = lags[1];
    // FatropVecBF dlam = lags[3];
    // FatropVecBF duinf = ux[4];

    // // FatropAlg fatropalg(ocpalg)
    // double obj_sc = 1.0;
    // ocpalg.ComputeScalings(
    //     obj_sc,
    //     scalesx,
    //     scaleslam,
    //     scaleslam);
    // ocpalg.EvalConstraintViolation(ux[0], ux[0], cv);
    // ocpalg.EvalGrad(1.0, ux[0], grad);
    // ocpalg.EvalHess(1.0, ux[0], lags[0]);
    // ocpalg.EvalJac(ux[0], ux[0]);
    // // grad.print();
    // ocpalg.OCPInitializer_.AdaptKKTInitial(&ocpalg.ocpkktmemory_, grad, grad);
    // ocpalg.ComputeSD(0.0, 0.0, 0.0,0.0, dux, dlam, dlam, dlam, dlam, dlam, dlam, dlam, dlam, dlam);
    // cout << "dlam" << endl;
    // dlam.print();
    // // grad.print();
    // // ocpalg.duinfevaluator_.DuInfEval(&ocpalg.ocpkktmemory_, 1.0, dlam, grad, duinf);

    // // cout << "du inf Linf" << endl;
    // // cout << Linf(duinf) << endl;

    // // }
    // // double el = blasfeo_toc(&timer);
    // // cout << "time elapsed " << el / N << endl;
    // RefCountPtr<Sparse_OCP> ocplssparse = new Sparse_OCP(ocptempladapter->GetOCPDims(), ocpalg.ocpkktmemory_);
    // ocplssparse->computeSD(&ocpalg.ocpkktmemory_, 0.0, 0.0, 0.0,0.0, dux, dlam, dlam, dlam, dlam, dlam, dlam, dlam, dlam, dlam);
    // cout << "dlam" << endl;
    // dlam.print();
    // // ocpalg.duinfevaluator_.DuInfEval(&ocpalg.ocpkktmemory_, 1.0, dlam, grad, duinf);
    // // cout << "du inf Linf" << endl;
    // cout << Linf(duinf) << endl;
    // // ocplssparse->print();
}