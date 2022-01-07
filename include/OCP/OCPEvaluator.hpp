#ifndef OCPEVALUATORINCLUDED
#define OCPEVALUATORINCLUDED
#include "FatropOCPKKT.hpp"
#include "FUNCTION_EVALUATION/FunctionEvaluation.hpp"
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
namespace fatrop
{
    class OCP_evaluator
    {
        int eval_hess(OCP_KKT *OCP, const fatrop_vector_bf &ux, const fatrop_vector_bf &lam)
        {
            OCPMACRO(MAT *, RSQrqt, _p);
            return 0;
        }
        int eval_jac(OCP_KKT *OCP, const fatrop_vector_bf &ux, const fatrop_vector_bf &lam)
        {
            OCPMACRO(MAT *, RSQrqt, _p);
            OCPMACRO(MAT *, Ggt, _p);
            return 0;
        }

    public:
        vector<fatrop_eval_base *> BAbtf;
        vector<fatrop_eval_base *> RSQrqtf;
        vector<fatrop_eval_base *> Ggtf;
    };
} // namespace fatrop

#endif // OCPEVALUATORINCLUDED