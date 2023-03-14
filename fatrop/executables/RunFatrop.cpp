#include <iostream>
#include <fstream>
#include <string>
#include <ocp/BasicOCPApplication.hpp>
using namespace fatrop;
int main(int argc, char **argv)
{

    class CustoMExpression : public StageExpression
    {
        void Eval(const double *u, const double *x, const double *global_params, const double *stage_params, double *res)
        {
            res[0] = u[0];
            res[1] = u[0] * cos(x[0]);
        };
        int n_rows()
        {
            return 2;
        }

        int n_cols()
        {
            return 1;
        }
    };

    if (argc == 3)
    {
        shared_ptr<BasicOCPApplication> app = BasicOCPApplicationBuilder::FromRockitInterface(argv[1], argv[2]);
        auto custom_evaluator = app->GetExprEvaluator(make_shared<CustoMExpression>())->at_control();
        app->SetNumericOption("tol", 1e-6);
        app->Optimize();
        vector<double> result = (app->LastBasicOCPSolution()).Eval(app->GetExprEvaluator("control_u")->at_t0());
        vector<double> result_custom = (app->LastBasicOCPSolution()).Eval(custom_evaluator);
    }
    else
    {
        cout << "run me as RunFatrop f.so f.json" << endl;
    }
}