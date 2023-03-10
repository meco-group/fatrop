#ifndef BASICOCPAPPLICATIONINCLUDED
#define BASICOCPAPPLICATIONINCLUDED
#include "ocp/BFOCPBasic.hpp"
#include "solver/AlgBuilder.hpp"
#include "ocp/BFOCPAdapter.hpp"
#include "solver/FatropAlg.hpp"
#include "ocp/FatropOCP.hpp"
#include "ocp/OCPLSRiccati.hpp"
#include "ocp/OCPNoScaling.hpp"
#include "ocp/FatropOCPBuilder.hpp"
#include <map>
#include "json/json.h"
#include <fstream>
#include <sstream>
#include <string>

namespace fatrop
{
    class BasicOCPApplication
    {
    public:
        BasicOCPApplication(const shared_ptr<BasicOCP> &ocp) : ocp_(ocp)
        {
        }
        void Build()
        {
            fatropparams_ = make_shared<FatropParams>();
            shared_ptr<FatropNLP> nlp(FatropOCPBuilder(ocp_, fatropparams_).Build());
            AlgBuilder algbuilder;
            algbuilder.BuildFatropAlgObjects(nlp, fatropparams_, fatropdata_, journaller_);
            fatropalg_ = algbuilder.BuildAlgorithm();
            //    algbuilder.BuildFatropAlgObjects(nlp, fatropdata, fatropparams, filter, linesearch, journaller);
            dirty = false;
        }
        int Optimize()
        {
            cout << "in optimize " << endl;
            assert(!dirty);
            cout << "dirty = " << dirty << endl;
            return fatropalg_->Optimize();
        }

        FatropVecBF &LastSolution()
        {
            assert(!dirty);
            return fatropdata_->x_curr;
        }

    private:
        bool dirty = true;
        const shared_ptr<BasicOCP> ocp_;
        shared_ptr<FatropData> fatropdata_;
        shared_ptr<Journaller> journaller_;
        shared_ptr<FatropAlg> fatropalg_;
        shared_ptr<FatropParams> fatropparams_;
        // map with all available samplers
    };

    class BasicOCPApplicationBuilder
    {
    public:
        static shared_ptr<BasicOCPApplication> FromRockitInterface(const string &functions, const string &json_spec_file)
        {
            shared_ptr<DLHandler> handle = make_shared<DLHandler>(functions);
            std::ifstream t(json_spec_file);
            std::stringstream buffer;
            buffer << t.rdbuf();
            json::jobject json_spec = json::jobject::parse(buffer.str());
            auto ocptemplatebasic = BasicOCPBuilder::FromRockitInterface(handle, json_spec);
            // instantiate the BasicOCPApplication
            auto result = make_shared<BasicOCPApplication>(ocptemplatebasic);
            // add all samplers
            // add all parameter setters
            return result;
        }
    };
}
#endif // BASICOCPAPPLICATIONINCLUDED