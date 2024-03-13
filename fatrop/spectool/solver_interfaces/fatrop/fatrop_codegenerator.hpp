#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <set>
#include "fatrop_codegenerator_template.hpp"
namespace fatrop
{
    namespace spectool
    {
        struct CodegenData
        {
            std::vector<int> nu;
            std::vector<int> nx;
            std::vector<int> np_stage;
            int np_global;
            std::vector<int> ng_eq;
            std::vector<int> ng_ineq;
            std::vector<int> nxp1;
            std::vector<std::string> eval_BAbt;
            std::vector<std::string> eval_RSQrqt;
            std::vector<std::string> eval_Ggt;
            std::vector<std::string> eval_Ggt_ineq;
            std::vector<std::string> eval_b;
            std::vector<std::string> eval_g;
            std::vector<std::string> eval_g_ineq;
            std::vector<std::string> eval_rq;
            std::vector<std::string> eval_L;
            std::vector<double> Lb;
            std::vector<double> Ub;
        };
        class Codegenerator
        {
            public:
            template <typename T>
            void add(const std::string &pre, const std::vector<T> &v, const std::string &post)
            {
                for (auto &el : v)
                {
                    res += pre + std::to_string(el) + post;
                }
            }
            void add(const std::string &s)
            {
                res += s;
            };
            void add_line(const std::string &s)
            {
                return add(s + "\n");
            }
            std::string get()
            {
                return res;
            }
            std::string res;
        };
        class FatropCodegenerator
        {
            public:
            static std::string generate_inst_eval(const std::string &func_name)
            {
                return "EvalCasGen(" + func_name + "_incref," + func_name + "_decref," + func_name + "_checkout," + func_name + "_release," + func_name + "_n_in_fcn," + func_name + "_n_out_fcn," + func_name + "_sp_in," + func_name + "_sp_out," + func_name + "_work_t work," + func_name + "_eval)";
            }
            static std::string generate_ustate_eval(const CodegenData &data)
            {
            
            }
            static std::string generate(const CodegenData &data)
            {
                Codegenerator gen;
                // all the following things should be part of one memory vector
                // instantiate all EvalCasGen's
                // instantiate all UStageEval's shared_ptrs
                // instantiate vector with UstageEvals shared ptrs
                gen.add(pre);
                gen.add(post);
            }
        };
    } // namespace spectool
} // namespace fatrop