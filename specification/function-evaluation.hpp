#pragma once
#include <casadi/casadi.hpp>
#include "casadi/core/code_generator.hpp"
#include "casadi/core/importer.hpp"
#include "casadi/core/function_internal.hpp"
#include <array>
#define JIT_HACKED_CASADI 1
namespace fatrop
{
    namespace specification
    {
        typedef int (*eval_t)(const double **arg, double **res,
                              long long int *iw, double *w, int);

        class eval_bf
        {
        public:
            eval_bf(){};
            eval_bf(const casadi::Function &funcin)
            {
                m = (int)funcin.size1_out(0);
                n = (int)funcin.size2_out(0);
                func = funcin;
                mem = 0;
                // allocate work vectors
                size_t sz_arg,
                    sz_res, sz_iw, sz_w;
                sz_arg = func.n_in();
                sz_res = func.n_out();
                func.sz_work(sz_arg, sz_res, sz_iw, sz_w);
                iw.resize(sz_iw);
                w.resize(sz_w);
                bufout.resize(func.nnz_out(0));
                bufdata.resize(sz_res);
                resdata.resize(sz_res);
                argdata.resize(sz_arg);
                n_in = func.n_in();
                fast_jit();

                // resdata = {nullptr};
                dirty = false;
            }
            void fast_jit()
            {
                casadi::FunctionInternal *func_internal = func.get();
                jit_name_ = func.name();
                jit_name_ = casadi::temporary_file(jit_name_, ".c");
                jit_name_ = std::string(jit_name_.begin(), jit_name_.begin() + jit_name_.size() - 2);
                if (func_internal->has_codegen())
                {
                    // this part is based on casadi/core/function_internal.cpp -- all rights reserved to the original authors
                    // JIT everything
                    casadi::Dict opts;
                    // Override the default to avoid random strings in the generated code
                    opts["prefix"] = "jit";
                    casadi::CodeGenerator gen(jit_name_, opts);
                    gen.add(func);
                    jit_options_ = casadi::Dict({{"flags", "-Ofast -march=native -ffast-math"}});
                    jit_directory = casadi::get_from_dict(jit_options_, "directory", std::string(""));
                    std::string compiler_plugin_ = "shell";
                    compiler_ = casadi::Importer(gen.generate(jit_directory), compiler_plugin_, jit_options_);
                    eval_ = (eval_t)compiler_.get_function(func.name());
                    casadi_assert(eval_ != nullptr, "Cannot load JIT'ed function.");
                    compiled_jit = true;
                }
                else
                {
                    std::cout << "jit compilation not possible for the provided functions" << std::endl;
                }
            }
            std::string jit_name_;
            casadi::Dict jit_options_;
            std::string jit_directory;
            bool compiled_jit = false;
            void operator=(const eval_bf &other)
            {
                // copy all member values
                m = other.m;
                n = other.n;
                n_in = other.n_in;
                func = other.func;
                mem = other.mem;
                bufout = other.bufout;
                bufdata = other.bufdata;
                resdata = other.resdata;
                argdata = other.argdata;
                iw = other.iw;
                w = other.w;
                dirty = other.dirty;
                fast_jit();
                // eval_ = other.eval_;
                // other.eval_ = nullptr;
            }
            void operator()(const double ** arg, MAT *res, double *buff)
            {
                if (dirty)
                    return;
                // inputs
                for (int j = 0; j < n_in; j++)
                    argdata[j] = arg[j];
                // outputs
                bufdata[0] = buff;
#ifdef JIT_HACKED_CASADI
                eval_(argdata.data(), bufdata.data(), iw.data(), w.data(), 0);
#else
                func(argdata.data(), bufdata.data(), iw.data(), w.data(), 0);
#endif
                // func(arg, bufdata);
                PACKMAT(m, n, buff, m, res, 0, 0);
            }
            // copy operator
            void operator()(const double** arg, MAT *res)
            {
                this->operator()(arg, res, bufout.data());
            }
            void operator()(const double ** arg, double *res)
            {
                if (dirty)
                    return;
                // inputs
                for (int j = 0; j < n_in; j++)
                    argdata[j] = arg[j];
                // outputs
                resdata[0] = res;
#ifdef JIT_HACKED_CASADI
                eval_(argdata.data(), resdata.data(), iw.data(), w.data(), 0);
#else
                func(argdata.data(), resdata.data(), iw.data(), w.data(), 0);
#endif
                // func(arg, resdata);
            }
            ~eval_bf()
            {
                if (compiled_jit)
                {
                    std::string jit_directory = get_from_dict(jit_options_, "directory", std::string(""));
                    std::string jit_name = jit_directory + jit_name_ + ".c";
                    if (remove(jit_name.c_str()))
                        casadi_warning("Failed to remove " + jit_name);
                }
            }
            int m;
            int n;
            int mem;
            int n_in;
            eval_t eval_;
            casadi::Importer compiler_;
            std::vector<double> bufout;
            std::vector<double *> bufdata;
            std::vector<double *> resdata;
            std::vector<const double *> argdata;
            std::vector<long long int> iw;
            std::vector<double> w;
            casadi::Function func;
            bool dirty = true;
            bool bfgs = true;
        };
    };
};