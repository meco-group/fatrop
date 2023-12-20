#include "casadi_jit.hpp"
namespace fatrop
{
    namespace spectool
    {
        CasadiFEJit::CasadiFEJit(const cs::Function &F, const cs::Dict &jit_options_)
        {
            cs::FunctionInternal *func_internal = F.get();
            if (func_internal->has_codegen())
            {
                jit_name_ = F.name();
                jit_name_ = casadi::temporary_file(jit_name_, ".c");
                jit_name_ = std::string(jit_name_.begin(), jit_name_.begin() + jit_name_.size() - 2);
                // JIT everything
                casadi::Dict opts;
                // Override the default to avoid random strings in the generated code
                opts["prefix"] = "jit";
                casadi::CodeGenerator gen(jit_name_, opts);
                gen.add(F);
                // jit_options_ = casadi::Dict({{"flags", "-Ofast -march=native -ffast-math"}});
                jit_directory = casadi::get_from_dict(jit_options_, "directory", std::string(""));
                std::string compiler_plugin_ = "shell";
                compiler_ = casadi::Importer(gen.generate(jit_directory), compiler_plugin_, jit_options_);
                eval_ = (eval_t)compiler_.get_function(F.name());
                // cache the function
                compiled_jit = true;
                casadi_assert(eval_ != nullptr, "Cannot load JIT'ed function.");
            }
            else
            {
                std::cout << "jit compilation not possible for the provided functions" << std::endl;
            }
        }
        CasadiFEJit::~CasadiFEJit()
        {
            std::cout << "destructor called" << std::endl;
            if (compiled_jit)
            {
                std::string jit_directory = get_from_dict(jit_options_, "directory", std::string(""));
                std::string jit_name = jit_directory + jit_name_ + ".c";
                if (remove(jit_name.c_str()))
                    casadi_warning("Failed to remove " + jit_name);
            }
        }

    }
}