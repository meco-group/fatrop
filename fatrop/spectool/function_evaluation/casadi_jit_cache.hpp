#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "casadi_jit.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        class CasadiJitCache: public std::map<size_t, std::shared_ptr<CasadiFEJit>>
        {
            void save(const std::string& filename)
            {
                // iterate over entries 
                // for(auto& entry: *this)
                // {
                //     // save the entry key 

                // }
            }
        }; 

    } // namespace spectrop
} // namespace fatrop