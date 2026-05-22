#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/pd_solver_orig.hpp"
#include "fatrop/dense/pd_solver_resto.hpp"
#include "fatrop/dense/pd_system_orig.hpp"
#include "fatrop/dense/pd_system_resto.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_search_dir.hxx"

using namespace fatrop;
// explicit template instantiation
template class fatrop::IpSearchDirImpl<PdSolverOrig<DenseType>, DenseType>;
template class fatrop::IpSearchDirImpl<PdSolverResto<DenseType>, DenseType>;
