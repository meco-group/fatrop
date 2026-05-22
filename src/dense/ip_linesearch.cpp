#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/pd_solver_orig.hpp"
#include "fatrop/dense/pd_solver_resto.hpp"
#include "fatrop/dense/pd_system_orig.hpp"
#include "fatrop/dense/pd_system_resto.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_linesearch.hxx"

using namespace fatrop;

template class fatrop::IpLinesearch<PdSolverOrig<DenseType>, DenseType>;
template class fatrop::IpLinesearch<PdSolverResto<DenseType>, DenseType>;
