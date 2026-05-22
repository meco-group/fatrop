#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_resto_phase_min_cl1.hxx"

template class fatrop::IpRestoPhaseMinCl1<fatrop::DenseType>;
