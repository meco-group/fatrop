#include "fatrop/ip_algorithm/ip_algorithm.hpp"
#include "fatrop/ip_algorithm/ip_algorithm.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"
#include "fatrop/ip_algorithm/ip_initializer.hpp"
#include "fatrop/ip_algorithm/ip_linesearch.hpp"
#include "fatrop/ip_algorithm/ip_mu_update.hpp"
#include "fatrop/ip_algorithm/ip_search_dir.hpp"
#include "fatrop/ocp/type.hpp"

using namespace fatrop;
IpAlgorithm::IpAlgorithm(const IpSearchDirSp &search_dir, const IpLineSearchSp &linesearch,
                         const IpInitializerSp &initializer, const IpMuUpdateSp &mu_update,
                         const IpEqMultInitializerSp &eq_mult_initializer)
    : search_dir_(search_dir), linesearch_(linesearch), initializer_(initializer),
      mu_update_(mu_update), eq_mult_initializer_(eq_mult_initializer)
{
}

void IpAlgorithm::reset()
{
    // todo who resets the ipdata?
    search_dir_->reset();
    linesearch_->reset();
    initializer_->reset();
    mu_update_->reset();
    eq_mult_initializer_->reset();
}

SolverReturnFlag IpAlgorithm::optimize(const bool is_resto) { return SolverReturnFlag::Success; }
