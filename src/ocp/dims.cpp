#include "fatrop/ocp/dims.hpp"
#include "fatrop/common/exception.hpp"
using namespace fatrop;
OcpDims::OcpDims(int K, const std::vector<Index> &nu, const std::vector<Index> &nx,
                 const std::vector<Index> &ng, const std::vector<Index> &ng_ineq)
    : K(K), number_of_controls(nu), number_of_states(nx), number_of_eq_constraints(ng),
      number_of_ineq_constraints(ng_ineq)
{
  check_problem_dimensions();
}

OcpDims::OcpDims(int K, std::vector<Index> &&nu, std::vector<Index> &&nx, std::vector<Index> &&ng,
                 std::vector<Index> &&ng_ineq)
    : K(K), number_of_controls(std::move(nu)), number_of_states(std::move(nx)),
      number_of_eq_constraints(std::move(ng)), number_of_ineq_constraints(std::move(ng_ineq))
{
  check_problem_dimensions();
}
void OcpDims::check_problem_dimensions() const {
    // check if the number of controls, states, and constraints are of the correct size
    fatrop_assert_msg(number_of_controls.size() == K, "The number of controls is not of size K.");
    fatrop_assert_msg(number_of_states.size() == K, "The number of states is not of size K.");
    fatrop_assert_msg(number_of_eq_constraints.size() == K, "The number of equality constraints is not of size K.");
    fatrop_assert_msg(number_of_ineq_constraints.size() == K, "The number of inequality constraints is not of size K.");
    // iterate over every time step and do some checks
    for (Index i = 0; i < K; i++) {
        fatrop_assert_msg(number_of_controls[i] >= 0, "The number of controls must be non-negative.");
        fatrop_assert_msg(number_of_states[i] >= 0, "The number of states must be non-negative.");
        fatrop_assert_msg(number_of_eq_constraints[i] >= 0, "The number of equality constraints must be non-negative.");
        fatrop_assert_msg(number_of_ineq_constraints[i] >= 0, "The number of inequality constraints must be non-negative.");
        fatrop_assert_msg(number_of_eq_constraints[i] <= number_of_states[i] + number_of_controls[i],
                          "The number of eq constraints exceeds the number of variables.");
    }
}