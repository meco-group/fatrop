#include "fatrop/ocp/problem_info.hpp"
#include <algorithm>
#include <numeric>
using namespace fatrop;
namespace fatrop
{
    namespace internal
    {
        // arguments: begin and end iterators of a vector, an offset, and begin iterator of output
        // vector at end of the function the output vector will contain [offset, offset +
        // sum(v[:1]), offset + sum(v[:2]), ...]
        template <typename InputIt, typename OutputIt>
        void compute_offsets(InputIt begin, InputIt end, Index offset, OutputIt out)
        {
            Index K = std::distance(begin, end) + 1;
            out[0] = 0;
            std::partial_sum(begin, end - 1, out + 1);
            // add the offset to each element
            std::transform(out, out + K - 1, out, [offset](Index x) { return x + offset; });
        }
    }
}
ProblemInfo<OcpType>::ProblemInfo(const OcpDims &dims)
    : dims(dims), offsets_primal_u(dims.K), offsets_primal_x(dims.K), offsets_g_eq_dyn(dims.K - 1),
      offsets_g_eq_path(dims.K), offsets_g_eq_slack(dims.K)
{
    using namespace internal;
    // compute the number of primal variables
    number_of_primal_variables =
        std::accumulate(dims.number_of_states.begin(), dims.number_of_states.end(), 0) +
        std::accumulate(dims.number_of_controls.begin(), dims.number_of_controls.end(), 0);
    // compute the number of slack variables
    number_of_slack_variables = std::accumulate(dims.number_of_ineq_constraints.begin(),
                                                dims.number_of_ineq_constraints.end(), 0);
    // compute the number of equality constraints
    number_of_eq_constraints =
        std::accumulate(dims.number_of_states.begin() + 1, dims.number_of_states.end(), 0) +
        std::accumulate(dims.number_of_eq_constraints.begin(), dims.number_of_eq_constraints.end(),
                        0) +
        std::accumulate(dims.number_of_ineq_constraints.begin(),
                        dims.number_of_ineq_constraints.end(), 0);
    // set up the offsets for the primal variables
    // the number of variables per stage is the sum of the number of states and controls
    number_of_stage_variables.resize(dims.K);
    std::transform(dims.number_of_states.begin(), dims.number_of_states.end(),
                   dims.number_of_controls.begin(), number_of_stage_variables.begin(),
                   std::plus<Index>());
    offset_primal = 0;
    compute_offsets(number_of_stage_variables.begin(), number_of_stage_variables.end(),
                    offset_primal, offsets_primal_u.begin());
    // the offsets for the states for a stage is the sum of the offsets for the control + the number
    // of controls
    std::transform(offsets_primal_u.begin(), offsets_primal_u.end(),
                   dims.number_of_controls.begin(), offsets_primal_x.begin(), std::plus<Index>());
    // set up the offsets for the slack variables
    offset_slack = 0;
    offsets_slack = std::vector<Index>(dims.K, 0);
    compute_offsets(dims.number_of_ineq_constraints.begin(), dims.number_of_ineq_constraints.end(),
                    offset_slack, offsets_slack.begin());
    offsets_eq = std::vector<Index>(dims.K, 0);
    compute_offsets(dims.number_of_eq_constraints.begin(), dims.number_of_eq_constraints.end(), 0,
                    offsets_eq.begin());
    // compute the number of (dyn, path, slack)-equality constraints
    number_of_g_eq_dyn =
        std::accumulate(dims.number_of_states.begin() + 1, dims.number_of_states.end(), 0);
    number_of_g_eq_path = std::accumulate(dims.number_of_eq_constraints.begin(),
                                          dims.number_of_eq_constraints.end(), 0);
    number_of_g_eq_slack = std::accumulate(dims.number_of_ineq_constraints.begin(),
                                           dims.number_of_ineq_constraints.end(), 0);
    // set up the offsets the equality constraints are ordered as a concatenation of [g_eq_dyn,
    // g_eq_path, g_eq_slack]
    offset_g_eq_dyn = 0;
    offset_g_eq_path = number_of_g_eq_dyn;
    offset_g_eq_slack = offset_g_eq_path + number_of_g_eq_path;
    compute_offsets(dims.number_of_states.begin() + 1, dims.number_of_states.end(), offset_g_eq_dyn,
                    offsets_g_eq_dyn.begin());
    compute_offsets(dims.number_of_eq_constraints.begin(), dims.number_of_eq_constraints.end(),
                    offset_g_eq_path, offsets_g_eq_path.begin());
    compute_offsets(dims.number_of_ineq_constraints.begin(), dims.number_of_ineq_constraints.end(),
                    offset_g_eq_slack, offsets_g_eq_slack.begin());
}