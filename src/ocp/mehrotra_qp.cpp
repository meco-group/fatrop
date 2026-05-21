//
// Copyright (C) 2026 Lander Vanroye, KU Leuven
//
// Explicit instantiation of the Mehrotra predictor-corrector QP solver for
// the OCP problem type, so downstream users only need to include the public
// declaration header (fatrop/qp/fatrop_qp.hpp).
//
#include "fatrop/qp/fatrop_qp.hxx"

#include "fatrop/ocp/aug_system_solver.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/type.hpp"

namespace fatrop
{
    template class MehrotraQpAlgorithm<OcpType>;
    template class MehrotraQpBuilder<OcpType>;
} // namespace fatrop
