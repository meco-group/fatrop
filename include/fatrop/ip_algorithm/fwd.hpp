//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_fwd_hpp__
#define __fatrop_ip_algorithm_fwd_hpp__

namespace fatrop
{
    class IpAlgBuilder;
    template <typename ProblemType> struct IpData;
    template <typename ProblemType> struct IpIterate;
    template <typename ProblemType, typename LinearSystemType, typename LinearSolverType>
    class IpSearchDirImpl;
    template <typename ProblemType> class PdSolverOrig;
    template <typename ProblemType> class PdSystemType;

} // namespace fatrop

#endif // __fatrop_ip_algorithm_fwd_hpp__
