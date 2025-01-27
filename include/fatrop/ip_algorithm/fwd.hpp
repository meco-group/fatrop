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
    template <typename ProblemType>
    class IpSearchDirImpl;
    template <typename ProblemType> class PdSolverOrig;
    template <typename ProblemType> class PdSystemType;
    template <typename ProblemType> struct IpIterate;
    template <typename ProblemType>
    class IpEqMultInitializer;
    class IpSearchDirBase;
    template <typename ProblemType>
    class IpSearchDirImpl;
    template <typename ProblemType>
    class IpNlpOrig;
    template <typename ProblemType>
    class IpAlgorithm;
    class IpFilterData;
    class IpFilter;

} // namespace fatrop

#endif // __fatrop_ip_algorithm_fwd_hpp__
