//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_search_dir_hpp__
#define __fatrop_ip_algorithm_ip_search_dir_hpp__
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/common/fwd.hpp"
#include <memory>

namespace fatrop
{
    /**
     * @brief Base class for search direction computation in interior point methods.
     */
    class IpSearchDirBase
    {
    public:
        /**
         * @brief Reset the search direction computation to its initial state.
         */
        virtual void reset() = 0;
        virtual void register_options(OptionRegistry& registry) = 0;

        /**
         * @brief Compute the search direction.
         * @return LinsolReturnFlag Indicating the success or failure of the computation.
         */
        virtual LinsolReturnFlag compute_search_dir() = 0;
    protected:
        virtual ~IpSearchDirBase() = default;
    };

    /**
     * @brief Implementation of search direction computation for a specific problem type.
     * 
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType>
    class IpSearchDirImpl : public IpSearchDirBase
    {
        typedef IpIterate<ProblemType> IpIterateType;
        typedef std::shared_ptr<IpIterateType> IpIterateSp;
        typedef IpData<ProblemType> IpDataType;
        typedef std::shared_ptr<IpDataType> IpDataSp;
        typedef LinearSolver<PdSolverOrig<ProblemType>, PdSystemType<ProblemType>> LinearSolverType;
        typedef std::shared_ptr<LinearSolverType> LinearSolverSp;

    public:
        /**
         * @brief Construct a new IpSearchDirImpl object.
         * 
         * @param ipdata Shared pointer to the interior point algorithm data.
         * @param linear_solver Shared pointer to the linear solver.
         */
        IpSearchDirImpl(const IpDataSp &ipdata, const LinearSolverSp &linear_solver);

        void reset() override;
        LinsolReturnFlag compute_search_dir() override;

    private:
        IpDataSp ipdata_;
        LinearSolverSp linear_solver_;
        VecRealAllocated rhs_x_;
        VecRealAllocated rhs_s_;
        VecRealAllocated rhs_g_;
        VecRealAllocated rhs_cl_;
        VecRealAllocated rhs_cu_;
        Scalar delta_w_last_ = 0.;
        Scalar delta_w0_ = 1e-4;
        Scalar delta_wmin_ = 1e-20;
        Scalar kappa_wmin_ = 1. / 3.;
        Scalar kappa_wplus_ = 8.;
        Scalar kappa_wplusem_ = 100.;
        Scalar kappa_c_ = 0.25;
        Scalar delta_c_stripe_ = 1e-6;

    public:
        // Setter methods for options
        void set_delta_w0(const Scalar& value) { delta_w0_ = value; }
        void set_delta_wmin(const Scalar& value) { delta_wmin_ = value; }
        void set_kappa_wmin(const Scalar& value) { kappa_wmin_ = value; }
        void set_kappa_wplus(const Scalar& value) { kappa_wplus_ = value; }
        void set_kappa_wplusem(const Scalar& value) { kappa_wplusem_ = value; }
        void set_kappa_c(const Scalar& value) { kappa_c_ = value; }
        void set_delta_c_stripe(const Scalar& value) { delta_c_stripe_ = value; }

        // Register options
        void register_options(OptionRegistry& registry);
    };
} // namespace fatrop

#endif
