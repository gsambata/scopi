#include "scopi/problems/DryWithFrictionFixedPoint.hpp"
#include <utility>

namespace scopi
{
    ProblemParams<DryWithFrictionFixedPoint>::ProblemParams()
    : mu(0.)
    , tol_fixed_point(1e-2)
    , max_iter_fixed_point(20)
    {}

    ProblemParams<DryWithFrictionFixedPoint>::ProblemParams(const ProblemParams<DryWithFrictionFixedPoint>& params)
    : mu(params.mu)
    , tol_fixed_point(params.tol_fixed_point)
    , max_iter_fixed_point(params.max_iter_fixed_point)
    {}

    DryWithFrictionFixedPoint::DryWithFrictionFixedPoint(std::size_t nparticles, double dt, const ProblemParams<DryWithFrictionFixedPoint>& params)
    : DryWithFrictionBase(nparticles, dt, params.mu) 
    , m_params(params)
    {}

    bool DryWithFrictionFixedPoint::should_solve_optimization_problem()
    {
        return (xt::linalg::norm(m_s_old - m_s)/(xt::linalg::norm(m_s)+1.) > m_params.tol_fixed_point && m_nb_iter < m_params.max_iter_fixed_point);
    }

}
