#include <cstddef>
#include "doctest/doctest.h"
#include <random>

#include "test_common.hpp"
#include "utils.hpp"

#include <scopi/objects/types/sphere.hpp>
#include <scopi/vap/vap_fpd.hpp>
#include <scopi/container.hpp>
#include <scopi/solver.hpp>
#include <tuple>

namespace scopi {

    TEST_CASE_TEMPLATE("two spheres asymetrical friction", SolverType, SOLVER_WITH_CONTACT_FRICTION(2, contact_kdtree, vap_fixed), SOLVER_WITH_CONTACT_FRICTION(2, contact_brute_force, vap_fixed))
    {
        static constexpr std::size_t dim = 2;
        double dt = .005;
        std::size_t total_it = 1000;
        double mu = 1./2.;

        sphere<dim> s1({{-0.2, -0.05}}, 0.1);
        sphere<dim> s2({{ 0.2,  0.05}}, 0.1);
        auto p = property<dim>().mass(1.).moment_inertia(0.1);

        scopi_container<dim>particles;
        particles.push_back(s1, p.desired_velocity({{0.25, 0}}));
        particles.push_back(s2, p.desired_velocity({{-0.25, 0}}));

        SolverType solver(particles, dt);
        solver.set_coeff_friction(mu);
        solver.solve(total_it);

        CHECK(diffFile("./Results/scopi_objects_0999.json", "../test/references/two_spheres_asymmetrical_friction.json", tolerance));
    }

    TEST_CASE_TEMPLATE("critical 2d spheres friction", SolverType, SOLVER_WITH_CONTACT_FRICTION(2, contact_kdtree, vap_fixed), SOLVER_WITH_CONTACT_FRICTION(2, contact_brute_force, vap_fixed))
    {
        static constexpr std::size_t dim = 2;
        double dt = .01;
        std::size_t total_it = 100;
        double mu = 1./2.;
        scopi_container<dim> particles;

        int n = 2; // 2*n*n particles
        std::minstd_rand0 generator(0);
        std::uniform_real_distribution<double> distrib_r(0.2, 0.4);
        std::uniform_real_distribution<double> distrib_move_x(-0.1, 0.1);
        std::uniform_real_distribution<double> distrib_move_y(-0.1, 0.1);
        std::uniform_real_distribution<double> distrib_velocity(2., 5.);
        auto prop = property<dim>().mass(1.).moment_inertia(0.1);

        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                double r = distrib_r(generator);
                double x = (i + 0.5) + distrib_move_x(generator);
                double y = (j + 0.5) + distrib_move_y(generator);
                double velocity = distrib_velocity(generator);
                sphere<dim> s1({{x, y}}, r);
                particles.push_back(s1, prop.desired_velocity({{velocity, 0.}}));

                r = distrib_r(generator);
                x = (n + i + 0.5) + distrib_move_x(generator);
                y = (j + 0.5) + distrib_move_y(generator);
                velocity = distrib_velocity(generator);
                sphere<dim> s2({{x, y}}, r);
                particles.push_back(s2, prop.desired_velocity({{-velocity, 0.}}));
            }
        }

        SolverType solver(particles, dt);
        solver.set_coeff_friction(mu);
        solver.solve(total_it);

        CHECK(diffFile("./Results/scopi_objects_0099.json", "../test/references/2d_case_spheres_friction.json", tolerance));
    }

    TEST_CASE_TEMPLATE("sphere inclined plan", SolverType, SOLVER_WITH_CONTACT_FRICTION(2, contact_kdtree, vap_fpd), SOLVER_WITH_CONTACT_FRICTION(2, contact_brute_force, vap_fpd))
    {
        std::tuple<double, std::size_t, double, double, double, double> data;
        std::vector<std::tuple<double, std::size_t, double, double, double, double>>
            data_container({std::make_tuple(0.005, 2000, 0.1, PI/6., 0.000508866, 0.000485907),
                            std::make_tuple(0.005, 2000, 0.1, PI/4., 0.000513437, 0.000464869),
                            std::make_tuple(0.005, 2000, 0.1, PI/3., 0.000516023, 0.000428305),
                            std::make_tuple(0.005, 2000, 0.5, PI/6., 0.00049975, 0.000499785),
                            std::make_tuple(0.005, 2000, 0.5, PI/4., 0.000499728, 0.000500104),
                            std::make_tuple(0.005, 2000, 0.5, PI/3., 0.000554791, 0.00044186), 
                            std::make_tuple(0.005, 2000, 1., PI/6., 0.000499751, 0.000499752), 
                            std::make_tuple(0.005, 2000, 1., PI/4., 0.000499751, 0.00049976),
                            std::make_tuple(0.005, 2000, 1., PI/3., 0.000499735, 0.000499765)});

        DOCTEST_VALUE_PARAMETERIZED_DATA(data, data_container);

        static constexpr std::size_t dim = 2;
        double radius = 1.;
        double g = 1.;
        double dt = std::get<0>(data);
        std::size_t total_it = std::get<1>(data);
        double mu = std::get<2>(data);
        double alpha = std::get<3>(data);

        auto prop = property<dim>().mass(1.).moment_inertia(1.*radius*radius/2.);
        plan<dim> p({{-radius*std::cos(PI/2.-alpha), -radius*std::sin(PI/2.-alpha)}}, PI/2.-alpha);
        sphere<dim> s({{0., 0.}}, radius);

        scopi_container<dim> particles;
        particles.push_back(p, property<dim>().deactivate());
        particles.push_back(s, prop.force({{0., -g}}));

        SolverType solver(particles, dt);
        solver.set_coeff_friction(mu);
        solver.solve(total_it);

        auto pos = particles.pos();
        auto omega = particles.omega();
        auto sol = scopi::analytical_solution_sphere_plan(alpha, mu, dt*(total_it+1), radius, g);
        auto pos_analytical = sol.first;
        auto omega_analytical = sol.second;
        double err_pos = xt::linalg::norm(pos(1) - pos_analytical) / xt::linalg::norm(pos_analytical);
        double err_omega = std::abs((omega(1)-omega_analytical)/omega_analytical);

        REQUIRE(err_pos == doctest::Approx(std::get<4>(data)));
        REQUIRE(err_omega == doctest::Approx(std::get<5>(data)));
    }

}
