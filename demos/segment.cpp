#include <random>

#include <scopi/objects/types/segment.hpp>
#include <scopi/objects/types/sphere.hpp>
#include <scopi/solver.hpp>

int main(int argc, char** argv)
{
    constexpr std::size_t dim = 2;
    double dt                 = .001;
    const double max_radius   = 0.06;
    std::size_t total_it      = 2000;
    std::size_t n_parts       = 200;//50
    
    scopi::initialize("spheres passing between two segments");

    auto& app = scopi::get_app();
    app.add_option("--nparts", n_parts, "Number of particles")->capture_default_str();
    app.add_option("--nite", total_it, "Number of iterations")->capture_default_str();
    app.add_option("--dt", dt, "Time step")->capture_default_str();
    scopi::scopi_container<dim> particles;
    scopi::ScopiSolver<dim> solver(particles);
    SCOPI_PARSE(argc, argv);

    // Parametric test
    std::size_t ntrain = 2;
    std::default_random_engine generator;
    for (std::size_t itrain = 0; itrain < ntrain; ++itrain)
    {
        
        std::uniform_real_distribution<double> distrib_lexit(0.16, 0.24); // lexit reference=0.2
        auto lexit = distrib_lexit(generator); // varying parameter: exit width

        double length=1;//0.4
        double l_mur=0.05;
        double ymur=1.0;
        double start=0.;
        
        scopi::segment<dim> seg1(scopi::type::position_t<dim>{start, ymur}, scopi::type::position_t<dim>{length, ymur});
        scopi::segment<dim> seg2(scopi::type::position_t<dim>{length+lexit, ymur}, scopi::type::position_t<dim>{length+length+lexit, ymur});
        scopi::segment<dim> seg3(scopi::type::position_t<dim>{length, ymur}, scopi::type::position_t<dim>{length, ymur+l_mur});
        scopi::segment<dim> seg4(scopi::type::position_t<dim>{length+lexit, ymur}, scopi::type::position_t<dim>{length+lexit, ymur+l_mur});
        scopi::segment<dim> seg5(scopi::type::position_t<dim>{start, ymur+l_mur}, scopi::type::position_t<dim>{length, ymur+l_mur});
        scopi::segment<dim> seg6(scopi::type::position_t<dim>{length+lexit, ymur+l_mur}, scopi::type::position_t<dim>{length+length+lexit, ymur+l_mur});
        scopi::segment<dim> seg7(scopi::type::position_t<dim>{start, ymur}, scopi::type::position_t<dim>{start, ymur+l_mur});
        scopi::segment<dim> seg8(scopi::type::position_t<dim>{length+length+lexit, ymur}, scopi::type::position_t<dim>{length+length+lexit,ymur+l_mur});

        particles.push_back(seg1, scopi::property<dim>().deactivate());
        particles.push_back(seg2, scopi::property<dim>().deactivate());
        particles.push_back(seg3, scopi::property<dim>().deactivate());
        particles.push_back(seg4, scopi::property<dim>().deactivate());
        particles.push_back(seg5, scopi::property<dim>().deactivate());
        particles.push_back(seg6, scopi::property<dim>().deactivate());
        particles.push_back(seg7, scopi::property<dim>().deactivate());
        particles.push_back(seg8, scopi::property<dim>().deactivate());

        std::uniform_real_distribution<double> distrib_x(0.5, 2*length+lexit-0.5);
        std::uniform_real_distribution<double> distrib_y(0.2, ymur-0.1);//0.7
        //std::uniform_real_distribution<double> distrib_r(0.01, max_radius);//r is constant for me

        for (std::size_t i = 0; i < n_parts; ++i)
        {
            auto x      = distrib_x(generator);
            auto y      = distrib_y(generator);
            auto radius = 0.02;//distrib_r(generator);//0.025
            auto prop   = scopi::property<dim>()
                            .desired_velocity({
                                {start+length+0.5*lexit - x, ymur - y}
            })
                            .mass(1.)
                            .moment_inertia(0.1);

            particles.push_back(scopi::sphere<dim>(
                                    {
                                        {x, y}
            },
                                    radius),
                                prop);
        }

        auto params                                 = solver.get_params();
        params.contact_method_params.dmax           = 2 * dt;
        params.contact_method_params.kd_tree_radius = 4 * max_radius;
        params.solver_params.write_velocity         = true;
        params.solver_params.write_lagrange_multiplier = true;
        params.solver_params.path                   ="test_dir";
        params.solver_params.filename               = "test"+ std::to_string(itrain); // file in which write the solutions// parametric test: + std::to_string(itrain)
        
        solver.run(dt, total_it);

    }
    return 0;
}
