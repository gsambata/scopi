#include <random>

#include <scopi/objects/types/segment.hpp>
#include <scopi/objects/types/sphere.hpp>
#include <scopi/solver.hpp>

int main(int argc, char** argv)
{
    constexpr std::size_t dim = 2;
    double dt                 = .001;
    double max_radius         = 0.06;
    std::size_t total_it      = 1000; // 1000, 500
    std::size_t n_parts       = 50;  // 100,20

    scopi::initialize("spheres passing between two segments");
    auto& app = scopi::get_app();
    app.add_option("--nparts", n_parts, "Number of particles")->capture_default_str();
    app.add_option("--nite", total_it, "Number of iterations")->capture_default_str();
    app.add_option("--dt", dt, "Time step")->capture_default_str();

    scopi::scopi_container<dim> particles;
    scopi::ScopiSolver<dim> solver(particles);
    SCOPI_PARSE(argc, argv);

    // Parametric test
    //std::size_t ntrain = 2;//train=10, test=2
    //for (std::size_t itrain = 0; itrain < ntrain; ++itrain)
    //{
        // let the obstacles position and the desidered velocity vary according to a normal distribution
        std::default_random_engine generator;

        //The spontaneous velocity generator should be around
        std::uniform_real_distribution<double> distrib_lexit(0.16, 0.24); // lexit reference=0.2

        auto lexit = distrib_lexit(generator); // varying parameter: exit width
        // compute the position of the obstacles according to lexit
        double xseg1      = 0.4;
        double lmur       = 1.4;
        double xleftseg1  = -1.;
        double xrightseg2 = xseg1 + lexit + lmur;
        double ybottom    = 1.;
        double ytop       = 1.25;

       
        scopi::segment<dim> seg1(scopi::type::position_t<dim>{xleftseg1, ybottom}, scopi::type::position_t<dim>{xseg1, ybottom});
        scopi::segment<dim> seg2(scopi::type::position_t<dim>{xseg1 + lexit, ybottom}, scopi::type::position_t<dim>{xrightseg2, ybottom});
        // other segments to form the obstacles
        // obstacle number 1
        scopi::segment<dim> seg3(scopi::type::position_t<dim>{xseg1, ybottom}, scopi::type::position_t<dim>{xseg1, ytop});
        scopi::segment<dim> seg5(scopi::type::position_t<dim>{xseg1, ytop}, scopi::type::position_t<dim>{xleftseg1, ytop});
        scopi::segment<dim> seg7(scopi::type::position_t<dim>{xleftseg1, ytop}, scopi::type::position_t<dim>{xleftseg1, ybottom});
        // obstacle number 2
        scopi::segment<dim> seg4(scopi::type::position_t<dim>{xseg1 + lexit, ybottom}, scopi::type::position_t<dim>{xseg1 + lexit, ytop});
        scopi::segment<dim> seg6(scopi::type::position_t<dim>{xseg1 + lexit, ytop}, scopi::type::position_t<dim>{xrightseg2, ytop});
        scopi::segment<dim> seg8(scopi::type::position_t<dim>{xrightseg2, ytop}, scopi::type::position_t<dim>{xrightseg2, ybottom});

        
        particles.push_back(seg1, scopi::property<dim>().deactivate());
        particles.push_back(seg2, scopi::property<dim>().deactivate());
        // other segments (to form the barrer obstacles)
        particles.push_back(seg3, scopi::property<dim>().deactivate());
        particles.push_back(seg4, scopi::property<dim>().deactivate());
        particles.push_back(seg5, scopi::property<dim>().deactivate());
        particles.push_back(seg6, scopi::property<dim>().deactivate());
        particles.push_back(seg7, scopi::property<dim>().deactivate());
        particles.push_back(seg8, scopi::property<dim>().deactivate());

        // std::uniform_real_distribution<double> distrib_r(0.01, max_radius);//dont'need the distribution of ratios in this test
        std::uniform_real_distribution<double> distrib_x(0.1, 0.9);
        std::uniform_real_distribution<double> distrib_y(0, 0.7);
        for (std::size_t i = 0; i < n_parts; ++i)
        {
            auto x      = distrib_x(generator);
            auto y      = distrib_y(generator);
            auto radius = 0.05; // distrib_r(generator);
            auto prop   = scopi::property<dim>()
                            .desired_velocity({
                                {0.5 - x, 2 - y}
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

       
        auto params = solver.get_params();
        // params.optim_params.alpha           = 0.2 / (dt * dt);
        //params.optim_params.dynamic_descent         = true;
        params.contact_method_params.dmax           = 2 * dt;
        params.contact_method_params.kd_tree_radius = 4 * max_radius;
        //params.solver_params.path                   = "test_dir_valid";//"test_dir";
        //params.solver_params.filename               = "test"+ std::to_string(itrain); // file in which write the solutions// parametric test: + std::to_string(itrain)
        params.solver_params.write_velocity         = true;
        params.solver_params.write_lagrange_multiplier = true;

        // solver.init_options(app);
        // CLI11_PARSE(app, argc, argv);

        // run the solver for different geometric configurations!
        solver.run(dt, total_it);
    //}
    return 0;
}
