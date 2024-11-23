#include <random>

#include <scopi/objects/types/segment.hpp>
#include <scopi/objects/types/sphere.hpp>//to create spheres
#include <scopi/solver.hpp>//to create the solver and run it

// To write a json file
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>


int main(int argc, char** argv)
{
    constexpr std::size_t dim = 2;
    double dt                 = .001;
    const double max_radius   = 0.05;
    std::size_t total_it      = 500;//500 total time iterations
    std::size_t n_parts       = 20;//100
    
    scopi::initialize("spheres passing between two segments");//just adds a title to your command line option

    auto& app = scopi::get_app();
    app.add_option("--nparts", n_parts, "Number of particles")->capture_default_str();
    app.add_option("--nite", total_it, "Number of iterations")->capture_default_str();
    app.add_option("--dt", dt, "Time step")->capture_default_str();
        

    // Parametric test:train
    // std::size_t ntrain = 10;
    // std::size_t nsim=ntrain;
    // std::string flag_test="train";
    // unsigned int seed= 1;

    //Parametric test:valid
    std::size_t nvalid = 5;
    std::size_t nsim=nvalid;
    std::string flag_test="valid";
    unsigned int seed = 54321;
    
    std::default_random_engine generator;
    generator.seed(seed);

    // Crée un objet JSON vide
    nlohmann::json json_obj;

    //Save the setting of the simulation
    json_obj["Tmax"].push_back(total_it);
    json_obj["Na"].push_back(n_parts);
    json_obj["deltat"].push_back(dt);

    //Distance function
    double computeDistance(double x1, double y1, double x2, double y2) {
    // Using the distance formula
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    }

    for (std::size_t itrain = 0; itrain < nsim; ++itrain)
    {
        scopi::scopi_container<dim> particles;//container of all the particles: why has the scopi container hpp not been included?
        scopi::ScopiSolver<dim> solver(particles);
        SCOPI_PARSE(argc, argv);//allows to get access to internal options

        //Distribution for lexit
        std::uniform_real_distribution<double> distrib_lexit(0.16, 0.24); // lexit reference=0.2
        auto lexit = distrib_lexit(generator); // varying parameter: exit width

        //Distribution for spontaneous velocity
        double sbar=10.0
        std::uniform_real_distribution<double> distrib_s(sbar-0.1*s, 10+0.1*sbar);
        auto s=distrib_s(generator); //For now, the spontaneous velocity is the same for all the particles. It only varies for different geometric configurations

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

        particles.push_back(seg1, scopi::property<dim>().deactivate());//these objects are only obstacles
        particles.push_back(seg2, scopi::property<dim>().deactivate());
        particles.push_back(seg3, scopi::property<dim>().deactivate());
        particles.push_back(seg4, scopi::property<dim>().deactivate());
        particles.push_back(seg5, scopi::property<dim>().deactivate());
        particles.push_back(seg6, scopi::property<dim>().deactivate());
        particles.push_back(seg7, scopi::property<dim>().deactivate());
        particles.push_back(seg8, scopi::property<dim>().deactivate());

        std::uniform_real_distribution<double> distrib_x(0.5, 2*length+lexit-0.5);
        std::uniform_real_distribution<double> distrib_y(0.5, ymur-0.1);//0.2 min
        //std::uniform_real_distribution<double> distrib_r(0.01, max_radius);//r is constant for me

        
        for (std::size_t i = 0; i < n_parts; ++i)
        {
            auto x      = distrib_x(generator);
            auto y      = distrib_y(generator);
            auto radius = 0.05;//distrib_r(generator);//0.05
            double dist_factor=computeDistance(start+length+0.5*lexit,2, x,y)
            auto prop   = scopi::property<dim>()
                            .desired_velocity({
                                {(start+length+0.5*lexit - x)*s/dist_factor, (2 - y)*s/dist_factor}// ymur+l_mur
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
        if (flag_test=="train"){
            params.solver_params.path                   ="test_dir";
        }
            else if (flag_test=="valid"){
                params.solver_params.path                   ="test_dir_valid";   
        }
        params.solver_params.filename               = "test"+ std::to_string(itrain); // file in which write the solutions// parametric test: + std::to_string(itrain)
        
        solver.run(dt, total_it);

        
        // Ajouter des données à l'objet JSON
        json_obj["lexit"].push_back(lexit);//add an element to the end of a container such that vector or list
        //add another field for the spontaneous velocity
        // Sauvegarder l'objet JSON dans un fichier
        if (flag_test=="train"){
            std::ofstream file("lexit_train.json");//ofstream to open the file for writing
            if (file.is_open()) {
                file << json_obj.dump(4); // The updated JSON object is serialized and written back to the file using json_obj.dump(4) (with an indentation of 4 spaces for readability).
                file.close();
                //std::cout << "Fichier JSON créé avec succès !" << std::endl;
            } else {
                std::cerr << "Erreur lors de l'ouverture du fichier !" << std::endl;
        }
        }
            else if (flag_test=="valid"){
                std::ofstream file("lexit_valid.json");//ofstream to open the file for writing
                if (file.is_open()) {
                    file << json_obj.dump(4); // The updated JSON object is serialized and written back to the file using json_obj.dump(4) (with an indentation of 4 spaces for readability).
                    file.close();
                    //std::cout << "Fichier JSON créé avec succès !" << std::endl;
                } else {
                    std::cerr << "Erreur lors de l'ouverture du fichier !" << std::endl;
                }
            }
    

    }
    

    return 0;
}
