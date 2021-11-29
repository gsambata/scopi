#include <xtensor/xmath.hpp>
#include <scopi/objects/types/sphere.hpp>
#include <scopi/objects/types/superellipsoid.hpp>
#include <scopi/solvers/mosek.hpp>

int main()
{
    constexpr std::size_t dim = 2;
    double PI = xt::numeric_constants<double>::PI;
    scopi::scopi_container<dim> particles;

    scopi::sphere<dim> s1({{-0.2, -0.05}}, 0.1);
    scopi::sphere<dim> s2({{ 0.2,  0.05}}, 0.2);
    scopi::superellipsoid<dim> s3({{-0.2, 0.}}, {scopi::quaternion(PI/4)}, {{.1, .05}}, {{1, 1}});
    particles.push_back(s1, {{0, 0}}, {{0.25, 0}}, 0, 0, {{0, 0}});
    particles.push_back(s2, {{0, 0}}, {{-0.25, 0}}, 0, 0, {{0, 0}});
    particles.push_back(s3, {{0, 0}}, {{-0.25, 0}}, 0, 0, {{0, 0}});

    particles.push_back(0, {{1., 1.}});
    for(std::size_t i = 0; i < particles.size(); ++i)
    {
        particles[i]->print();
    }
    std::cout << std::endl;

    particles.reset_periodic();
    for(std::size_t i = 0; i < particles.size(); ++i)
    {
        particles[i]->print();
    }

    return 0;
}
