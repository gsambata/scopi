#include "scopi/particles.hpp"


Particles::~Particles() {
}

Particles::Particles(const std::string &name) : _name(name) {
  _idmax = 0;
}

void Particles::add_particles_in_ball(
  std::size_t n,
  double rmin,
  double rmax,
  double xc,
  double yc,
  double zc,
  double rball
) {

  const double pi = std::acos(-1);

  xt::xtensor<double, 1> phi = xt::random::rand<double>({n})*2*pi;
  xt::xtensor<double, 1> theta = xt::acos(-1+xt::random::rand<double>({n})*2);
  xt::xtensor<double, 1> rb = rball*xt::random::rand<double>({n});

  xt::xtensor<std::size_t, 1> id = xt::arange<std::size_t>(_idmax,_idmax+n);
  _idmax += n;
  xt::xtensor<double, 1> x = xc+rb * xt::sin( theta) * xt::cos( phi );
  xt::xtensor<double, 1> y = yc+rb * xt::sin( theta) * xt::sin( phi );
  xt::xtensor<double, 1> z = zc+rb * xt::cos( theta );
  xt::xtensor<double, 1> r = rmin+xt::random::rand<double>({n})*(rmax-rmin);

  id = xt::concatenate(xt::xtuple(data.id,id));
  x = xt::concatenate(xt::xtuple(data.x,x));
  y = xt::concatenate(xt::xtuple(data.y,y));
  z = xt::concatenate(xt::xtuple(data.z,z));
  r = xt::concatenate(xt::xtuple(data.r,r));

  data.resize(x.size());

  data.id = id;
  data.x = x;
  data.y = y;
  data.z = z;
  data.r = r;

}

void Particles::add_particles_in_box(
  std::size_t n,
  double rmin,
  double rmax,
  double xmin,
  double xmax,
  double ymin,
  double ymax,
  double zmin,
  double zmax
) {

  xt::xtensor<std::size_t, 1> id = xt::arange<std::size_t>(_idmax,_idmax+n);
  _idmax += n;
  xt::xtensor<double, 1> x = xmin+xt::random::rand<double>({n})*(xmax-xmin);
  xt::xtensor<double, 1> y = ymin+xt::random::rand<double>({n})*(ymax-ymin);
  xt::xtensor<double, 1> z = zmin+xt::random::rand<double>({n})*(zmax-zmin);
  xt::xtensor<double, 1> r = rmin+xt::random::rand<double>({n})*(rmax-rmin);

  id = xt::concatenate(xt::xtuple(data.id,id));
  x = xt::concatenate(xt::xtuple(data.x,x));
  y = xt::concatenate(xt::xtuple(data.y,y));
  z = xt::concatenate(xt::xtuple(data.z,z));
  r = xt::concatenate(xt::xtuple(data.r,r));

  data.resize(x.size());

  data.id = id;
  data.x = x;
  data.y = y;
  data.z = z;
  data.r = r;

}


xt::xtensor<double, 2> Particles::get_data() const{
  return xt::stack(
    xt::xtuple(
      xt::flatten(data.x),
      xt::flatten(data.y),
      xt::flatten(data.z),
      xt::flatten(data.r)
    ),
    1
  );
}

xt::xtensor<double,1> Particles::get_x() const{
  return data.x;
}

xt::xtensor<double,1> Particles::get_y() const{
  return data.y;
}

xt::xtensor<double,1> Particles::get_z() const{
  return data.z;
}

xt::xtensor<double,1> Particles::get_r() const{
  return data.r;
}



void Particles::print() const{

  std::cout<<"\n-- C++ -- Particles : name = "<<_name<<std::endl;
  auto print_data = [](const auto &p) {
    std::cout << "-- C++ -- Particles : id,x,y,z,r = " << p.id << ", "<< p.x << ", " << p.y  << ", " << p.z  << ", " << p.r << "\n";
  };
  std::for_each(data.begin(), data.end(), print_data);
  std::cout<<"-- C++ -- Particles : data memory address = "<< &data.x[0] << std::endl;

}
