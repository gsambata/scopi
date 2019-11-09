#ifndef DEF_Projection
#define DEF_Projection

#include <iostream>
#include <vector>
#include <typeinfo>

#include "mkl_service.h"
#include "mkl_spblas.h"

#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xdynamic_view.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor-python/pyarray.hpp>


/// @class Projection
/// @bried This class
class Projection
{
public:

  /// @brief Constructor
  /// Instantiate
  Projection(
    const std::size_t maxiter,
    const double rho,
    const double dmin,
    const double tol,
    const double dt
  );

  /// @brief Destructor
  ~Projection();


  // xt::pyarray<double> run(
  void run(
      xt::pyarray<double> &xyzr,
    xt::pyarray<double> &contacts,
    xt::pyarray<double> &V,
    xt::pyarray<double> &D,
    xt::pyarray<double> &invM
  );

  /// @brief Print the contacts
  void print();

private:

  std::size_t _nc; // number of contacts

  const std::size_t _maxiter;
  const double _rho;
  const double _dmin;
  const double _tol;
  const double _dt;

  std::vector<double> _B_coef;
  std::vector<int> _B_col, _B_index;
  matrix_descr _descB;
  sparse_matrix_t _mklB;
  sparse_status_t _stat_B_csr;

  matrix_descr _descBt;
  sparse_matrix_t _mklBt;
  sparse_status_t _stat_Bt_csr;

};

#endif
