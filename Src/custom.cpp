#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "ComplexTraits.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Hamiltonian.hpp"

template <typename T, unsigned D>
void Hamiltonian<T, D>::assign_local_potential(
  const unsigned index_,
  const unsigned io_
)
{
  Eigen::Matrix<double, D, 1> pos;
  Coordinates<std::ptrdiff_t, D + 1> global(r.Lt);
  Coordinates<std::ptrdiff_t, D + 1> local(r.Ld);
  Eigen::Map<Eigen::Matrix<std::ptrdiff_t, D, 1>> v_global(global.coord);
  local.set_coord(index_);
  r.convertCoordinates(global, local);
  pos = r.rLat * v_global.template cast<double>();

  if constexpr (D == 2) {
    const Eigen::Matrix<double, D, 1> origin(0.5 * r.Lt[0], 0.5 * r.Lt[1]);
    const Eigen::Matrix<double, D, 1> rel_pos = pos - r.rLat * origin;
    const value_type radius = rel_pos.norm();
    const value_type boundary = r.Lt[0] / 8;
    if (radius < boundary)
      custom_pot = 0.0;
    else {
      const value_type sigma = r.Lt[0] / 8;
      const value_type Vmax = 2.0;
      const value_type dr = radius - boundary;
      custom_pot = Vmax * (1.0 - std::exp(-0.5 * dr * dr / (sigma * sigma)));
    }
  }
  // General Purpose
  // if constexpr(D == Dimension of the system) {
  //   custom_pot = (...)
  // }
}

#define instantiate(type, dim) template class Hamiltonian<type, dim>;
#include "instantiate.hpp"
