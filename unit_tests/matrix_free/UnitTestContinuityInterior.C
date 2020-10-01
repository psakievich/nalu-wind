// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include "Kokkos_Core.hpp"
#include <Kokkos_CopyViews.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <algorithm>
#include <stk_simd/Simd.hpp>
#include <type_traits>

#include "matrix_free/ContinuityInterior.h"
#include "matrix_free/LobattoQuadratureRule.h"
#include "matrix_free/LinearAreas.h"
#include "matrix_free/LinearDiffusionMetric.h"
#include "matrix_free/LinearAdvectionMetric.h"
#include "matrix_free/KokkosViewTypes.h"
#include "matrix_free/MakeRCP.h"
#include "gtest/gtest.h"
#include "mpi.h"

namespace sierra {
namespace nalu {
namespace matrix_free {
namespace test_continuity {

static constexpr int order = 1;
static constexpr int nodes_per_elem = (order + 1) * (order + 1) * (order + 1);
static constexpr int num_elems = 1;

Teuchos::RCP<const Tpetra::Map<>>
make_map()
{
  return make_rcp<Tpetra::Map<>>(
    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
    num_elems * nodes_per_elem, 1,
    make_rcp<Teuchos::MpiComm<int>>(MPI_COMM_WORLD));
}

template <int p>
void
set_aux_fields(elem_offset_view<p> offsets, vector_view<p> coordinates)
{
  constexpr auto nodes = GLL<p>::nodes;
  Kokkos::parallel_for(
    num_elems, KOKKOS_LAMBDA(int index) {
      for (int k = 0; k < p + 1; ++k) {
        const auto cz = nodes[k];
        for (int j = 0; j < p + 1; ++j) {
          const auto cy = nodes[j];
          for (int i = 0; i < p + 1; ++i) {
            const auto cx = nodes[i];
            coordinates(index, k, j, i, 0) = cx;
            coordinates(index, k, j, i, 1) = cy;
            coordinates(index, k, j, i, 2) = cz;
            offsets(index, 0, k, j, i) = 1 + index * simd_len * nodes_per_elem +
                                         k * (order + 1) * (order + 1) +
                                         j * (order + 1) + i;
          }
        }
      }
    });
}

} // namespace test_continuity
class ContinuityResidualFixture : public ::testing::Test
{
public:
  ContinuityResidualFixture()
  {
    scalar_view<order> pressure{"pressure", num_elems};
    Kokkos::deep_copy(pressure, 0);

    vector_view<order> coordinates{"coordinates", num_elems};
    test_continuity::set_aux_fields<order>(offsets, coordinates);

    scalar_view<order> density{"density", num_elems};
    Kokkos::deep_copy(density, 1.0);

    vector_view<order> velocity{"velocity", num_elems};
    Kokkos::deep_copy(velocity, 1.0);

    vector_view<order> gp{"gp", num_elems};
    Kokkos::deep_copy(gp, 0.0);

    auto areas = geom::linear_areas<order>(coordinates);
    metric = geom::diffusion_metric<order>(coordinates);

    geom::linear_advection_metric<order>(
      1., areas, metric, density, velocity, gp, pressure, mdot);
  }
  static constexpr int order = 1;
  static constexpr int nodes_per_elem = (order + 1) * (order + 1) * (order + 1);
  static constexpr int num_elems = 1;

  elem_offset_view<order> offsets{"offsets", num_elems};
  scs_scalar_view<order> mdot{"mdot", num_elems};
  scs_vector_view<order> metric{"metric", num_elems};
  Tpetra::MultiVector<> delta{test_continuity::make_map(), 1};
  Tpetra::MultiVector<> rhs{test_continuity::make_map(), 1};
};

TEST_F(ContinuityResidualFixture, residual_executes)
{
  decltype(rhs.getLocalViewDevice()) shared_rhs("empty_rhs", 1, 1);

  rhs.putScalar(0.);
  continuity_residual<order>(1., offsets, mdot, rhs.getLocalViewDevice());
}

TEST_F(ContinuityResidualFixture, linearized_residual_executes)
{
  delta.putScalar(0.);
  rhs.putScalar(0.);

  rhs.putScalar(0.);
  continuity_linearized_residual<order>(
    offsets, metric, delta.getLocalViewDevice(), rhs.getLocalViewDevice());
}
} // namespace matrix_free
} // namespace nalu
} // namespace sierra