// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef ACTUATORFUNCTORSFAST_H_
#define ACTUATORFUNCTORSFAST_H_

#include <actuator/ActuatorNGP.h>
#include <actuator/ActuatorBulkFAST.h>
#include <actuator/ActuatorFunctors.h>

namespace sierra {
namespace nalu {

namespace actfast {
// tags
struct ZeroArrays{};
struct ComputeLocations{};
struct AssignVelocities{};
struct ComputeForces{};
struct SpreadForces{};
struct ComputeThrust{};
}

// typedefs
using ActuatorNgpFAST = Actuator<ActuatorMetaFAST, ActuatorBulkFAST>;
using ActFastZero = ActuatorFunctor<ActuatorBulkFAST, actfast::ZeroArrays, ActuatorExecutionSpace>;
using ActFastUpdatePoints = ActuatorFunctor<ActuatorBulkFAST, actfast::ComputeLocations, Kokkos::DefaultHostExecutionSpace>;
using ActFastAssignVel = ActuatorFunctor<ActuatorBulkFAST, actfast::AssignVelocities, Kokkos::DefaultHostExecutionSpace>;
using ActFastComputeForce = ActuatorFunctor<ActuatorBulkFAST,actfast::ComputeForces, Kokkos::DefaultHostExecutionSpace>;
using ActFastSpreadForce = ActuatorFunctor<ActuatorBulkFAST, actfast::SpreadForces, ActuatorExecutionSpace>;

// declarations
template<>
ActFastZero::ActuatorFunctor(ActuatorBulkFAST& actBulk);

template <>
ActFastUpdatePoints::ActuatorFunctor(ActuatorBulkFAST& actBulk);

template<>
ActFastAssignVel::ActuatorFunctor(ActuatorBulkFAST& actBulk);

template<>
ActFastComputeForce::ActuatorFunctor(ActuatorBulkFAST& actBulk);

template<>
ActFastSpreadForce::ActuatorFunctor(ActuatorBulkFAST& actBulk);

template <>
void
ActuatorNgpFAST::execute()
{
  auto velReduce   = actBulk_.velocity_.template      view<ActuatorFixedMemSpace>();
  auto forceReduce = actBulk_.actuatorForce_.template view<ActuatorFixedMemSpace>();

  Kokkos::parallel_for("zeroQuantitiesActuatorNgpFAST", numActPoints_, ActFastZero(actBulk_));

  // set range policy to only operating over points owned by local fast turbine
  auto fast_range_policy = actBulk_.local_range_policy(actMeta_);
  Kokkos::parallel_for("updatePointLocationsActuatorNgpFAST", fast_range_policy, ActFastUpdatePoints(actBulk_));

  actBulk_.stk_search_act_pnts(actMeta_);

  Kokkos::parallel_for("interpolateVelocitiesActuatorNgpFAST", numActPoints_, InterpolateActVel(actBulk_));

  actBulk_.reduce_view_on_host(velReduce);

  Kokkos::parallel_for("assignFastVelActuatorNgpFAST", fast_range_policy, ActFastAssignVel(actBulk_));

  actBulk_.step_fast();

  Kokkos::parallel_for("computeForcesActuatorNgpFAST", fast_range_policy, ActFastComputeForce(actBulk_));

  actBulk_.reduce_view_on_host(forceReduce);
  // TODO(psakiev) spread forces to nodes
  Kokkos::parallel_for("spreadForcesActuatorNgpFAST", numActPoints_, ActFastSpreadForce(actBulk_));
  // TODO(psakiev) compute thrust
}

} /* namespace nalu */
} /* namespace sierra */

#endif /* ACTUATORFUNCTORSFAST_H_ */
