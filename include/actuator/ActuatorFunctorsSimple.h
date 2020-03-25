// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef ACTUATORFUNCTORSSIMPLE_H_
#define ACTUATORFUNCTORSSIMPLE_H_

#include <actuator/ActuatorNGP.h>
#include <actuator/ActuatorBulkSimple.h>
#include <actuator/ActuatorFunctors.h>
#include <NaluEnv.h>

namespace sierra {
namespace nalu {

namespace actsimple {
// tags
struct ZeroArrays{};
struct ComputeLocations{};
struct AssignVelocities{};
struct ComputeForces{};
}

// typedefs
using ActSimpleZero = ActuatorFunctor<ActuatorBulkSimple, actsimple::ZeroArrays, ActuatorExecutionSpace>;
using ActSimpleUpdatePoints = ActuatorFunctor<ActuatorBulkSimple, actsimple::ComputeLocations, Kokkos::DefaultHostExecutionSpace>;
using ActSimpleAssignVel = ActuatorFunctor<ActuatorBulkSimple, actsimple::AssignVelocities, Kokkos::DefaultHostExecutionSpace>;
using ActSimpleComputeForce = ActuatorFunctor<ActuatorBulkSimple,actsimple::ComputeForces, Kokkos::DefaultHostExecutionSpace>;

// declarations
template<>
ActSimpleZero::ActuatorFunctor(ActuatorBulkSimple& actBulk);

template<>
void ActSimpleZero::operator()(const int& index) const;

template <>
ActSimpleUpdatePoints::ActuatorFunctor(ActuatorBulkSimple& actBulk);

template<>
void ActSimpleUpdatePoints::operator()(const int& index) const;

template<>
ActSimpleAssignVel::ActuatorFunctor(ActuatorBulkSimple& actBulk);

template<>
void ActSimpleAssignVel::operator()(const int& index) const;

template<>
ActSimpleComputeForce::ActuatorFunctor(ActuatorBulkSimple& actBulk);

template<>
void ActSimpleComputeForce::operator()(const int& index) const;

struct ActSimpleSetUpThrustCalc{
  using execution_space = ActuatorFixedExecutionSpace;

  ActSimpleSetUpThrustCalc(ActuatorBulkSimple& actBulk);

  void operator()(int index) const;

  ActuatorBulkSimple& actBulk_;
};

struct ActSimpleComputeThrust{
  using execution_space = ActuatorFixedExecutionSpace;

  ActSimpleComputeThrust(ActuatorBulkSimple& actBulk, stk::mesh::BulkData& stkBulk);

  void operator()(int index) const;

  ActuatorBulkSimple& actBulk_;
  stk::mesh::BulkData& stkBulk_;
};

} /* namespace nalu */
} /* namespace sierra */

#endif /* ACTUATORFUNCTORSSIMPLE_H_ */
