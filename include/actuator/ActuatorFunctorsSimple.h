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

#include <actuator/ActuatorTypes.h>
#include <actuator/ActuatorBulkSimple.h>
#include <actuator/ActuatorFunctors.h>
#include <NaluEnv.h>

namespace sierra {
namespace nalu {

struct ActSimpleUpdatePoints
{
  using execution_space = ActuatorFixedExecutionSpace;

  ActSimpleUpdatePoints(ActuatorBulkSimple& actBulk);
  void operator()(int index) const;

  ActDualViewHelper<ActuatorFixedMemSpace> helper_;
  ActFixVectorDbl points_;
  ActFixScalarInt offsets_;
  const int turbId_;
  fast::OpenFAST& fast_;
};

struct ActSimpleAssignVel
{
  using execution_space = ActuatorFixedExecutionSpace;

  ActSimpleAssignVel(ActuatorBulkSimple& actBulk);
  void operator()(int index) const;

  ActDualViewHelper<ActuatorFixedMemSpace> helper_;
  ActFixVectorDbl velocity_;
  ActFixScalarInt offset_;
  const int turbId_;
  fast::OpenFAST& fast_;
};

struct ActSimpleComputeForce
{
  using execution_space = ActuatorFixedExecutionSpace;

  ActSimpleComputeForce(ActuatorBulkSimple& actBulk);
  void operator()(int index) const;

  ActDualViewHelper<ActuatorFixedMemSpace> helper_;
  ActFixVectorDbl force_;
  ActFixScalarInt offset_;
  const int turbId_;
  fast::OpenFAST& fast_;
};

struct ActSimpleSetUpThrustCalc
{
  using execution_space = ActuatorFixedExecutionSpace;

  ActSimpleSetUpThrustCalc(ActuatorBulkSimple& actBulk);

  void operator()(int index) const;

  ActuatorBulkSimple& actBulk_;
};

struct ActSimpleStashOrientationVectors
{
  using execution_space = ActuatorFixedExecutionSpace;

  ActSimpleStashOrientationVectors(ActuatorBulkSimple& actBulk);

  void operator()(int index) const;

  ActDualViewHelper<ActuatorFixedMemSpace> helper_;
  ActFixTensorDbl orientation_;
  ActFixScalarInt offset_;
  const int turbId_;
  fast::OpenFAST& fast_;
};

struct ActSimpleComputeThrustInnerLoop
{

  ActSimpleComputeThrustInnerLoop(ActuatorBulkSimple& actBulk) : actBulk_(actBulk)
  {
  }
  void operator()(
    const uint64_t pointId,
    const double* nodeCoords,
    double* sourceTerm,
    const double dualNodalVolume,
    const double scvIp) const;
  void preloop() {}

  ActuatorBulkSimple& actBulk_;
};

struct ActSimpleSpreadForceWhProjInnerLoop
{
  ActSimpleSpreadForceWhProjInnerLoop(ActuatorBulkSimple& actBulk)
    : actBulk_(actBulk)
  {
  }

  void operator()(
    const uint64_t pointId,
    const double* nodeCoords,
    double* sourceTerm,
    const double dualNodalVolume,
    const double scvIp) const;
  void preloop();

  ActuatorBulkSimple& actBulk_;
};

using ActSimpleComputeThrust = GenericLoopOverCoarseSearchResults<
  ActuatorBulkSimple,
  ActSimpleComputeThrustInnerLoop>;

using ActSimpleSpreadForceWhProjection = GenericLoopOverCoarseSearchResults<
  ActuatorBulkSimple,
  ActSimpleSpreadForceWhProjInnerLoop>;

} /* namespace nalu */
} /* namespace sierra */

#endif /* ACTUATORFUNCTORSSIMPLE_H_ */
