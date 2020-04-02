// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef ACTUATORBULKSIMPLE_H_
#define ACTUATORBULKSIMPLE_H_

#include <actuator/ActuatorBulk.h>
#include "OpenFAST.H"

namespace sierra {
namespace nalu {

struct ActuatorMetaSimple : public ActuatorMeta
{
  ActuatorMetaSimple(const ActuatorMeta& actMeta);

  // HOST ONLY
  fast::fastInputs fastInputs_;
  std::vector<std::string> turbineNames_;
  std::vector<std::string> turbineOutputFileNames_;
  bool filterLiftLineCorrection_;
  bool isotropicGaussian_;
  bool is_disk();
  int get_fast_index(
    fast::ActuatorNodeType type,
    int turbId,
    int index = 0,
    int bladeNum = 0) const;

  // TODO(SAKIEVICH) not certain all these need to be dual views
  int maxNumPntsPerBlade_;
  ActVectorDblDv epsilon_;
  ActVectorDblDv epsilonChord_;
  ActVectorDblDv epsilonTower_;
  ActFixScalarBool useUniformAziSampling_;
  ActFixScalarInt nPointsSwept_;
  ActFixScalarInt nBlades_;
};

struct ActuatorBulkSimple : public ActuatorBulk
{
  ActuatorBulkSimple(const ActuatorMetaSimple& actMeta, double naluTimeStep);

  Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace> local_range_policy();

  void interpolate_velocities_to_fast();
  void step_fast();
  bool fast_is_time_zero();
  void output_torque_info();
  void init_openfast(const ActuatorMetaSimple& actMeta, double naluTimeStep);
  void init_epsilon(const ActuatorMetaSimple& actMeta);
  virtual void zero_open_fast_views();

  virtual ~ActuatorBulkSimple();

  ActFixVectorDbl turbineThrust_;
  ActFixVectorDbl turbineTorque_;
  ActFixVectorDbl hubLocations_;
  ActFixVectorDbl hubOrientation_;

  ActVectorDblDv epsilonOpt_;
  ActTensorDblDv orientationTensor_;

  // TODO(SAKIEVICH) this kill lambdas that are pass by value (KOKKOS_LAMBDA)
  // may need to rethink functor/bulk design.  Perhaps have an internal object
  // in bulk for gpu data and pass that into the actuatorFunctors.
  fast::OpenFAST openFast_;
  const int localTurbineId_;
  const int tStepRatio_;
  ActDualViewHelper<ActuatorMemSpace> dvHelper_;
};

// helper functions to
// squash calls to std::cout from TPL's aka OpenFAST
inline
void squash_fast_output(std::function<void()>func)
{
  std::stringstream buffer;
  std::streambuf* sHoldCout = std::cout.rdbuf();
  std::cout.rdbuf(buffer.rdbuf());
  func();
  std::cout.rdbuf(sHoldCout);
}

} // namespace nalu
} // namespace sierra

#endif /* ACTUATORBULKSIMPLE_H_ */
