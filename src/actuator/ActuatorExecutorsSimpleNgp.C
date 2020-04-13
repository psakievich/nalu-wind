// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//
#include <actuator/ActuatorExecutorsSimpleNgp.h>

namespace sierra {
namespace nalu {

ActuatorLineSimpleNGP::ActuatorLineSimpleNGP(
  const ActuatorMetaSimple& actMeta,
  ActuatorBulkSimple& actBulk,
  stk::mesh::BulkData& stkBulk)
  : actMeta_(actMeta),
    actBulk_(actBulk),
    stkBulk_(stkBulk),
    numActPoints_(actMeta_.numPointsTotal_)
{
}

void
ActuatorLineSimpleNGP::operator()()
{
  auto forceReduce = actBulk_.actuatorForce_.view_host();

  actBulk_.zero_source_terms(stkBulk_);

  //LCCDELETE
  // if (actBulk_.fast_is_time_zero()) {
  //   update();
  // }
  // actBulk_.interpolate_velocities_to_fast(); //LCCDELETE

  update();

  //throw std::runtime_error("ActuatorLineSimpleNGP::update stop");  // LCCSTOP

  //actBulk_.step_fast(); // DELETE

  // set range policy to only operating over points owned by local fast turbine
  auto fastRangePolicy = actBulk_.local_range_policy();

  //throw std::runtime_error("ActuatorLineSimpleNGP::compute start");  // LCCSTOP

  Kokkos::parallel_for(
    "computeForcesActuatorNgpSimple", fastRangePolicy,
    ActSimpleComputeForce(actBulk_, actMeta_));
  //throw std::runtime_error("ActuatorLineSimpleNGP::compute forces");  // LCCSTOP

  actuator_utils::reduce_view_on_host(forceReduce);

  const int localSizeCoarseSearch =
    actBulk_.coarseSearchElemIds_.view_host().extent_int(0);

  if (actMeta_.isotropicGaussian_) {
    Kokkos::parallel_for(
      "spreadForcesActuatorNgpSimple", localSizeCoarseSearch,
      SpreadActuatorForce(actBulk_, stkBulk_));
  } else {
    const int rank = NaluEnv::self().parallel_rank();
    Kokkos::deep_copy(actBulk_.orientationTensor_.view_host(),0.0);
    Kokkos::parallel_for(
      "gatherBladeOrientations", fastRangePolicy,
      ActSimpleStashOrientationVectors(actBulk_));

    actuator_utils::reduce_view_on_host(
      actBulk_.orientationTensor_.view_host());

    Kokkos::parallel_for(
      "spreadForceUsingProjDistance", localSizeCoarseSearch,
      ActSimpleSpreadForceWhProjection(actBulk_, stkBulk_));
  }

  actBulk_.parallel_sum_source_term(stkBulk_);
  //throw std::runtime_error("ActuatorLineSimpleNGP:: parallel_sum_source_term");  // LCCSTOP

  // IGNORE THE THRUST CALC FOR NOW
  bool run_thrustcalc=false;
  if (run_thrustcalc) {
    // Always run compute the thrust and torque
    Kokkos::parallel_for(
			 "setUpTorqueCalc", actMeta_.numberOfActuators_,
			 ActSimpleSetUpThrustCalc(actBulk_));
    //NaluEnv::self().naluOutput() << "setUpTorqueCalc done"<< std::endl; //LCCOUT

    actuator_utils::reduce_view_on_host(actBulk_.hubLocations_);
    actuator_utils::reduce_view_on_host(actBulk_.hubOrientation_);
  
    Kokkos::parallel_for(
			 "computeTorque", localSizeCoarseSearch,
			 ActSimpleComputeThrust(actBulk_, stkBulk_));
    NaluEnv::self().naluOutput() << "computeTorque done"<< std::endl; //LCCOUT
    actuator_utils::reduce_view_on_host(actBulk_.turbineThrust_);
    actuator_utils::reduce_view_on_host(actBulk_.turbineTorque_);
    actBulk_.output_torque_info();
  }
  //throw std::runtime_error("ActuatorLineSimpleNGP:: computeTorque");  // LCCSTOP
}

void
ActuatorLineSimpleNGP::update()
{
  NaluEnv::self().naluOutputP0()  // LCCOUT
    << "Blade: " <<actBulk_.localTurbineId_ << std::endl;

  auto velReduce = actBulk_.velocity_.view_host();
  auto pointReduce = actBulk_.pointCentroid_.view_host();
  actBulk_.zero_open_fast_views();

  // set range policy to only operating over points owned by local fast turbine
  auto fastRangePolicy = actBulk_.local_range_policy();

  // Get p1 and p2 for blade geometry
  std::vector<double> p1(3);
  std::vector<double> p2(3);
  for (int j=0; j<3; j++) { 
    p1[j] = actMeta_.p1_.h_view(actBulk_.localTurbineId_,j);
    p2[j] = actMeta_.p2_.h_view(actBulk_.localTurbineId_,j);
  }
  int Npts=actMeta_.num_force_pts_blade_.h_view(actBulk_.localTurbineId_);

  // if (actBulk_.debug_output_)
  //   NaluEnv::self().naluOutput()  // LCCOUT
  //     << "Blade: " <<actBulk_.localTurbineId_ 
  //     << " p1: "<<p1[0]<<" "<<p1[1]<<" "<<p1[2]
  //     << " p2: "<<p2[0]<<" "<<p2[1]<<" "<<p2[2]<< std::endl;

  Kokkos::parallel_for(
    "updatePointLocationsActuatorNgpSimple", fastRangePolicy,
    ActSimpleUpdatePoints(actBulk_, p1, p2, Npts));
  actuator_utils::reduce_view_on_host(pointReduce);
  actBulk_.stk_search_act_pnts(actMeta_, stkBulk_);

  Kokkos::parallel_for(
    "interpolateVelocitiesActuatorNgpSimple", numActPoints_,
    InterpActuatorVel(actBulk_, stkBulk_));
  actuator_utils::reduce_view_on_host(velReduce);

  // if (actBulk_.debug_output_)
  //   NaluEnv::self().naluOutput() << " Starting density interpolation " <<std::endl;  // LCCOUT

  Kokkos::parallel_for(
    "interpolateDensityActuatorNgpSimple", numActPoints_,
    InterpActuatorDensity(actBulk_, stkBulk_));
  auto rhoReduce = actBulk_.density_.view_host();
  actuator_utils::reduce_view_on_host(rhoReduce);

  // This is for output purposes
  Kokkos::parallel_for(
    "assignSimpleVelActuatorNgpSimple", fastRangePolicy,
    ActSimpleAssignVel(actBulk_));

}

  /*
ActuatorDiskFastNGP::ActuatorDiskFastNGP(
  const ActuatorMetaFAST& actMeta,
  ActuatorBulkDiskFAST& actBulk,
  stk::mesh::BulkData& stkBulk)
  : actMeta_(actMeta),
    actBulk_(actBulk),
    stkBulk_(stkBulk),
    numActPoints_(actMeta_.numPointsTotal_)
{
}
  */

  /*
void
ActuatorDiskFastNGP::operator()()
{
  auto velReduce = actBulk_.velocity_.view_host();
  auto pointReduce = actBulk_.pointCentroid_.view_host();

  if (!actBulk_.searchExecuted_) {
    actBulk_.stk_search_act_pnts(actMeta_, stkBulk_);
  }

  actBulk_.zero_source_terms(stkBulk_);
  actBulk_.zero_open_fast_views();

  Kokkos::parallel_for(
    "interpolateVelocitiesActuatorNgpFAST", numActPoints_,
    InterpActuatorVel(actBulk_, stkBulk_));

  actuator_utils::reduce_view_on_host(velReduce);

  auto fastRangePolicy = actBulk_.local_range_policy();

  Kokkos::parallel_for(
    "assignFastVelActuatorNgpFAST", fastRangePolicy,
    ActFastAssignVel(actBulk_));

  auto forceReduce = actBulk_.actuatorForce_.view_host();

  actBulk_.interpolate_velocities_to_fast();

  actBulk_.step_fast();

  Kokkos::parallel_for(
    "computeForcesActuatorNgpFAST", fastRangePolicy,
    ActFastComputeForce(actBulk_));

  actuator_utils::reduce_view_on_host(forceReduce);

  actBulk_.spread_forces_over_disk(actMeta_);

  const int localSizeCoarseSearch =
    actBulk_.coarseSearchElemIds_.view_host().extent_int(0);

  Kokkos::parallel_for(
    "spreadForcesActuatorNgpFAST", localSizeCoarseSearch,
    SpreadActuatorForce(actBulk_, stkBulk_));

  actBulk_.parallel_sum_source_term(stkBulk_);

  if (actBulk_.openFast_.isDebug()) {
    Kokkos::parallel_for(
      "setUpTorqueCalc", actMeta_.numberOfActuators_,
      ActFastSetUpThrustCalc(actBulk_));

    actuator_utils::reduce_view_on_host(actBulk_.hubLocations_);
    actuator_utils::reduce_view_on_host(actBulk_.hubOrientation_);

    Kokkos::parallel_for(
      "computeTorque", localSizeCoarseSearch,
      ActFastComputeThrust(actBulk_, stkBulk_));
    actuator_utils::reduce_view_on_host(actBulk_.turbineThrust_);
    actuator_utils::reduce_view_on_host(actBulk_.turbineTorque_);
    actBulk_.output_torque_info();
  }
}
  */
} // namespace nalu
} // namespace sierra
