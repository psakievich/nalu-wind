// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include <actuator/ActuatorFunctorsSimple.h>
#include <actuator/UtilitiesActuator.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <NaluEnv.h>
#include <FieldTypeDef.h>
#include "utils/LinearInterpolation.h"
#include <cmath>

namespace sierra {
namespace nalu {

InterpActuatorDensity::InterpActuatorDensity(
  ActuatorBulk& actBulk, stk::mesh::BulkData& stkBulk)
  : actBulk_(actBulk),
    stkBulk_(stkBulk),
    coordinates_(stkBulk_.mesh_meta_data().get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "coordinates")),
    density_(stkBulk_.mesh_meta_data().get_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "density"))
{
  actBulk_.density_.sync_host();
  actBulk_.density_.modify_host();
}

void
InterpActuatorDensity::operator()(int index) const
{
  auto rho = actBulk_.density_.view_host();
  auto localCoord = actBulk_.localCoords_;

  if (actBulk_.pointIsLocal_(index)) {

    stk::mesh::Entity elem = stkBulk_.get_entity(
      stk::topology::ELEMENT_RANK, actBulk_.elemContainingPoint_(index));

    const int nodesPerElem = stkBulk_.num_nodes(elem);

    // just allocate for largest expected size (hex27)
    double ws_coordinates[81], ws_density[81];

    actuator_utils::gather_field(
      3, &ws_coordinates[0], *coordinates_, stkBulk_.begin_nodes(elem),
      nodesPerElem);

    actuator_utils::gather_field_for_interp(
      1, &ws_density[0], *density_, stkBulk_.begin_nodes(elem), nodesPerElem);

    actuator_utils::interpolate_field(
      1, elem, stkBulk_, &(localCoord(index, 0)), &ws_density[0],
      &(rho(index)));
    rho(index) /= actBulk_.localParallelRedundancy_(index);
    // for (int i = 0; i < 3; i++) {
    //   vel(index, i) /= actBulk_.localParallelRedundancy_(index);
    // }
  }
}

ActSimpleUpdatePoints::ActSimpleUpdatePoints(ActuatorBulkSimple& actBulk, 
					     std::vector<double> p1,
					     std::vector<double> p2,
					     int numpts)
  : points_(helper_.get_local_view(actBulk.pointCentroid_)),
    offsets_(helper_.get_local_view(actBulk.turbIdOffset_)),
    turbId_(actBulk.localTurbineId_),
    fast_(actBulk.openFast_),
    p1_(p1),
    p2_(p2),
    numpoints_(numpts)
{
  helper_.touch_dual_view(actBulk.pointCentroid_);
}

void
ActSimpleUpdatePoints::operator()(int index) const
{

  ThrowAssert(turbId_ >= 0);
  const int pointId = index - offsets_(turbId_);
  auto point = Kokkos::subview(points_, index, Kokkos::ALL);
  std::vector<double> dx(3, 0.0);
  double denom = (double)numpoints_;

  for (int i=0; i<3; i++) {
    dx[i] = (p2_[i] - p1_[i])/denom; 
  }

  for (int i=0; i<3; i++) {
    point(i) = p1_[i] + 0.5*dx[i] + dx[i]*(float)pointId;
  }

  // NaluEnv::self().naluOutputP0()  // LCCOUT
  //   << "pointId: " << pointId
  //   << " point: "<<point(0)<<" "<<point(1)<<" "<<point(2)<<" "<< std::endl;

  //fast_.getForceNodeCoordinates(point.data(), pointId, turbId_);
}

ActSimpleAssignVel::ActSimpleAssignVel(ActuatorBulkSimple& actBulk)
  : velocity_(helper_.get_local_view(actBulk.velocity_)),
    density_(helper_.get_local_view(actBulk.density_)),
    points_(helper_.get_local_view(actBulk.pointCentroid_)),
    offset_(helper_.get_local_view(actBulk.turbIdOffset_)),
    turbId_(actBulk.localTurbineId_),
    debug_output_(actBulk.debug_output_),
    fast_(actBulk.openFast_)
{
}

void
ActSimpleAssignVel::operator()(int index) const
{

  const int pointId = index - offset_(turbId_);
  auto vel = Kokkos::subview(velocity_, index, Kokkos::ALL);
  auto rho = Kokkos::subview(density_, index);

  // Use this to double check the velocities and point positions
  auto point = Kokkos::subview(points_, index, Kokkos::ALL);
  if (debug_output_)
    NaluEnv::self().naluOutput() 
      << "Blade "<< turbId_  // LCCOUT
      << " pointId: " << pointId << std::scientific<< std::setprecision(5)
      << " point: "<<point(0)<<" "<<point(1)<<" "<<point(2)<<" "
      << " vel: "<<vel(0)<<" "<<vel(1)<<" "<<vel(2)<<" "
      << " rho: "<< *rho.data() //*rhoptr
      << std::endl;
  // Do nothing otherwise

  //fast_.setVelocityForceNode(vel.data(), pointId, turbId_);
}

ActSimpleComputeForce::ActSimpleComputeForce(ActuatorBulkSimple& actBulk,
					     const ActuatorMetaSimple& actMeta)
  : force_(helper_.get_local_view(actBulk.actuatorForce_)),
    velocity_(helper_.get_local_view(actBulk.velocity_)),
    density_(helper_.get_local_view(actBulk.density_)),
    offset_(helper_.get_local_view(actBulk.turbIdOffset_)),
    turbId_(actBulk.localTurbineId_),
    debug_output_(actBulk.debug_output_),
    fast_(actBulk.openFast_) // DELETE THIS LATER
{
  // Used to initialize via
    // aoatable_(actMeta.aoa_polartable_[actBulk.localTurbineId_]),
    // cltable_(actMeta.cl_polartable_[actBulk.localTurbineId_]),
    // cdtable_(actMeta.cd_polartable_[actBulk.localTurbineId_]),
    // twist_table_(actMeta.twist_table_[actBulk.localTurbineId_]),
    // elem_area_(actMeta.elem_area_[actBulk.localTurbineId_]),

  helper_.touch_dual_view(actBulk.actuatorForce_);
  if (NaluEnv::self().parallel_rank() == turbId_) {
    // Set up the polar table arrays
    size_t Npolartable=actMeta.aoa_polartable_[actBulk.localTurbineId_].size();
    aoatable_.resize(Npolartable);
    cltable_.resize(Npolartable);
    cdtable_.resize(Npolartable);
    // Copy over the polar tables
    for (size_t i=0; i<Npolartable; i++) {
      aoatable_[i] = actMeta.aoa_polartable_[actBulk.localTurbineId_][i];
      cltable_[i]  = actMeta.cl_polartable_[actBulk.localTurbineId_][i];
      cdtable_[i]  = actMeta.cd_polartable_[actBulk.localTurbineId_][i];
    }
    // Copy over the twist/area tables
    size_t Npts = actMeta.twist_table_[actBulk.localTurbineId_].size();
    twist_table_.resize(Npts);
    elem_area_.resize(Npts);
    for (size_t i=0; i<Npts; i++) {
      twist_table_[i] = actMeta.twist_table_[actBulk.localTurbineId_][i];
      elem_area_[i] = actMeta.elem_area_[actBulk.localTurbineId_][i];
    }
    

    // Set up the directions
    p1zeroalphadir_.x_ = actMeta.p1zeroalphadir_.h_view(turbId_, 0);
    p1zeroalphadir_.y_ = actMeta.p1zeroalphadir_.h_view(turbId_, 1);
    p1zeroalphadir_.z_ = actMeta.p1zeroalphadir_.h_view(turbId_, 2);

    chordnormaldir_.x_ = actMeta.chordnormaldir_.h_view(turbId_, 0);
    chordnormaldir_.y_ = actMeta.chordnormaldir_.h_view(turbId_, 1);
    chordnormaldir_.z_ = actMeta.chordnormaldir_.h_view(turbId_, 2);

    spandir_.x_        = actMeta.spandir_.h_view(turbId_, 0);
    spandir_.y_        = actMeta.spandir_.h_view(turbId_, 1);
    spandir_.z_        = actMeta.spandir_.h_view(turbId_, 2);

  // if (debug_output_)
  //   NaluEnv::self().naluOutput() 
  //     << "Blade "<< turbId_  // LCCOUT
  //     << std::scientific<< std::setprecision(7)
  //     << " p1zeroAOA: "
  //     <<p1zeroalphadir_.x_<<" "<<p1zeroalphadir_.y_<<" "<<p1zeroalphadir_.z_
  //     << " chordnorm: "
  //     <<chordnormaldir_.x_<<" "<<chordnormaldir_.y_<<" "<<chordnormaldir_.z_
  //     << " spandir: "
  //     <<spandir_.x_<<" "<<spandir_.y_<<" "<<spandir_.z_<<" "<<std::endl;
  }
}

void
ActSimpleComputeForce::operator()(int index) const
{

  auto pointForce = Kokkos::subview(force_, index, Kokkos::ALL);
  const int localId = index - offset_(turbId_);

  auto vel     = Kokkos::subview(velocity_, index, Kokkos::ALL);
  auto density = Kokkos::subview(density_, index);

  if (NaluEnv::self().parallel_rank() == turbId_) {

  double twist = twist_table_[localId]; 
  Coordinates ws; // Total wind speed
  ws.x_ = vel(0);
  ws.y_ = vel(1);
  ws.z_ = vel(2);

  // Calculate the angle of attack (AOA)
  double alpha;
  Coordinates ws2D;
  AirfoilTheory2D::calculate_alpha(ws, p1zeroalphadir_, 
				   spandir_, chordnormaldir_, twist, 
				   ws2D, alpha);
  // Calculate Cl and Cd
  double cl;
  double cd;
  AirfoilTheory2D::calculate_cl_cd(alpha, aoatable_, cltable_, cdtable_,
				   cl, cd);

  // Magnitude of wind speed
  double ws2Dnorm = sqrt(ws2D.x_*ws2D.x_ + 
			 ws2D.y_*ws2D.y_ +
			 ws2D.z_*ws2D.z_);
  
  // Calculate lift and drag forces
  double rho  = *density.data();
  rho = 1.0;
  double area = elem_area_[localId];
  double Q    = 0.5*rho*ws2Dnorm*ws2Dnorm;
  double lift = cl*Q*area;
  double drag = cd*Q*area;

  // Set the directions
  Coordinates ws2Ddir;  // Direction of drag force
  if (ws2Dnorm > 0.0) {
    ws2Ddir.x_ = ws2D.x_/ws2Dnorm;
    ws2Ddir.y_ = ws2D.y_/ws2Dnorm;
    ws2Ddir.z_ = ws2D.z_/ws2Dnorm;
  } else {
    ws2Ddir.x_ = 0.0; 
    ws2Ddir.y_ = 0.0; 
    ws2Ddir.z_ = 0.0; 
  }
  Coordinates liftdir;  // Direction of lift force
  if (ws2Dnorm > 0.0) {
    liftdir.x_ = ws2Ddir.y_*spandir_.z_ - ws2Ddir.z_*spandir_.y_; 
    liftdir.y_ = ws2Ddir.z_*spandir_.x_ - ws2Ddir.x_*spandir_.z_; 
    liftdir.z_ = ws2Ddir.x_*spandir_.y_ - ws2Ddir.y_*spandir_.x_; 
  } else {
    liftdir.x_ = 0.0; 
    liftdir.y_ = 0.0; 
    liftdir.z_ = 0.0; 
  }

  // Set the pointForce
  pointForce(0) = -(lift*liftdir.x_ + drag*ws2Ddir.x_);
  pointForce(1) = -(lift*liftdir.y_ + drag*ws2Ddir.y_);
  pointForce(2) = -(lift*liftdir.z_ + drag*ws2Ddir.z_);

  if (debug_output_)
    NaluEnv::self().naluOutput() 
      << "Blade "<< turbId_  // LCCOUT // << std::scientific
      << " pointId: " << localId << std::setprecision(5)
      << " alpha: "<<alpha
      << " ws2D: "<<ws2D.x_<<" "<<ws2D.y_<<" "<<ws2D.z_<<" "
      << " Cl, Cd: "<<cl<<" "<<cd
      << " lift, drag = "<<lift<<" "<<drag
      //<< " pForce: "<<pointForce(0)<<" "<<pointForce(1)<<" "<<pointForce(2)
      << std::endl;
  }
  //fast_.getForce(pointForce.data(), localId, turbId_);
}

void 
AirfoilTheory2D::calculate_alpha(
    Coordinates ws, 
    Coordinates zeroalphadir, 
    Coordinates spandir,
    Coordinates chordnormaldir, 
    double twist, 
    Coordinates &ws2D, 
    double &alpha) 
{
  // Project WS onto 2D plane defined by zeroalpahdir and chordnormaldir
  double WSspan = ws.x_*spandir.x_ + ws.y_*spandir.y_ + ws.z_*spandir.z_;
  ws2D.x_ = ws.x_ - WSspan*spandir.x_;
  ws2D.y_ = ws.y_ - WSspan*spandir.y_;
  ws2D.z_ = ws.z_ - WSspan*spandir.z_;

  // Project WS2D onto zeroalphadir and chordnormaldir
  double WStan = 
    ws2D.x_*zeroalphadir.x_ + 
    ws2D.y_*zeroalphadir.y_ +  
    ws2D.z_*zeroalphadir.z_ ;
  
  double WSnormal = 
    ws2D.x_*chordnormaldir.x_ + 
    ws2D.y_*chordnormaldir.y_ + 
    ws2D.z_*chordnormaldir.z_ ;
  

  double alphaNoTwist = atan2(WSnormal, WStan)*180.0/M_PI;

  // NaluEnv::self().naluOutput() 
  //   << "alphaNoTwist = "<<alphaNoTwist
  //   << " WSnormal = "<<WSnormal<<" WStan = "<<WStan 
  //   << " cn = "
  //   << chordnormaldir.x_ <<" "<< chordnormaldir.y_<<" "<<chordnormaldir.z_<< " "
  //   <<std::endl;
  alpha = alphaNoTwist + twist;  
}


void 
AirfoilTheory2D::calculate_cl_cd(
    double alpha,
    std::vector<double> aoatable,
    std::vector<double> cltable,
    std::vector<double> cdtable,
    double &cl,
    double &cd)
{
  // Get cl and cd from the tables
  utils::linear_interp(aoatable, cltable, alpha, cl);
  utils::linear_interp(aoatable, cdtable, alpha, cd);

  // Do another other processing needed on cl/cd
  // [..nothing at this time..]
}

ActSimpleSetUpThrustCalc::ActSimpleSetUpThrustCalc(ActuatorBulkSimple& actBulk)
  : actBulk_(actBulk)
{
}

void
ActSimpleSetUpThrustCalc::operator()(int index) const
{
  auto hubLoc = Kokkos::subview(actBulk_.hubLocations_, index, Kokkos::ALL);
  auto hubOri = Kokkos::subview(actBulk_.hubOrientation_, index, Kokkos::ALL);
  auto thrust = Kokkos::subview(actBulk_.turbineThrust_, index, Kokkos::ALL);
  auto torque = Kokkos::subview(actBulk_.turbineTorque_, index, Kokkos::ALL);

  for (int i = 0; i < 3; i++) {
    thrust(i) = 0.0;
    torque(i) = 0.0;
  }

  for (int j = 0; j < 3; j++) {
    hubLoc(j) = 0.0;
    hubOri(j) = 0.0;
  }
  /*
  if (actBulk_.localTurbineId_ == index) {
    actBulk_.openFast_.getHubPos(hubLoc.data(), index);
    actBulk_.openFast_.getHubShftDir(hubOri.data(), index);
  } else {
    for (int j = 0; j < 3; j++) {
      hubLoc(j) = 0.0;
      hubOri(j) = 0.0;
    }
  }
  */
}

void
ActSimpleComputeThrustInnerLoop::operator()(
  const uint64_t pointId,
  const double* nodeCoords,
  double* sourceTerm,
  const double dual_vol,
  const double scvIp) const
{

  auto offsets = actBulk_.turbIdOffset_.view_host();

  // shouldn't thrust and torque contribs only come from blades?
  // probably not worth worrying about since this is just a debug calculation

  // determine turbine
  /*
  int turbId = 0;
  const int nPointId = static_cast<int>(pointId);
  for (; turbId < offsets.extent_int(0); turbId++) {
    if (nPointId >= offsets(turbId)) {
      break;
    }
  }
  */
  if (NaluEnv::self().parallel_rank()<actBulk_.num_blades_) {
    int turbId = NaluEnv::self().parallel_rank();
  //auto hubLoc = Kokkos::subview(actBulk_.hubLocations_, turbId, Kokkos::ALL);
  //auto hubOri = Kokkos::subview(actBulk_.hubOrientation_, turbId, Kokkos::ALL);
  auto thrust = Kokkos::subview(actBulk_.turbineThrust_, turbId, Kokkos::ALL);
  auto torque = Kokkos::subview(actBulk_.turbineTorque_, turbId, Kokkos::ALL);

  double r[3], rPerpShaft[3], forceTerm[3];

  for (int i = 0; i < 3; i++) {
    //forceTerm[i] = sourceTerm[i] * scvIp;
    forceTerm[i] = sourceTerm[i]*scvIp;
    //r[i] = nodeCoords[i] - hubLoc(i);
    thrust(i) += forceTerm[i];
  }

  /*
  double rDotHubOri = 0;
  for (int i = 0; i < 3; i++) {
    rDotHubOri += r[i] * hubOri(i);
  }
  

  for (int i = 0; i < 3; i++) {
    rPerpShaft[i] = r[i] - rDotHubOri * hubOri(i);
  }

  torque(0) += (rPerpShaft[1] * forceTerm[2] - rPerpShaft[2] * forceTerm[1]);
  torque(1) += (rPerpShaft[2] * forceTerm[0] - rPerpShaft[0] * forceTerm[2]);
  torque(2) += (rPerpShaft[0] * forceTerm[1] - rPerpShaft[1] * forceTerm[0]);
  */
  torque(0) = 0.0;
  torque(1) = 0.0;
  torque(2) = 0.0;
  }
}

ActSimpleStashOrientationVectors::ActSimpleStashOrientationVectors(
  ActuatorBulkSimple& actBulk)
  : orientation_(helper_.get_local_view(actBulk.orientationTensor_)),
    offset_(helper_.get_local_view(actBulk.turbIdOffset_)),
    turbId_(actBulk.localTurbineId_),
    fast_(actBulk.openFast_)
{
  helper_.touch_dual_view(actBulk.orientationTensor_);
  actBulk.turbIdOffset_.sync_host();
}

void
ActSimpleStashOrientationVectors::operator()(int index) const
{
  const int pointId = index - offset_(turbId_);
  auto localOrientation = Kokkos::subview(orientation_, index, Kokkos::ALL);
  //if (fast_.getForceNodeType(turbId_, pointId) == fast::BLADE) {
  if(pointId>0){
    fast_.getForceNodeOrientation(localOrientation.data(), pointId, turbId_);

    // swap columns of matrix since openfast stores data
    // as (thick, chord, span) and we want (chord, thick, span)
    double colSwapTemp;
    for (int i = 0; i < 9;i+=3) {
      colSwapTemp = localOrientation(i);
      localOrientation(i) = localOrientation(i+1);
      localOrientation(i+1) = colSwapTemp;
    }
  } else {
    // identity matrix
    // (all other terms should have already been set to zero)
    localOrientation(0) = 1.0;
    localOrientation(3) = 1.0;
    localOrientation(6) = 1.0;
  }

}

void
ActSimpleSpreadForceWhProjInnerLoop::preloop()
{
  actBulk_.actuatorForce_.sync_host();
}

void
ActSimpleSpreadForceWhProjInnerLoop::operator()(
  const uint64_t pointId,
  const double* nodeCoords,
  double* sourceTerm,
  const double dual_vol,
  const double scvIp) const
{

  auto pointCoords =
    Kokkos::subview(actBulk_.pointCentroid_.view_host(), pointId, Kokkos::ALL);

  auto pointForce =
    Kokkos::subview(actBulk_.actuatorForce_.view_host(), pointId, Kokkos::ALL);

  auto epsilon =
    Kokkos::subview(actBulk_.epsilon_.view_host(), pointId, Kokkos::ALL);

  auto orientation = Kokkos::subview(
    actBulk_.orientationTensor_.view_host(), pointId, Kokkos::ALL);

  double distance[3]={0, 0, 0};
  double projectedDistance[3]={0, 0, 0};
  double projectedForce[3]={0, 0, 0};

  actuator_utils::compute_distance(
    3, nodeCoords, pointCoords.data(), &distance[0]);

  // transform distance from Cartesian to blade coordinate system
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      projectedDistance[i] += distance[j] * orientation(i+j*3);
    }
  }

  const double gauss = actuator_utils::Gaussian_projection(
    3, &projectedDistance[0], epsilon.data());

  for (int j = 0; j < 3; j++) {
    projectedForce[j] = gauss * pointForce(j);
  }

  for (int j = 0; j < 3; j++) {
    sourceTerm[j] += projectedForce[j] * scvIp / dual_vol;
  }

}

} /* namespace nalu */
} /* namespace sierra */
