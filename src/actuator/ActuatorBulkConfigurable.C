// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include <actuator/ActuatorBulkConfigurable.h>

namespace sierra {
namespace nalu {
ActuatorBulkConfigurable::ActuatorBulkConfigurable(
  const ActuatorMetaFAST& actMeta, const double naluTimeStep)
  : ActuatorBulkFAST(actMeta, naluTimeStep)
{
  if (actMeta.numberOfActuators_ != 1)
    throw std::runtime_error(
      "ActConfigurable only supports 1 actuator for now");
}

// TODO need to account for normal if turbine isn't lined up with z
VecBoundBox
ActuatorBulkConfigurable::discretize_disk_section_with_boxes(
  const double r,
  const double dr,
  const double dTheta,
  const double thetaStart,
  const double zRef,
  const double dz)
{
  const int nBoxes = std::round(dTheta * r / dr);
  ActFixVectorDbl points("the points", nBoxes);
  ActFixVectorDbl dX("spacing", nBoxes);
  const double deltaT = dTheta / nBoxes;
  for (int i = 0; i < nBoxes; i++) {
    points(i, 0) = r * std::cos(thetaStart + deltaT * i);
    points(i, 1) = r * std::sin(thetaStart + deltaT * i);
    points(i, 2) = zRef;
    dX(i, 0) = dr;
    dX(i, 1) = dr;
    dX(i, 2) = dz;
  }
  return create_bounding_boxes(points, dX);
}

} // namespace nalu
} // namespace sierra