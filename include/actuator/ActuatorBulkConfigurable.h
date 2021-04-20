// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//
#ifndef ACTUATORBULKCONFIGURABLE_H_
#define ACTUATORBULKCONFIGURABLE_H_

#include <actuator/ActuatorBulkFAST.h>
#include <actuator/ActuatorSpace.h>

namespace sierra {
namespace nalu {

struct ActuatorBulkConfigurable : public ActuatorBulkFAST
{
  ActuatorBulkConfigurable(
    const ActuatorMetaFAST& actMeta, const double naluTimeStep);
  VecBoundBox discretize_disk_section_with_boxes(
    const double r,
    const double dr,
    const double dTheta,
    const double thetaStart,
    const double zRef,
    const double dz);
  std::vector<actuator::ActuatorSpace> actuatorSpaces_;
};
} // namespace nalu
} // namespace sierra

#endif /* ACTUATORBULKCONFIGURABLE_H_ */
