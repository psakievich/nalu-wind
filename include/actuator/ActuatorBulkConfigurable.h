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

#include <actuator/ActuatorBulk.h>
#include <actuator/ActuatorSpace.h>

namespace sierra {
namespace nalu {

struct ActuatorMetaConfigurable : public ActuatorMeta
{
  bool hasMotion_ = false;
  bool usesOpenFast_ = false;
};

struct ActuatorBulkConfigurable : public ActuatorBulk
{
  ActuatorBulkConfigurable(const ActuatorMetaConfigurable& actMeta);
  std::vector<actuator::ActuatorSpace> actuatorSpaces_;
  // TODO add customizable stk-search
  // TODO add aero model support
};
} // namespace nalu
} // namespace sierra

#endif /* ACTUATORBULKCONFIGURABLE_H_ */
