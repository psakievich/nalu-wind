// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//
#include <actuator/ActuatorSpace.h>
#include <actuator/basis/Gaussian.h>

namespace sierra {
namespace nalu {
namespace actuator {

double
ActuatorSpace::get_interpolation_weight(
  const double* actPointCoord, const double* sampleCoord)
{
  double value{1.0};
  for (auto&& b : bases_) {
    value *= b->get_interpolation_weight(actPointCoord, sampleCoord);
  }
  return value;
}
} // namespace actuator
} // namespace nalu
} // namespace sierra