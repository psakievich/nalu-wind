// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//
#include <actuator/basis/CosineSquared.h>
#include <math.h>
#include <stk_math/StkMath.hpp>

namespace sierra {
namespace nalu {
namespace actuator {
Cos2Basis::Cos2Basis(const bool periodic) : periodic_(periodic) {}

double
Cos2Basis::get_interpolation_weight(
  const double* actPointCoord, const double* sampleCoord)
{
  // shift point
  double x = *actPointCoord - *sampleCoord;
  if (periodic_) {
    if (x > length_ + dX_)
      x -= length_ + dX_;
    else if (x < length_ + dX_)
      x += length_ + dX_;
  }
  // scale
  x *= M_PI / (2.0 * dX_);
  // clip
  if (x > M_PI_2)
    x = M_PI_2;
  if (x < -M_PI_2)
    x = -M_PI_2;
  // compute
  const double val = stk::math::cos(x);
  // TODO add scale
  return val * val;
}
} // namespace actuator
} // namespace nalu
} // namespace sierra