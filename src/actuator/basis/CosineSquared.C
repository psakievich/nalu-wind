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
#include <iostream>

namespace sierra {
namespace nalu {
namespace actuator {
Cos2Basis::Cos2Basis(const int N) : n_(N), dX_(2 * M_PI / N) {}

double
Cos2Basis::get_interpolation_weight(
  const double* actPointCoord, const double* sampleCoord)
{
  // shift point
  double x = *actPointCoord - *sampleCoord;
  if (x > M_PI)
    x -= 2 * M_PI;
  if (x < -M_PI)
    x += 2 * M_PI;
  // scale
  x *= M_PI / (2.0 * dX_); // x*N/4 for periodic
  // clip
  if (x > dX_ || x < -dX_)
    return 0.0;
  // compute
  const double val = stk::math::cos(x);
  // normalizes the integral dx = 2 pi /N
  return val * val / dX_;
}
} // namespace actuator
} // namespace nalu
} // namespace sierra