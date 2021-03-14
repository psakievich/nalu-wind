// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include <actuator/basis/Gaussian.h>
#include <math.h>
#include <stk_math/StkMath.hpp>
namespace sierra {
namespace nalu {
namespace actuator {

IsotropicGaussianBasis::IsotropicGaussianBasis(
  const int dimension, const double epsilon)
  : nDim_(dimension),
    denominator_(stk::math::pow(epsilon * stk::math::sqrt(M_PI), nDim_)),
    epsilon2_(epsilon * epsilon)
{
}
double
IsotropicGaussianBasis::get_interpolation_weight_impl(
  const double* actPointCoord, const double* sampleCoord)
{
  double distance = 0;

  for (int i = 0; i < nDim_; i++) {
    distance +=
      (actPointCoord[i] - sampleCoord[i]) * (actPointCoord[i] - sampleCoord[i]);
  }

  return stk::math::exp(-distance / epsilon2_) / denominator_;
}

} // namespace actuator
} // namespace nalu
} // namespace sierra