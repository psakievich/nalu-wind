// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//
#include "gtest/gtest.h"
#include <stk_math/StkMath.hpp>
#include <actuator/basis/Gaussian.h>

namespace sierra {
namespace nalu {
namespace actuator {

namespace {

TEST(ActuatorIsotropicGaussian, NGP_unityGaussianTest)
{
  const double epsilon = stk::math::pow(M_PI, -0.5);
  for (int i = 1; i <= 3; i++) {
    IsotropicGaussianBasis b(i, epsilon);
    const std::vector<double> actPointCoord(i, 1.0);
    const std::vector<double> sampleCoord(i, 1.0);
    EXPECT_DOUBLE_EQ(
      1.0, b.get_interpolation_weight(actPointCoord.data(), sampleCoord.data()))
      << "ndim: " << i;
  }
}

TEST(ActuatorIsotropicGaussian, NGP_gaussianSampling)
{
  const double epsilon = stk::math::pow(M_PI, -0.5);
  IsotropicGaussianBasis b(1, epsilon);
  const double actPointCoord = 0.0;
  const int N = 1000;
  const double dX = 2.0 / (N - 1);
  double integral = 0.0;
  for (int i = 0; i < N; ++i) {
    const double x = -1.0 + i * dX;
    integral += b.get_interpolation_weight(&actPointCoord, &x) * dX;
  }
  EXPECT_NEAR(1.0, integral, 1e-1);
}

} // namespace

} // namespace actuator
} // namespace nalu
} // namespace sierra
