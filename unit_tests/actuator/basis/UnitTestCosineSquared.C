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
#include <cmath>
#include <actuator/basis/CosineSquared.h>

namespace sierra {
namespace nalu {
namespace actuator {
namespace {

TEST(ActuatorBasisCos2, DISABLED_intgralPeriodic3Points)
{
  const int nActPoints = 3;
  const int nSamplePoints = 1000;
  double dtheta = 2 * M_PI / nActPoints;
  double actPointCoord[nActPoints];
  double integral[nActPoints];
  for (int i = 0; i < nActPoints; i++) {
    actPointCoord[i] = i * dtheta;
    integral[i] = 0.0;
  }
  double dx = 2 * M_PI / (nSamplePoints);

  Cos2Basis b(nActPoints);

  EXPECT_DOUBLE_EQ(
    0.0, b.get_interpolation_weight(&actPointCoord[1], &actPointCoord[2]));
  EXPECT_DOUBLE_EQ(
    0.5 * dtheta,
    b.get_interpolation_weight(&actPointCoord[1], &actPointCoord[1]));

  for (int i = 0; i < nSamplePoints; i++) {
    const double x = i * dx;
    for (int j = 0; j < nActPoints; j++) {
      integral[j] += b.get_interpolation_weight(&actPointCoord[j], &x);
    }
  }

  const double exactIntegral = 1.0 / nActPoints;
  for (int i = 0; i < nActPoints; i++) {
    EXPECT_NEAR(exactIntegral, integral[i] / nSamplePoints, 1e-10)
      << integral[i] / exactIntegral;
  }
}

TEST(ActuatorBasisCos2, DISABLED_intgralNonPeriodic3Points)
{
  const int nActPoints = 3;
  const int nSamplePoints = 10;
  double dtheta = 1.0 / (nActPoints - 1);
  double actPointCoord[nActPoints];
  double integral[nActPoints];
  for (int i = 0; i < nActPoints; i++) {
    actPointCoord[i] = i * dtheta;
    integral[i] = 0.0;
  }
  double dx = 1.0 / (nSamplePoints - 1);

  Cos2Basis b(nActPoints);
  for (int i = 0; i < nSamplePoints; i++) {
    const double x = i * dx;
    for (int j = 0; j < nActPoints; j++) {
      integral[j] += b.get_interpolation_weight(&actPointCoord[j], &x);
    }
  }

  const double exactIntegral = 1.0 / nActPoints;
  for (int i = 0; i < nActPoints; i++) {
    EXPECT_NEAR(exactIntegral, integral[i], 1e-10)
      << integral[i] / exactIntegral;
  }
}
} // namespace
} // namespace actuator
} // namespace nalu
} // namespace sierra