// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include "gtest/gtest.h"
#include <actuator/ActuatorSpace.h>
#include <actuator/ActuatorBasis.h>

namespace sierra {
namespace nalu {
namespace {

TEST(ActuatorSpaceTest, construction)
{
  actuator::NullBasis b;
  actuator::ActuatorSpace s;
  s.add_basis<actuator::NullBasis>(b);
  double x = 1.0;
  EXPECT_EQ(1.0, s.get_interpolation_weight(&x, &x));
  s.add_basis<actuator::NullBasis>(b);
  EXPECT_EQ(1.0, s.get_interpolation_weight(&x, &x));
  s.add_basis<actuator::NullBasis>(b);
  EXPECT_EQ(1.0, s.get_interpolation_weight(&x, &x));
}
} // namespace
} // namespace nalu
} // namespace sierra