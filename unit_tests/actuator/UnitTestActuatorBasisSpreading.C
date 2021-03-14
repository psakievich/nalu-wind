// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include <actuator/ActuatorFunctors.h>
#include <actuator/ActuatorParsing.h>
#include <actuator/ActuatorInfo.h>
#include <actuator/UtilitiesActuator.h>
#include <UnitTestUtils.h>
#include <yaml-cpp/yaml.h>
#include <gtest/gtest.h>

namespace sierra {
namespace nalu {

/*
Current idea is to use a GenericLoopOverCoarseSearch and the inner loop functor
will be constructed using the spreading basis

so SpreadingFunctor(rBasis,thetaBasis,zBasis)
where the basis can be of specific types

r - Gaussian and linear
theta - Gaussian, Fourier, cos2
z - Guassian and linear

so the basis will need to also have it's own data for caching the interpolation
since this is a stationary disk althought that data will match the coarse
search...

association of search is currently elem a is within search tolerance of actuator
point b

so we can also add neighbors to b to get the bounding box of the actuator
space's discretization infact we may need to redefine the search

*/
TEST(ActuatorBasis1D, )

} // namespace nalu
} // namespace sierra