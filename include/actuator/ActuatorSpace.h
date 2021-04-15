// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//
#ifndef ACTUATORSPACE_H_
#define ACTUATORSPACE_H_

#include <actuator/ActuatorBasis.h>
#include <memory>
#include <vector>
namespace sierra {
namespace nalu {
namespace actuator {

// use polymorphism for development. probably should figure out how to do this
// with templates later for gpu compatibility
class ActuatorSpace
{
public:
  ActuatorSpace(){};
  template <typename T>
  inline void add_basis(T b)
  {
    // add via copy constructor
    bases_.push_back(std::make_unique<T>(b));
  }
  double get_interpolation_weight(
    const double* actPointCoord, const double* sampleCoord);

private:
  std::vector<std::unique_ptr<Basis>> bases_;
};
} // namespace actuator
} // namespace nalu
} // namespace sierra

#endif /* ACTUATORSPACE_H_ */
