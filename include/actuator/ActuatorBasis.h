// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef ACTUATORBASIS_H_
#define ACTUATORBASIS_H_

namespace sierra {
namespace nalu {
namespace actuator {

class Basis
{
public:
  virtual inline double get_interpolation_weight(
    const double* actPointCoord, const double* sampleCoord) = 0;
  virtual ~Basis() = default;
};

class NullBasis : public Basis
{
public:
  NullBasis() = default;

  inline double get_interpolation_weight(const double*, const double*) override
  {
    return 1.0;
  }
};

} // namespace actuator
} // namespace nalu
} // namespace sierra

#endif /* ACTUATORBASIS_H_ */
