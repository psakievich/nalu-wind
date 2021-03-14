// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef COSINESQUARED_H_
#define COSINESQUARED_H_

namespace sierra {
namespace nalu {
namespace actuator {
class Cos2Basis
{
private:
  const bool periodic_;
  const double length_;
  const double dX_;

public:
  Cos2Basis(const bool periodic);
  double get_interpolation_weight(
    const double* actPointCoord, const double* sampleCoord);
};
} // namespace actuator
} // namespace nalu
} // namespace sierra

#endif /* COSINESQUARED_H_ */
