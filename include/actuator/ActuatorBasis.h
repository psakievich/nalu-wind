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

/**
 * @brief Interface class for all bases using CRTP
 *
 * Thinking CRTP can be used for a pointer since template combinatorics may
 * explode as we add more basis. Right now CRTP not seeming necessary so may
 * remove with a refactor.
 */
template <typename Derived>
class Basis
{
public:
  inline double get_interpolation_weight(
    const double* actPointCoord, const double* sampleCoord)
  {
    return static_cast<Derived*>(this)->get_interpolation_weight_impl(
      actPointCoord, sampleCoord);
  };
};

class NullBasis : public Basis<NullBasis>
{
public:
  NullBasis() = default;

  inline double get_interpolation_weight_impl(const double*, const double*)
  {
    return 1.0;
  }
};

} // namespace actuator
} // namespace nalu
} // namespace sierra

#endif /* ACTUATORBASIS_H_ */
