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

namespace sierra {
namespace nalu {

class ActuatorSpaceContainer
{
  virtual ~ActuatorSpaceContainer() = default;
};

template <typename...>
class ActuatorSpace;

template <typename B0, typename B1, typename B2>
class ActuatorSpace<B0, B1, B2>
{
public:
  ActuatorSpace(B0 b0, B1 b1, B2 b2) : basis0_(b0), basis1_(b1), basis2_(b2) {}
  inline double get_interpolation_weight(
    const double* actPointCoord, const double* sampleCoord)
  {
    return basis0_.get_interpolation_weight(actPointCoord, sampleCoord) *
           basis1_.get_interpolation_weight(actPointCoord, sampleCoord) *
           basis2_.get_interpolation_weight(actPointCoord, sampleCoord);
  }

private:
  B0 basis0_;
  B1 basis1_;
  B2 basis2_;
};

template <typename B0, typename B1>
class ActuatorSpace<B0, B1>
{
public:
  ActuatorSpace(B0 b0, B1 b1) : basis0_(b0), basis1_(b1) {}
  inline double get_interpolation_weight(
    const double* actPointCoord, const double* sampleCoord)
  {
    return basis0_.get_interpolation_weight(actPointCoord, sampleCoord) *
           basis1_.get_interpolation_weight(actPointCoord, sampleCoord);
  }

private:
  B0 basis0_;
  B1 basis1_;
};

template <typename B0>
class ActuatorSpace<B0>
{
public:
  ActuatorSpace(B0 b0) : basis0_(b0) {}
  inline double get_interpolation_weight(
    const double* actPointCoord, const double* sampleCoord)
  {
    return basis0_.get_interpolation_weight(actPointCoord, sampleCoord);
  }

private:
  B0 basis0_;
};

} // namespace nalu
} // namespace sierra

#endif /* ACTUATORSPACE_H_ */
