// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef INCLUDE_ACTUATOR_ACTUATORTYPES_H_
#define INCLUDE_ACTUATOR_ACTUATORTYPES_H_

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

namespace sierra{
namespace nalu{


#ifdef ACTUATOR_ON_DEVICE
using ActuatorMemSpace = Kokkos::CudaSpace;
using ActuatorMemLayout = Kokkos::LayoutRight;
using ActuatorExecutionSpace = Kokkos::DefaultExecutionSpace;
#else
using ActuatorMemSpace = Kokkos::HostSpace;
using ActuatorMemLayout = Kokkos::LayoutLeft;
using ActuatorExecutionSpace = Kokkos::DefaultHostExecutionSpace;
#endif
using ActuatorFixedMemSpace = Kokkos::HostSpace;
using ActuatorFixedMemLayout = Kokkos::LayoutLeft;
using ActuatorFixedExecutionSpace = Kokkos::DefaultHostExecutionSpace;


using ActScalarIntDv = Kokkos::DualView<int*,     ActuatorMemLayout, ActuatorMemSpace>;
using ActScalarDblDv = Kokkos::DualView<double*,  ActuatorMemLayout, ActuatorMemSpace>;
using ActVectorDblDv = Kokkos::DualView<double**, ActuatorMemLayout, ActuatorMemSpace>;

using ActFixRangePolicy = Kokkos::RangePolicy<ActuatorFixedExecutionSpace>;
using ActFixScalarInt = Kokkos::View<int*,     ActuatorFixedMemLayout, ActuatorFixedMemSpace>;
using ActFixScalarDbl = Kokkos::View<double*,  ActuatorFixedMemLayout, ActuatorFixedMemSpace>;
using ActFixVectorDbl = Kokkos::View<double*[3], ActuatorFixedMemLayout, ActuatorFixedMemSpace>;
using ActFixElemIds   = Kokkos::View<uint64_t*, ActuatorFixedMemLayout, ActuatorFixedMemSpace>;
} //namespace nalu
} //namespace sierra

#endif // INCLUDE_ACTUATOR_ACTUATORTYPES_H_
