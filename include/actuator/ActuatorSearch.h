// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#ifndef ACTUATORSEARCH_H_
#define ACTUATORSEARCH_H_

#include <actuator/ActuatorTypes.h>
#include <stk_mesh/base/BulkData.hpp>
#include <Kokkos_Core.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_search/CoarseSearch.hpp>

// common type defs
using theKey = stk::search::IdentProc<uint64_t, int>;
using Point = stk::search::Point<double>;
using Sphere = stk::search::Sphere<double>;
using Box = stk::search::Box<double>;
using boundingSphere = std::pair<Sphere, theKey>;
using boundingElementBox = std::pair<Box, theKey>;
using boundingBox = std::pair<Box, theKey>;
using VecBoundSphere = std::vector<boundingSphere>;
using VecBoundElemBox = std::vector<boundingElementBox>;
using VecBoundBox = std::vector<boundingBox>;
using SearchKeyPair = std::pair<theKey, theKey>;
using VecSearchKeyPair = std::vector<SearchKeyPair>;

namespace sierra {
namespace nalu {

VecBoundSphere
create_bounding_spheres(ActFixVectorDbl points, ActFixScalarDbl searchRadius);

VecBoundBox create_bounding_boxes(ActFixVectorDbl points, ActFixVectorDbl dX);

VecBoundElemBox create_element_boxes(
  stk::mesh::BulkData& stkBulk, std::vector<std::string> partNameList);

template<typename T>
void execute_coarse_search(
  T& geoEntities,
  VecBoundElemBox& elemBoxes,
  ActScalarU64Dv& coarsePointIds,
  ActScalarU64Dv& coarseElemIds,
  stk::search::SearchMethod searchMethod)
{
    VecSearchKeyPair searchKeyPair;
  stk::search::coarse_search(
    geoEntities, elemBoxes, searchMethod, MPI_COMM_SELF, searchKeyPair);

  // sort by actuator point id auto can be used with c++ 14
  std::sort(
    searchKeyPair.begin(), searchKeyPair.end(),
    [](const auto& left, const auto& right) {
      return left.second < right.second;
    });

    const std::size_t numLocalMatches = searchKeyPair.size();

  coarsePointIds.resize(numLocalMatches);
  coarseElemIds.resize(numLocalMatches);

  coarsePointIds.modify_host();
  coarseElemIds.modify_host();

  for (std::size_t i = 0; i < numLocalMatches; i++) {
    coarsePointIds.h_view(i) = searchKeyPair[i].first.id();
    coarseElemIds.h_view(i) = searchKeyPair[i].second.id();
  }
}

void execute_fine_search(
  stk::mesh::BulkData& stkBulk,
  ActScalarU64Dv coarsePointIds,
  ActScalarU64Dv coarseElemIds,
  ActFixVectorDbl points,
  ActFixElemIds matchElemIds,
  ActFixVectorDbl localCoords,
  ActFixScalarBool isLocalPoint,
  ActFixScalarInt localParallelRedundancy);

} // namespace nalu
} // namespace sierra

#endif // ACTUATORSEARCH_H_
