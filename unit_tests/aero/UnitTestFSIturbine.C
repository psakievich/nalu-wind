#include "gtest/gtest.h"

#include "stk_util/environment/WallTime.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/NgpField.hpp"
#include "stk_mesh/base/GetNgpField.hpp"

#include "UnitTestUtils.h"

#include "KokkosInterface.h"
#include "aero/fsi/FSIturbine.h"

namespace {

// using block_1 and surface_1 for everything is no doubt silly...
// as we ramp up the FSI testing we will no doubt need to use
// more sophisticated input. But this gets us off the
// ground, so to speak.
const std::string fsiInputs = "tower_parts: [block_1] \n"
                              "hub_parts: [block_1]\n"
                              "nacelle_parts: [block_1] \n"
                              "blade_parts:\n"
                              "  - [block_1]\n"
                              "  - [block_1]\n"
                              "deflection_ramping:\n"
                              "  span_ramp_distance: 10.0\n"
                              "  zero_theta_ramp_angle: 180.0\n"
                              "  theta_ramp_span: 15.0\n"
                              "  temporal_ramp_start: 0\n"
                              "  temporal_ramp_end: 10\n"
                              "tower_boundary_parts: [surface_1] \n"
                              "hub_boundary_parts: [surface_1]\n"
                              "nacelle_boundary_parts: [surface_1] \n"
                              "blade_boundary_parts:\n"
                              "  - [surface_1]\n";

YAML::Node
create_fsi_yaml_node()
{
  YAML::Node yamlNode = YAML::Load(fsiInputs);
  return yamlNode;
}

TEST_F(CylinderMesh, construct_FSIturbine)
{
  const double innerRadius = 1.0;
  const double outerRadius = 2.0;
  fill_mesh_and_initialize_test_fields(20, 20, 20, innerRadius, outerRadius);

  YAML::Node yamlNode = create_fsi_yaml_node();
  sierra::nalu::fsiTurbine fsiTurb(0, yamlNode);
  EXPECT_NO_THROW(fsiTurb.setup(bulk));
}

} // namespace
