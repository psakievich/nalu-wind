#ifndef TIOGASTKIFACE_H
#define TIOGASTKIFACE_H

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include "overset/TiogaOptions.h"
#include "overset/OversetFieldData.h"

#include <vector>
#include <memory>
#include <array>

namespace YAML {
class Node;
}

namespace TIOGA {
class tioga;
}

namespace sierra {
namespace nalu {

class OversetManagerTIOGA;
}
} // namespace sierra

namespace tioga_nalu {

class TiogaBlock;

/** Nalu interface to TIOGA (Topology Independent Overset Grid Assembly)
 *
 *  This class provides a two-way data transfer interface for TIOGA library and
 *  provides overset connectivity capability for Nalu.
 */
class TiogaSTKIface
{
public:
  /**
   *  @param oversetManager Reference to Nalu OversetManager object
   *  @param node YAML node containing overset inputs
   */
  TiogaSTKIface(
    sierra::nalu::OversetManagerTIOGA&, const YAML::Node&, const std::string&);

  ~TiogaSTKIface();

  /** Setup block structure information (steps before mesh creation)
   */
  void setup(stk::mesh::PartVector&);

  /** Initialize mesh data structure (steps after mesh creation)
   */
  void initialize();

  /** Determine overset connectivity by calling into TIOGA API
   *
   *  This method performs several steps: updates coordinates (if necessary,
   *  during mesh motion), registers the mesh blocks to TIOGA, calculate mesh
   *  connectivity information (hole, fringe, and field point determination),
   *  update the "overset inactive part" for hole elements, create the {fringe
   *  node, donor element} mapping pair data structures for overset simulations.
   */
  void execute(const bool isDecoupled);

  void register_mesh();

  void post_connectivity_work(const bool isDecoupled = true);

  int register_solution(const std::vector<sierra::nalu::OversetFieldData>&);

  void update_solution(const std::vector<sierra::nalu::OversetFieldData>&);

  virtual void
  overset_update_fields(const std::vector<sierra::nalu::OversetFieldData>&);

  virtual void overset_update_field(
    stk::mesh::FieldBase* field,
    const int nrows = 1,
    const int ncols = 1,
    const bool doFinalSyncToDevice = true);

private:
  TiogaSTKIface() = delete;
  TiogaSTKIface(const TiogaSTKIface&) = delete;

  /** Process the input parameters and initialize all data structures necessary
   * to call TIOGA.
   */
  void load(const YAML::Node&);

  /** Ghost donor elements to receptor MPI ranks
   */
  void update_ghosting();

  /** Reset all connectivity data structures when recomputing connectivity
   */
  void reset_data_structures();

  /** Gather receptor-donor pair information
   *
   *  Extract the set of fringe nodes on participating meshes and prepare for
   *  reconciliation for shared nodes.
   */
  void get_receptor_info();

  /** Populate the overset info data structure
   *
   *  Reconcile shared/owned fringe status and populate the overset fringe
   *  vector for later use by algorithms.
   */
  void populate_overset_info();

  //! Synchronize fields before performing overset connectivity
  void pre_connectivity_sync();

  //! Synchronize modified fields after performing overset connectivity
  void post_connectivity_sync();

  //! Reference to Nalu OversetManager object
  sierra::nalu::OversetManagerTIOGA& oversetManager_;

  //! Reference to the STK MetaData object
  stk::mesh::MetaData& meta_;

  //! Reference to the STK BulkData object
  stk::mesh::BulkData& bulk_;

  TiogaOptions tiogaOpts_;

  //! List of TIOGA data structures for each mesh block participating in overset
  //! connectivity
  std::vector<std::unique_ptr<TiogaBlock>> blocks_;

  //! Reference to the TIOGA API interface
  TIOGA::tioga& tg_;

  //! Work array used to hold donor elements that require ghosting to receptor
  //! MPI ranks
  stk::mesh::EntityProcVec elemsToGhost_;

  //! List of receptor nodes that are shared entities across MPI ranks. This
  //! information is used to synchronize the field vs. fringe point status for
  //! these shared nodes across processor boundaries.
  std::vector<stk::mesh::EntityId> receptorIDs_;

  //! Donor elements corresponding to TiogaSTKIface::receptorIDs_ that must be
  //! ghosted to another MPI rank to ensure that owned and shared nodes are
  //! consistent.
  std::vector<stk::mesh::EntityId> donorIDs_;

  //! Name of the coordinates field (for moving mesh simulations)
  std::string coordsName_;
};

} // namespace tioga_nalu

#endif /* TIOGASTKIFACE_H */
