// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include <GammaEquationSystem.h>
#include <AlgorithmDriver.h>
#include <AssembleScalarNonConformalSolverAlgorithm.h>
#include <AssembleNodeSolverAlgorithm.h>
#include <AssembleNodalGradNonConformalAlgorithm.h>
#include <AuxFunctionAlgorithm.h>
#include <ConstantAuxFunction.h>
#include <CopyFieldAlgorithm.h>
#include <DirichletBC.h>
#include <EquationSystem.h>
#include <EquationSystems.h>
#include <Enums.h>
#include <FieldFunctions.h>
#include <LinearSolvers.h>
#include <LinearSolver.h>
#include <LinearSystem.h>
#include <NaluEnv.h>
#include <NaluParsing.h>
#include <Realm.h>
#include <Realms.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <SolverAlgorithmDriver.h>

// template for supp algs
#include <AlgTraits.h>
#include <kernel/KernelBuilder.h>
#include <kernel/KernelBuilderLog.h>

// consolidated
#include <AssembleElemSolverAlgorithm.h>

// edge kernels
#include <edge_kernels/ScalarEdgeSolverAlg.h>
#include <edge_kernels/ScalarOpenEdgeKernel.h>

// node kernels
#include <node_kernels/NodeKernelUtils.h>
#include <node_kernels/ScalarMassBDFNodeKernel.h>
#include <node_kernels/ScalarGclNodeKernel.h>
#include <node_kernels/BLTGammaM2015NodeKernel.h>

// ngp
#include "ngp_utils/NgpFieldBLAS.h"
#include "ngp_algorithms/NodalGradEdgeAlg.h"
#include "ngp_algorithms/NodalGradElemAlg.h"
#include "ngp_algorithms/NodalGradBndryElemAlg.h"
#include "ngp_algorithms/EffSSTDiffFluxCoeffAlg.h"
#include "utils/StkHelpers.h"

#include <overset/UpdateOversetFringeAlgorithmDriver.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra {
namespace nalu {

//==========================================================================
// Class Definition
//==========================================================================
// GammaEquationSystem - manages Gamma equation in SST transition model
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
GammaEquationSystem::GammaEquationSystem(EquationSystems& eqSystems)
  : EquationSystem(eqSystems, "GammaEQS", "gamma_transition"),
    managePNG_(realm_.get_consistent_mass_matrix_png("gamma_transition")),
    gamma_(NULL),
    dgamdx_(NULL),
    gamTmp_(NULL),
    minDistanceToWall_(NULL),
    visc_(NULL),
    tvisc_(NULL),
    evisc_(NULL),
    nodalGradAlgDriver_(realm_, "dgamdx")
{
  dofName_ = "gamma_transition";

  // extract solver name and solver object
  std::string solverName =
    realm_.equationSystems_.get_solver_block_name("gamma_transition");
  LinearSolver* solver = realm_.root()->linearSolvers_->create_solver(
    solverName, realm_.name(), EQ_GAMMA_TRANS);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // determine nodal gradient form
  set_nodal_gradient("gamma_transition");
  NaluEnv::self().naluOutputP0()
    << "Edge projected nodal gradient for gamma_transition: "
    << edgeNodalGradient_ << std::endl;

  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // create projected nodal gradient equation system
  if (managePNG_)
    throw std::runtime_error(
      "GammaEquationSystem::Error managePNG is not complete");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
GammaEquationSystem::~GammaEquationSystem() = default;

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::register_nodal_fields(stk::mesh::Part* part)
{

  stk::mesh::MetaData& meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  gamma_ = &(meta_data.declare_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "gamma_transition", numStates));
  stk::mesh::put_field_on_mesh(*gamma_, *part, nullptr);
  realm_.augment_restart_variable_list("gamma_transition");

  dgamdx_ = &(meta_data.declare_field<VectorFieldType>(
    stk::topology::NODE_RANK, "dgamdx"));
  stk::mesh::put_field_on_mesh(*dgamdx_, *part, nDim, nullptr);

  // delta solution for linear solver; share delta since this is a split system
  gamTmp_ = &(meta_data.declare_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "gamTmp"));
  stk::mesh::put_field_on_mesh(*gamTmp_, *part, nullptr);

  visc_ = &(meta_data.declare_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field_on_mesh(*visc_, *part, nullptr);

  tvisc_ = &(meta_data.declare_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "turbulent_viscosity"));
  stk::mesh::put_field_on_mesh(*tvisc_, *part, nullptr);

  evisc_ = &(meta_data.declare_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "effective_viscosity_gamma"));
  stk::mesh::put_field_on_mesh(*evisc_, *part, nullptr);

  // make sure all states are properly populated (restart can handle this)
  if (
    numStates > 2 &&
    (!realm_.restarted_simulation() || realm_.support_inconsistent_restart())) {
    ScalarFieldType& gammaN = gamma_->field_of_state(stk::mesh::StateN);
    ScalarFieldType& gammaNp1 = gamma_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm* theCopyAlg = new CopyFieldAlgorithm(
      realm_, part, &gammaNp1, &gammaN, 0, 1, stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::register_interior_algorithm(stk::mesh::Part* part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType& gammaNp1 = gamma_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType& dgamdxNone = dgamdx_->field_of_state(stk::mesh::StateNone);

  if (edgeNodalGradient_ && realm_.realmUsesEdges_)
    nodalGradAlgDriver_.register_edge_algorithm<ScalarNodalGradEdgeAlg>(
      algType, part, "gamma_nodal_grad", &gammaNp1, &dgamdxNone);
  else
    nodalGradAlgDriver_.register_elem_algorithm<ScalarNodalGradElemAlg>(
      algType, part, "gamma_nodal_grad", &gammaNp1, &dgamdxNone,
      edgeNodalGradient_);

  // solver; interior contribution (advection + diffusion)
  if (!realm_.solutionOptions_->useConsolidatedSolverAlg_) {

    std::map<AlgorithmType, SolverAlgorithm*>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if (itsi == solverAlgDriver_->solverAlgMap_.end()) {
      SolverAlgorithm* theAlg = NULL;
      if (realm_.realmUsesEdges_) {
        const bool useAvgMdot = (realm_.solutionOptions_->turbulenceModel_ ==
                                 TurbulenceModel::SST_AMS)
                                  ? true
                                  : false;
        theAlg = new ScalarEdgeSolverAlg(
          realm_, part, this, gamma_, dgamdx_, evisc_, useAvgMdot);
      } else {
        throw std::runtime_error(
          "GAMMAEQS: Attempt to use non-NGP element solver algorithm");
      }
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;

      // look for fully integrated source terms
      std::map<std::string, std::vector<std::string>>::iterator isrc =
        realm_.solutionOptions_->elemSrcTermsMap_.find("gamma_transition");
      if (isrc != realm_.solutionOptions_->elemSrcTermsMap_.end()) {
        throw std::runtime_error(
          "GammaElemSrcTerms::Error can not use element source "
          "terms for an edge-based scheme");
      }
    } else {
      itsi->second->partVec_.push_back(part);
    }

    // Check if the user has requested CMM or LMM algorithms; if so, do not
    // include Nodal Mass algorithms
    //
    NaluEnv::self().naluOutputP0() << "register gamma interior: " << std::endl;

    std::vector<std::string> checkAlgNames = {
      "gamma_transition_time_derivative",
      "lumped_gamma_transition_time_derivative"};
    bool elementMassAlg = supp_alg_is_requested(checkAlgNames);
    if (elementMassAlg) {
      throw std::runtime_error(
        "consistent mass integration of gamma time-derivative unavailable");
    }

    auto& solverAlgMap = solverAlgDriver_->solverAlgMap_;
    process_ngp_node_kernels(
      solverAlgMap, realm_, part, this,
      [&](AssembleNGPNodeSolverAlgorithm& nodeAlg) {
        nodeAlg.add_kernel<ScalarMassBDFNodeKernel>(realm_.bulk_data(), gamma_);

        NaluEnv::self().naluOutputP0()
          << "call BLTGammaM2015NodeKernel: " << std::endl;

        nodeAlg.add_kernel<BLTGammaM2015NodeKernel>(realm_.meta_data());
      },
      [&](AssembleNGPNodeSolverAlgorithm& nodeAlg, std::string& srcName) {
        if (srcName == "gcl") {
          nodeAlg.add_kernel<ScalarGclNodeKernel>(realm_.bulk_data(), gamma_);
          NaluEnv::self().naluOutputP0() << " - " << srcName << std::endl;
        } else
          throw std::runtime_error("SDREqSys: Invalid source term: " + srcName);
      });
  } else {
    throw std::runtime_error("GAMMAEQS: Element terms not supported");
  }

  // effective diffusive flux coefficient alg
  if (!effDiffFluxAlg_) {
    effDiffFluxAlg_.reset(new EffSSTDiffFluxCoeffAlg(
      realm_, part, visc_, tvisc_, evisc_, 1.0, 1.0));
  } else {
    effDiffFluxAlg_->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::register_inflow_bc(
  stk::mesh::Part* part,
  const stk::topology& /*theTopo*/,
  const InflowBoundaryConditionData& inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType& gammaNp1 = gamma_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType& dgamdxNone = dgamdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData& meta_data = realm_.meta_data();

  // register boundary data; gamma_bc
  ScalarFieldType* theBcField = &(meta_data.declare_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "gamma_bc"));
  stk::mesh::put_field_on_mesh(*theBcField, *part, nullptr);

  // extract the value for user specified tke and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  GammaInf gamma = userData.gamma_;
  std::vector<double> userSpec(1);
  userSpec[0] = gamma.gamma_;

  // new it
  ConstantAuxFunction* theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm* auxAlg = new AuxFunctionAlgorithm(
    realm_, part, theBcField, theAuxFunc, stk::topology::NODE_RANK);

  // how to populate the field?
  if (userData.externalData_) {
    // xfer will handle population; only need to populate the initial value
    realm_.initCondAlg_.push_back(auxAlg);
  } else {
    // put it on bcData
    bcDataAlg_.push_back(auxAlg);
  }

  // copy gamma_bc to gamma_transition np1...
  CopyFieldAlgorithm* theCopyAlg = new CopyFieldAlgorithm(
    realm_, part, theBcField, &gammaNp1, 0, 1, stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg);

  // non-solver; dgamdx; allow for element-based shifted
  nodalGradAlgDriver_.register_face_algorithm<ScalarNodalGradBndryElemAlg>(
    algType, part, "gamma_nodal_grad", &gammaNp1, &dgamdxNone,
    edgeNodalGradient_);

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm*>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if (itd == solverAlgDriver_->solverDirichAlgMap_.end()) {
    DirichletBC* theAlg =
      new DirichletBC(realm_, this, part, &gammaNp1, theBcField, 0, 1);
    solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
  } else {
    itd->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::register_open_bc(
  stk::mesh::Part* part,
  const stk::topology& /* partTopo */,
  const OpenBoundaryConditionData& /* openBCData */)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType& gammaNp1 = gamma_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType& dgamdxNone = dgamdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dgamdx; allow for element-based shifted
  nodalGradAlgDriver_.register_face_algorithm<ScalarNodalGradBndryElemAlg>(
    algType, part, "gamma_nodal_grad", &gammaNp1, &dgamdxNone,
    edgeNodalGradient_);
}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::register_wall_bc(
  stk::mesh::Part* part,
  const stk::topology& /*theTopo*/,
  const WallBoundaryConditionData& /* wallBCData */)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType& gammaNp1 = gamma_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType& dgamdxNone = dgamdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dgamdx; allow for element-based shifted
  nodalGradAlgDriver_.register_face_algorithm<ScalarNodalGradBndryElemAlg>(
    algType, part, "gamma_nodal_grad", &gammaNp1, &dgamdxNone,
    edgeNodalGradient_);
}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::register_symmetry_bc(
  stk::mesh::Part* part,
  const stk::topology& /*theTopo*/,
  const SymmetryBoundaryConditionData& /* symmetryBCData */)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  // np1
  ScalarFieldType& gammaNp1 = gamma_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType& dgamdxNone = dgamdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dgamdx; allow for element-based shifted
  nodalGradAlgDriver_.register_face_algorithm<ScalarNodalGradBndryElemAlg>(
    algType, part, "gamma_nodal_grad", &gammaNp1, &dgamdxNone,
    edgeNodalGradient_);
}

//--------------------------------------------------------------------------
//-------- register_overset_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(gamma_);

  equationSystems_.register_overset_field_update(gamma_, 1, 1);
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::reinitialize_linear_system()
{
  // If this is decoupled overset simulation and the user has requested that the
  // linear system be reused, then do nothing
  if (decoupledOverset_ && linsys_->config().reuseLinSysIfPossible())
    return;

  // delete linsys
  delete linsys_;

  // create new solver
  std::string solverName =
    realm_.equationSystems_.get_solver_block_name("gamma_transition");
  LinearSolver* solver = realm_.root()->linearSolvers_->reinitialize_solver(
    solverName, realm_.name(), EQ_GAMMA_TRANS);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- assemble_nodal_gradient() ---------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::assemble_nodal_gradient()
{
  const double timeA = -NaluEnv::self().nalu_time();
  nodalGradAlgDriver_.execute();
  timerMisc_ += (NaluEnv::self().nalu_time() + timeA);
}

//--------------------------------------------------------------------------
//-------- compute_effective_flux_coeff() ----------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::compute_effective_diff_flux_coeff()
{
  const double timeA = -NaluEnv::self().nalu_time();
  effDiffFluxAlg_->execute();
  timerMisc_ += (NaluEnv::self().nalu_time() + timeA);
}

//--------------------------------------------------------------------------
//-------- predict_state() -------------------------------------------------
//--------------------------------------------------------------------------
void
GammaEquationSystem::predict_state()
{
  const auto& ngpMesh = realm_.ngp_mesh();
  const auto& fieldMgr = realm_.ngp_field_manager();
  const auto& gammaN = fieldMgr.get_field<double>(
    gamma_->field_of_state(stk::mesh::StateN).mesh_meta_data_ordinal());
  auto& gammaNp1 = fieldMgr.get_field<double>(
    gamma_->field_of_state(stk::mesh::StateNP1).mesh_meta_data_ordinal());

  const auto& meta = realm_.meta_data();
  const stk::mesh::Selector sel =
    (meta.locally_owned_part() | meta.globally_shared_part() |
     meta.aura_part()) &
    stk::mesh::selectField(*gamma_);
  nalu_ngp::field_copy(ngpMesh, sel, gammaNp1, gammaN);
  gammaNp1.modify_on_device();
}

} // namespace nalu
} // namespace sierra
