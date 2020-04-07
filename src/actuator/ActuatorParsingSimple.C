// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//
#include <actuator/ActuatorBulk.h>
#include <actuator/ActuatorBulkSimple.h>
#include <NaluParsing.h>
#include <actuator/ActuatorParsingSimple.h>
#include <NaluEnv.h>

namespace sierra {
namespace nalu {

namespace {

std::vector<double> 
extend_double_vector(std::vector<double> vec, const int N)
{
  if ((vec.size() != 1) && (vec.size() != N))
    throw std::runtime_error("Vector is not of size 1 or "+std::to_string(N));
  if (vec.size() == 1) 
    { // Extend the vector to size N
      std::vector<double> newvec(N, vec[0]);
      return newvec;
    }
  if (vec.size() == N) 
    return vec;
  return vec;  // Should not get here
}

void
readTurbineData(int iTurb, ActuatorMetaSimple& actMetaSimple, YAML::Node turbNode)
{
  fast::fastInputs& fi = actMetaSimple.fastInputs_;
  // Read turbine data for a given turbine using the YAML node
  get_required(turbNode, "turb_id", fi.globTurbineData[iTurb].TurbID);
  get_required(
    turbNode, "fast_input_filename",
    fi.globTurbineData[iTurb].FASTInputFileName);
  get_required(
    turbNode, "restart_filename",
    fi.globTurbineData[iTurb].FASTRestartFileName);

  get_required(
    turbNode, "turbine_base_pos", fi.globTurbineData[iTurb].TurbineBasePos);
  if (turbNode["turbine_hub_pos"]) {
    NaluEnv::self().naluOutputP0()
      << "WARNING::turbine_hub_pos is not used. "
      << "The hub location is computed in OpenFAST and is controlled by the "
         "ElastoDyn input file.";
  }
  get_required(
    turbNode, "num_force_pts_blade",
    fi.globTurbineData[iTurb].numForcePtsBlade);

  actMetaSimple.maxNumPntsPerBlade_ = std::max(
    actMetaSimple.maxNumPntsPerBlade_,
    fi.globTurbineData[iTurb].numForcePtsBlade);

  get_required(
    turbNode, "num_force_pts_tower", fi.globTurbineData[iTurb].numForcePtsTwr);

  get_if_present_no_default(
    turbNode, "nacelle_cd", fi.globTurbineData[iTurb].nacelle_cd);
  get_if_present_no_default(
    turbNode, "nacelle_area", fi.globTurbineData[iTurb].nacelle_area);
  get_if_present_no_default(
    turbNode, "air_density", fi.globTurbineData[iTurb].air_density);

  int* numBlades = &(actMetaSimple.nBlades_(iTurb));
  *numBlades = 3;
  get_if_present_no_default(turbNode, "num_blades", *numBlades);
  ThrowErrorMsgIf(
    *numBlades != 3 && *numBlades != 2,
    "ERROR::ActuatorParsingSimple::Currently only 2 and 3 bladed turbines are "
    "supported.");

  if (actMetaSimple.is_disk()) {
    get_if_present_no_default(
      turbNode, "num_swept_pts", actMetaSimple.nPointsSwept_(iTurb));
    actMetaSimple.useUniformAziSampling_(iTurb) =
      actMetaSimple.nPointsSwept_(iTurb) != 0;
    ThrowErrorMsgIf(
      *numBlades != 3, "The ActuatorDisk model requires a base 3 bladed "
                       "turbine, but a 2 bladed turbine was supplied.");
  }

  actMetaSimple.numPointsTurbine_.h_view(iTurb) =
    1 // hub
    + fi.globTurbineData[iTurb].numForcePtsTwr +
    fi.globTurbineData[iTurb].numForcePtsBlade * (*numBlades);
  actMetaSimple.numPointsTotal_ += actMetaSimple.numPointsTurbine_.h_view(iTurb);
}
} // namespace

ActuatorMetaSimple
actuator_Simple_parse(const YAML::Node& y_node, const ActuatorMeta& actMeta)
{
  ActuatorMetaSimple actMetaSimple(actMeta);

  NaluEnv::self().naluOutputP0() << "In actuator_Simple_parse() "<< std::endl; //LCCOUT

  const YAML::Node y_actuator = y_node["actuator"];
  ThrowErrorMsgIf(
		  !y_actuator, "actuator argument is "
                 "missing from yaml node passed to actuator_Simple_parse");

  // Load the debug option
  const YAML::Node debug_output = y_actuator["debug_output"];
  if (debug_output) 
    actMetaSimple.debug_output_ = debug_output.as<bool>();
  else
    actMetaSimple.debug_output_ = false;

  size_t n_simpleblades_;
  get_required(y_actuator, "n_simpleblades", n_simpleblades_);
  actMetaSimple.n_simpleblades_ =  n_simpleblades_;
  NaluEnv::self().naluOutputP0() << "N blade: " << actMetaSimple.n_simpleblades_<< std::endl; //LCCOUT
  if (actMetaSimple.n_simpleblades_ > 0) {
    for (int iBlade= 0; iBlade < n_simpleblades_; iBlade++) {
      //NaluEnv::self().naluOutputP0() << "Reading blade: " << iBlade<< std::endl; //LCCOUT
      const YAML::Node cur_blade =
	y_actuator["Blade" + std::to_string(iBlade)];

      size_t num_force_pts_blade;
      get_required(cur_blade, "num_force_pts_blade", num_force_pts_blade);
      actMetaSimple.num_force_pts_blade_.h_view(iBlade) = num_force_pts_blade;
      actMetaSimple.numPointsTurbine_.h_view(iBlade)    = num_force_pts_blade;

      if (actMetaSimple.debug_output_)
	NaluEnv::self().naluOutputP0() 
	  << "Reading blade: " << iBlade
	  << " num_force_pts_blade: "
	  << actMetaSimple.numPointsTurbine_.h_view(iBlade) << std::endl; //LCCOUT

      // Get the epsilon
      const YAML::Node epsilon_chord = cur_blade["epsilon_chord"];
      const YAML::Node epsilon = cur_blade["epsilon"];
      if (epsilon && epsilon_chord) {
	throw std::runtime_error(
	  "epsilon and epsilon_chord have both been specified for Blade " +
	  std::to_string(iBlade) + "\nYou must pick one or the other.");
      }
      std::vector<double> epsilonTemp(3);
      
      // only require epsilon
      if (epsilon.Type() == YAML::NodeType::Scalar) {
	double isotropicEpsilon;
	get_required(cur_blade, "epsilon", isotropicEpsilon);
	actMetaSimple.isotropicGaussian_ = true;
	for (int j = 0; j < 3; j++) {
	  actMetaSimple.epsilon_.h_view(iBlade, j) = isotropicEpsilon;
	}
      } else {
	get_required(cur_blade, "epsilon", epsilonTemp);
	for (int j = 0; j < 3; j++) {
	  actMetaSimple.epsilon_.h_view(iBlade, j) = epsilonTemp[j];
	}
	if (
	    epsilonTemp[0] == epsilonTemp[1] &&
	    epsilonTemp[1] == epsilonTemp[2]) {
	  actMetaSimple.isotropicGaussian_ = true;
	} 
      }
      // check epsilon values
      for (int j = 0; j < 3; j++) {
	if (actMetaSimple.epsilon_.h_view(iBlade, j) <= 0.0) {
	  throw std::runtime_error(
              "ERROR:: zero value for epsilon detected. "
              "All epsilon components must be greater than zero");
          }
      }
      // Handle blade properties
      // p1 
      std::vector<double> p1Temp(3);
      get_required(cur_blade, "p1", p1Temp);
      for (int j = 0; j < 3; j++) {
	actMetaSimple.p1_.h_view(iBlade, j) = p1Temp[j];
      }
      // p2
      std::vector<double> p2Temp(3);
      get_required(cur_blade, "p2", p2Temp);
      for (int j = 0; j < 3; j++) {
	actMetaSimple.p2_.h_view(iBlade, j) = p2Temp[j];
      }
      // p1zeroAOA
      std::vector<double> p1zeroAOATemp(3);
      Coordinates p1zeroAOA;
      get_required(cur_blade, "p1_zero_alpha_dir", p1zeroAOATemp);
      double p1zeroAOAnorm = sqrt(p1zeroAOATemp[0]*p1zeroAOATemp[0] + 
				  p1zeroAOATemp[1]*p1zeroAOATemp[1] + 
				  p1zeroAOATemp[2]*p1zeroAOATemp[2]);
      for (int j = 0; j < 3; j++) {
	actMetaSimple.p1zeroalphadir_.h_view(iBlade, j) = 
	  p1zeroAOATemp[j]/p1zeroAOAnorm;
      }
      p1zeroAOA.x_ = actMetaSimple.p1zeroalphadir_.h_view(iBlade, 0);
      p1zeroAOA.y_ = actMetaSimple.p1zeroalphadir_.h_view(iBlade, 1);
      p1zeroAOA.z_ = actMetaSimple.p1zeroalphadir_.h_view(iBlade, 2);

      // Calculate some stuff
      // Span direction
      Coordinates spandir;
      spandir.x_ = p2Temp[0] - p1Temp[0];
      spandir.y_ = p2Temp[1] - p1Temp[1];
      spandir.z_ = p2Temp[2] - p1Temp[2];
      double spandirnorm = sqrt(spandir.x_*spandir.x_ + spandir.y_*spandir.y_ +
			 spandir.z_*spandir.z_);
      spandir.x_ = spandir.x_/spandirnorm;
      spandir.y_ = spandir.y_/spandirnorm;
      spandir.z_ = spandir.z_/spandirnorm;
      actMetaSimple.spandir_.h_view(iBlade, 0) = spandir.x_;
      actMetaSimple.spandir_.h_view(iBlade, 1) = spandir.y_;
      actMetaSimple.spandir_.h_view(iBlade, 2) = spandir.z_;
      // Chord normal direction
      Coordinates chordnormaldir;
      chordnormaldir.x_ = p1zeroAOA.y_*spandir.z_ - p1zeroAOA.z_*spandir.y_;
      chordnormaldir.y_ = p1zeroAOA.z_*spandir.x_ - p1zeroAOA.x_*spandir.z_;
      chordnormaldir.z_ = p1zeroAOA.x_*spandir.y_ - p1zeroAOA.y_*spandir.x_;
      actMetaSimple.chordnormaldir_.h_view(iBlade, 0) = chordnormaldir.x_;
      actMetaSimple.chordnormaldir_.h_view(iBlade, 1) = chordnormaldir.y_;
      actMetaSimple.chordnormaldir_.h_view(iBlade, 2) = chordnormaldir.z_;

      // output directions
      if (actMetaSimple.debug_output_) {
	NaluEnv::self().naluOutputP0()  // LCCOUT
	  << "Blade: " << iBlade << " p1zeroAOA dir: "
	  <<p1zeroAOA.x_<<" "<<p1zeroAOA.y_<<" "<<p1zeroAOA.z_<< std::endl;
	NaluEnv::self().naluOutputP0()  // LCCOUT
	  << "Blade: " << iBlade << " Span dir: "
	  <<spandir.x_<<" "<<spandir.y_<<" "<<spandir.z_<< std::endl; 
	NaluEnv::self().naluOutputP0() // LCCOUT
	  << "Blade: " << iBlade 
	  << " chord norm dir: "<<std::setprecision(5)
	  <<chordnormaldir.x_<<" "<<chordnormaldir.y_<<" "<<chordnormaldir.z_<< std::endl; 
      }
		
      // Chord definitions
      const YAML::Node chord_table = cur_blade["chord_table"];
      std::vector<double> chordtemp;
      if (chord_table) 
	chordtemp = chord_table.as<std::vector<double>>();
      else 
	throw std::runtime_error("ActuatorSimpleNGP: missing chord_table");
      actMetaSimple.chord_table_.push_back(extend_double_vector(chordtemp, num_force_pts_blade)); 

      // twist definitions
      const YAML::Node twist_table = cur_blade["twist_table"];
      std::vector<double> twisttemp;
      if (twist_table)
	twisttemp = twist_table.as<std::vector<double>>();
      else
	throw std::runtime_error("ActuatorSimpleNGP: missing twist_table");
      actMetaSimple.twist_table_.push_back(extend_double_vector(twisttemp, num_force_pts_blade)); 

      // Polar tables
      // --- aoa ---
      const YAML::Node aoa_table = cur_blade["aoa_table"];
      std::vector<double> aoatemp;      
      if (aoa_table)
	aoatemp = aoa_table.as<std::vector<double>>();
      else
	throw std::runtime_error("ActuatorSimpleNGP: missing aoa_table");
      actMetaSimple.aoa_polartable_.push_back(aoatemp); 
      size_t polartableN = aoatemp.size();

      // --- cl ---
      const YAML::Node cl_table = cur_blade["cl_table"];
      std::vector<double> cltemp;      
      if (cl_table)
	cltemp = cl_table.as<std::vector<double>>();
      else
	throw std::runtime_error("ActuatorSimpleNGP: missing cl_table");
      actMetaSimple.cl_polartable_.push_back(extend_double_vector(cltemp, polartableN)); 

      // --- cd ---
      const YAML::Node cd_table = cur_blade["cd_table"];
      std::vector<double> cdtemp;      
      if (cd_table)
	cdtemp = cd_table.as<std::vector<double>>();
      else
	throw std::runtime_error("ActuatorSimpleNGP: missing cd_table");
      actMetaSimple.cd_polartable_.push_back(extend_double_vector(cdtemp, polartableN)); 

	
    } // End loop over blades
  }else {
      throw std::runtime_error("Number of simple blades <= 0 ");
  }

  throw std::runtime_error("ActuatorSimple: loading blades");
  // ---vvv--- FAST STUFF TO DELETE ---vvvv---
  bool INCLUDEFASTSTUFF=false;
  if (INCLUDEFASTSTUFF) {
  fast::fastInputs& fi = actMetaSimple.fastInputs_;
  fi.comm = NaluEnv::self().parallel_comm();
  fi.nTurbinesGlob = actMetaSimple.numberOfActuators_;

  if (fi.nTurbinesGlob > 0) {
    fi.dryRun = false;
    get_if_present(y_actuator, "debug", fi.debug, false);
    get_required(y_actuator, "t_start", fi.tStart);
    std::string simStartType = "na";
    get_required(y_actuator, "simStart", simStartType);
    if (simStartType == "init") {
      if (fi.tStart == 0) {
        fi.simStart = fast::init;
      } else {
        throw std::runtime_error(
          "actuators: simStart type not consistent with start time for FAST");
      }
    } else if (simStartType == "trueRestart") {
      fi.simStart = fast::trueRestart;
    } else if (simStartType == "restartDriverInitFAST") {
      fi.simStart = fast::restartDriverInitFAST;
    }
    get_required(y_actuator, "n_every_checkpoint", fi.nEveryCheckPoint);
    get_required(y_actuator, "dt_fast", fi.dtFAST);

    get_required(y_actuator, "t_max", fi.tMax);

    if (y_actuator["super_controller"]) {
      get_required(y_actuator, "super_controller", fi.scStatus);
      get_required(y_actuator, "sc_libFile", fi.scLibFile);
      get_required(y_actuator, "num_sc_inputs", fi.numScInputs);
      get_required(y_actuator, "num_sc_outputs", fi.numScOutputs);
    }

    fi.globTurbineData.resize(fi.nTurbinesGlob);

    for (int iTurb = 0; iTurb < fi.nTurbinesGlob; iTurb++) {
      if (y_actuator["Turbine" + std::to_string(iTurb)]) {

        const YAML::Node cur_turbine =
          y_actuator["Turbine" + std::to_string(iTurb)];

        get_required(
          cur_turbine, "turbine_name", actMetaSimple.turbineNames_[iTurb]);

        std::string turbFileName;
        get_if_present(
          cur_turbine, "file_to_dump_turb_pts",
          actMetaSimple.turbineOutputFileNames_[iTurb]);

        get_if_present_no_default(
          cur_turbine, "fllt_correction",
          actMetaSimple.filterLiftLineCorrection_);

        ThrowErrorMsgIf(
          actMetaSimple.filterLiftLineCorrection_,
          "Filtered lifting line correction has not been implemented in the NGP"
          " actuator models yet.  Please use ActLineFAST instead.");
        // The value epsilon / chord [non-dimensional]
        // This is a vector containing the values for:
        //   - chord aligned (x),
        //   - tangential to chord (y),
        //   - spanwise (z)
        const YAML::Node epsilon_chord = cur_turbine["epsilon_chord"];
        const YAML::Node epsilon = cur_turbine["epsilon"];
        if (epsilon && epsilon_chord) {
          throw std::runtime_error(
            "epsilon and epsilon_chord have both been specified for Turbine " +
            std::to_string(iTurb) + "\nYou must pick one or the other.");
        }
        if (epsilon && actMetaSimple.filterLiftLineCorrection_) {
          throw std::runtime_error(
            "epsilon and fllt_correction have both been specified for "
            "Turbine " +
            std::to_string(iTurb) +
            "\nepsilon_chord and epsilon_min should be used with "
            "fllt_correction.");
        }

        std::vector<double> epsilonTemp(3);
        if (
          actMeta.actuatorType_ == ActuatorType::ActLineFASTNGP ||
          actMeta.actuatorType_ == ActuatorType::ActDiskFASTNGP) {
          // only require epsilon
          if (epsilon.Type() == YAML::NodeType::Scalar) {
            double isotropicEpsilon;
            get_required(cur_turbine, "epsilon", isotropicEpsilon);
            actMetaSimple.isotropicGaussian_ = true;
            for (int j = 0; j < 3; j++) {
              actMetaSimple.epsilon_.h_view(iTurb, j) = isotropicEpsilon;
            }
          } else {
            get_required(cur_turbine, "epsilon", epsilonTemp);
            for (int j = 0; j < 3; j++) {
              actMetaSimple.epsilon_.h_view(iTurb, j) = epsilonTemp[j];
            }
            if (
              epsilonTemp[0] == epsilonTemp[1] &&
              epsilonTemp[1] == epsilonTemp[2]) {
              actMetaSimple.isotropicGaussian_ = true;
            } else if (actMeta.actuatorType_ == ActuatorType::ActDiskFASTNGP) {
              throw std::runtime_error("ActDiskFASTNGP does not currently "
                                       "support anisotropic epsilons.");
            }
          }
          // single value epsilon
          // multi value epsilon
        } else if (actMeta.actuatorType_ == ActuatorType::AdvActLineFASTNGP) {
          // require epsilon chord and epsilon min
          get_required(cur_turbine, "epsilon_chord", epsilonTemp);
          for (int j = 0; j < 3; j++) {
            if (epsilonTemp[j] <= 0.0) {
              throw std::runtime_error(
                "ERROR:: zero value for epsilon_chord detected. "
                "All epsilon components must be greater than zero");
            }
            actMetaSimple.epsilonChord_.h_view(iTurb, j) = epsilonTemp[j];
          }

          // Minimum epsilon allowed in simulation. This is required when
          //   specifying epsilon/chord
          get_required(cur_turbine, "epsilon_min", epsilonTemp);
          for (int j = 0; j < 3; j++) {
            actMetaSimple.epsilon_.h_view(iTurb, j) = epsilonTemp[j];
          }
        }
        // check epsilon values
        for (int j = 0; j < 3; j++) {
          if (actMetaSimple.epsilon_.h_view(iTurb, j) <= 0.0) {
            throw std::runtime_error(
              "ERROR:: zero value for epsilon detected. "
              "All epsilon components must be greater than zero");
          }
        }

        // An epsilon value used for the tower
        const YAML::Node epsilon_tower = cur_turbine["epsilon_tower"];
        // If epsilon tower is given store it.
        // If not, use the standard epsilon value
        if (epsilon_tower) {
          epsilonTemp = epsilon_tower.as<std::vector<double>>();
          for (int j = 0; j < 3; j++) {
            actMetaSimple.epsilonTower_.h_view(iTurb, j) = epsilonTemp[j];
          }
        } else {
          for (int j = 0; j < 3; j++) {
            actMetaSimple.epsilonTower_.h_view(iTurb, j) =
              actMetaSimple.epsilon_.h_view(iTurb, j);
          }
        }

        readTurbineData(iTurb, actMetaSimple, cur_turbine);
      } else {
        throw std::runtime_error(
          "Node for Turbine" + std::to_string(iTurb) +
          " not present in input file or I cannot read it");
      }
    }

  } else {
    throw std::runtime_error("Number of turbines <= 0 ");
  }
  } // INCLUDEFASTSTUFF
  return actMetaSimple;
}

} // namespace nalu
} // namespace sierra
