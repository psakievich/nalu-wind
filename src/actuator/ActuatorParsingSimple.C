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
    actMetaSimple.numPointsTotal_ = 0;
    for (int iBlade= 0; iBlade < n_simpleblades_; iBlade++) {

      const YAML::Node cur_blade =
	y_actuator["Blade" + std::to_string(iBlade)];

      size_t num_force_pts_blade;
      get_required(cur_blade, "num_force_pts_blade", num_force_pts_blade);
      actMetaSimple.num_force_pts_blade_.h_view(iBlade) = num_force_pts_blade;
      actMetaSimple.numPointsTurbine_.h_view(iBlade)    = num_force_pts_blade;
      actMetaSimple.numPointsTotal_                    += num_force_pts_blade;

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
      std::vector<double> chord_table_extended = 
	extend_double_vector(chordtemp, num_force_pts_blade);
      actMetaSimple.chord_table_.push_back(chord_table_extended); 

      // twist definitions
      const YAML::Node twist_table = cur_blade["twist_table"];
      std::vector<double> twisttemp;
      if (twist_table)
	twisttemp = twist_table.as<std::vector<double>>();
      else
	throw std::runtime_error("ActuatorSimpleNGP: missing twist_table");
      actMetaSimple.twist_table_.push_back(extend_double_vector(twisttemp, num_force_pts_blade)); 

      // Calculate elem areas
      std::vector<double> elemareatemp(num_force_pts_blade, 0.0);
      std::vector<double> dx(3,0.0);
      for (int i=0; i<3; i++) 
	dx[i] = (p2Temp[i] - p1Temp[i])/(double)num_force_pts_blade;
      double dx_norm = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      for (int i=0; i<num_force_pts_blade; i++)
	elemareatemp[i] = dx_norm*chord_table_extended[i];
      actMetaSimple.elem_area_.push_back(elemareatemp);

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

  NaluEnv::self().naluOutputP0() << " actMetaSimple.numPointsTotal_ = "
		      << actMetaSimple.numPointsTotal_<< std::endl; // LCCOUT
  NaluEnv::self().naluOutputP0() << "Done actuator_Simple_parse()"<< std::endl; // LCCOUT
  return actMetaSimple;
}

} // namespace nalu
} // namespace sierra
