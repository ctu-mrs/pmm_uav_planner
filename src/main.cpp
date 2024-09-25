/**
 * @file main.cpp
 * @author Krystof Teissing (k.teissing@gmail.com), Matej Novosad (novosma2@fel.cvut.cz)
 * @version 0.1
 * @date 2024-09-20
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <csignal>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <Eigen/Eigen>
#include <chrono>

#include "pmm_trajectory3d.hpp"
#include "yaml-cpp/yaml.h"
#include "pmm_mg_trajectory3d.hpp"
#include "common.hpp"

using namespace pmm;

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

void signal_callback(int sig) {
  std::cout << "Signal " << sig << " received" << std::endl;
  exit(sig);
}

int plan_pmm_trajectory(std::string planner_config_file, std::string waypoints_config_file) {

  std::cout << "--------------------------- LADING PATH CONFIG ---------------------------" << std::endl;

  YAML::Node planner_config = YAML::LoadFile(planner_config_file);
  YAML::Node waypoints_config = YAML::LoadFile(waypoints_config_file);

  // Drone parameters ----------------------------------------------------------------------------/
  Scalar max_acc_norm, max_vel_norm, g;
  Vector<3> max_velocity, D;
  planner_config["uav"]["maximum_acceleration_norm"] >> max_acc_norm;
  planner_config["uav"]["maximum_velocity"] >> max_vel_norm;
  planner_config["uav"]["gravitational_acceleration"] >> g;
  planner_config["uav"]["drag_coefficients"] >> D;

  G = g;
  GVEC = {0, 0, -g};
  D_coeffs = -(Eigen::Matrix3d()<< D(0), 0, 0, 0, D(1), 0, 0, 0, D(2)).finished();

  // Velocity Optimization parameters ------------------------------------------------------------/
  // Thrust decomposition
  int TD_max_iter;
  Scalar TD_acc_precision;
  bool drag;
  planner_config["thrust_decomposition"]["drag"] >> drag;
  planner_config["thrust_decomposition"]["precision"] >> TD_acc_precision;
  planner_config["thrust_decomposition"]["max_iter"] >> TD_max_iter;

  // Gradient Method optimization parameters
  // First round
  Scalar alpha;
  Scalar alpha_reduction_factor;
  Scalar alpha_min_threshold;
  int max_iter;
  planner_config["velocity_optimization"]["first_run"]["alpha"] >> alpha;
  planner_config["velocity_optimization"]["first_run"]["alpha_reduction_factor"] >> alpha_reduction_factor;
  planner_config["velocity_optimization"]["first_run"]["alpha_min_threshold"] >> alpha_min_threshold;
  planner_config["velocity_optimization"]["first_run"]["max_iter"] >> max_iter;

  // Second round
  Scalar alpha2;
  Scalar alpha_reduction_factor2;
  Scalar alpha_min_threshold2;
  int max_iter2;
  planner_config["velocity_optimization"]["second_run"]["alpha"] >> alpha2;
  planner_config["velocity_optimization"]["second_run"]["alpha_reduction_factor"] >> alpha_reduction_factor2;
  planner_config["velocity_optimization"]["second_run"]["alpha_min_threshold"] >> alpha_min_threshold2;
  planner_config["velocity_optimization"]["second_run"]["max_iter"] >> max_iter2;

  double dT_precision;
  planner_config["velocity_optimization"]["dT_precision"] >> dT_precision;
  bool second_round_opt = true;

  bool debug, export_trajectory;
  Scalar sampling_step;
  std::string sampled_trajectory_file;
  planner_config["debug"] >> debug;
  planner_config["export"]["sampled_trajectory"] >> export_trajectory;
  planner_config["export"]["sampling_step"] >> sampling_step;
  planner_config["export"]["sampled_trajectory_file"] >> sampled_trajectory_file;

  // Load waypoint data --------------------------------------------------------------------------/
  Vector<3> start_velocity;
  Vector<3> end_velocity;
  Vector<3> start_position;
  Vector<3> end_position;
  waypoints_config["start"]["velocity"] >> start_velocity;
  waypoints_config["end"]["velocity"] >> end_velocity;
  waypoints_config["start"]["position"] >> start_position;
  waypoints_config["end"]["position"] >> end_position;
  std::vector<Vector<3>> path_waypoints;

  if (!parseArrayParam<Vector<3>>(waypoints_config, "waypoints", path_waypoints))
    std::cerr << "can not load param gates" << std::endl;

  path_waypoints.insert(path_waypoints.begin(), start_position);
  path_waypoints.push_back(end_position);

  std::cout << "Number of waypoints: " << path_waypoints.size() << std::endl;

  std::cout << "-------------------------- COMPUTING TRAJECTORY --------------------------" << std::endl;

  // Compute trajectory, see definition for parameter description
  PMM_MG_Trajectory3D mp_tr(path_waypoints, start_velocity, end_velocity, max_acc_norm, max_vel_norm, dT_precision, max_iter, alpha,
                            alpha_reduction_factor, alpha_min_threshold, TD_max_iter, TD_acc_precision, second_round_opt, 
                            max_iter2, alpha2, alpha_reduction_factor2, alpha_min_threshold2, drag, debug);

  std::cout << "Optimized trajectory duration after second optimization run with LTD: " << mp_tr.duration()  << " s." << std::endl;

  if (export_trajectory) {
    mp_tr.sample_and_export_trajectory(sampling_step, "scripts/trajectory_data/" + sampled_trajectory_file);
  }

  std::cout << "Visualize results by running scripts/plot_results.py!"<< std::endl;
  std::cout << "#----------------------------------------- END ------------------------------------------#" << std::endl << std::endl;
  return 0;

}

int main(int argc, char** argv) {

  if (argc < 3) {
    std::cerr << "ERROR: Please specify planner config file and waypoints config file!" << std::endl;
    return -1;
  }

  std::string planner_config_file = argv[1];
  std::string waypoints_config_file = argv[2];
  
  plan_pmm_trajectory(planner_config_file, waypoints_config_file);

  return 0;
}