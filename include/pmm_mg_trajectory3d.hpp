/**
 * @file pmm_mg_trajectory3d.hpp
 * @author Krystof Teissing (teisskry@gmail.com), Matej Novosad (novosma2@fel.cvut.cz)
 * @version 0.1
 * @date 2024-09-20
 * 
 */

#pragma once

#include "pmm_trajectory3d.hpp"
#include "common.hpp"
#include <iostream>
#include <fstream>


namespace pmm{

// Enum to store the current state of the two-segment scope velocity optimization case
enum seg_optim_case {
    MIN_MIN,
    MIN_SYNC,
    SYNC_MIN,
    SYNC_SYNC,
};

class PMM_MG_Trajectory3D{
    public:
        /**
         * @brief Default constructor
         * 
         */
        PMM_MG_Trajectory3D();

        /**
         * @brief Construct a new 3D PMM multi-goal trajectory object
         * 
         * @param n_path_segments Number of segments with defined start and end configurations
         */
        PMM_MG_Trajectory3D(const int n_path_segments);

        /**
         * @brief Construct a new 3D PMM multi-goal trajectory object from provided initial conditions
         *        and compute the trajectory for all segments
         * 
         * @param p Vector of 3D goals
         * @param v Vector of initial 3D velocities; v[0] is the given start velocity and v[-1] is the 
         *          given end velocity, both are handled as fixed values. The remaining ones can be further
         *          optimized.
         * @param a_max Vector of symmetric acceleration limits for each axis
         * @param v_max Vector of symmetric velocity limits for each axis
         */
        PMM_MG_Trajectory3D(std::vector<Vector<3>> p, std::vector<Vector<3>> v, Vector<3> a_max, Vector<3> v_max);

        /**
         * @brief Construct a new 3D PMM multi-goal trajectory object from provided initial conditions
         *        and compute the trajectory for all segments
         * 
         * @param p Vector of 3D goals
         * @param v Vector of initial 3D velocities; v[0] is the given start velocity and v[-1] is the 
         *          given end velocity, both are handled as fixed values. The remaining ones can be further
         *          optimized.
         * @param a_max Vector of max acceleration limits for each axis
         * @param a_min Vector of min acceleration limits for each axis
         * @param v_max Vector of symmetric velocity limits for each axis
         */
        PMM_MG_Trajectory3D(std::vector<Vector<3>> p, std::vector<Vector<3>> v, Vector<3> a_min, Vector<3> a_max, Vector<3> v_max);

        /**
         * @brief Construct a new 3D PMM multi-goal trajectory object from provided initial conditions
         *        and compute the trajectory for all segments
         * 
         * @param p Vector of 3D goals
         * @param v Vector of initial 3D velocities; v[0] is the given start velocity and v[-1] is the 
         *          given end velocity, both are handled as fixed values. The remaining ones can be further
         *          optimized.
         * @param a_max Vector of max acceleration limits for each axis
         * @param a_min Vector of min acceleration limits for each axis
         * @param v_max Vector of symmetric velocity limits for each axis
         * @param equalize_axes Flag stating if axes should be equalized
         */
        PMM_MG_Trajectory3D(std::vector<Vector<3>> p, std::vector<Vector<3>> v, Vector<3> a_min, Vector<3> a_max, Vector<3> v_max, bool equalize_axes);

        /**
         * @brief Construct a new 3D PMM multi-goal trajectory object from provided initial conditions
         *        and compute the trajectory for all segments assuming limited acceleration norm
         * 
         * @param p Vector of 3D goals
         * @param v Vector of initial 3D velocities; v[0] is the given start velocity and v[-1] is the 
         *          given end velocity, both are handled as fixed values. The remaining ones can be further
         *          optimized.
         * @param max_acc_norm Maximal acceleration norm
         * @param max_acc_norm Maximal velocity norm
         * @param td_max_iter Maximum number of iterations for Thrust Decomposition algorithm
         * @param td_acc_norm_precision_limit Stopping criterion for TD
         * @param drag Flag determining if drag modeling is used
         * @param debug Flag determining if outputs should be shown in the console
         */
        PMM_MG_Trajectory3D(std::vector<Vector<3>> p, std::vector<Vector<3>> v, const Scalar max_acc_norm, const Scalar max_vel_norm, const int td_max_iter=10,
                            const Scalar td_acc_norm_precision_limit=1e-3, const bool drag=false, const bool debug=false);
        
        /**
         * @brief Initialize, compute and optimize trajectory given the provided initial conditions and optimization parameters
         * 
         * @param waypoints Waypoints which shall be connected with the trajectory
         * @param start_vel Start velocity - velocity vector in the first waypoint
         * @param end_vel End velocity - velocity in the last waypoint
         * @param max_acc_norm Norm of the maximal acceleration of the drone given its limited thrust
         * @param max_vel_norm Maximal velocity norm
         * @param dT_precision Stopping criterion, when the duration change of the updated trajectory is below this threshold, the optimization is terminated.
         * @param GD_max_iter Maximum number of iterations of the gradient descent optimization.
         * @param GD_alpha Gradient descent step.
         * @param GD_alpha_red_fctr Factor from <0,1> by which the gradient descent step should be reduced when close to singularity.
         * @param GD_alpha_min_thr Minimal value for gradient descent step, when the step size is reduced under this value, the corresponding axis is not updated any more.
         * @param TD_max_iter Maximum number of iterations for Thrust Decomposition algorithm
         * @param TD_acc_precision Stopping criterion for Thrust Decomposition algorithm, precision of the resulting maximal acceleration norm used in the trajectory
         * @param scnd_round_optim Flag, if second optimization round should be run, where the trajectory is optimized while the Thrust Decomposition algorithm is used
         * @param GD_max_iter_2nd_rnd Maximum number of iterations of the gradient descent optimization in the second round.
         * @param GD_alpha_2nd_rnd Gradient descent step in the second round.
         * @param GD_alpha_red_fctr_2nd_rnd Factor from <0,1> by which the gradient descent step should be reduced when close to singularity in the second optimization round.
         * @param GD_alpha_min_thr_2nd_rnd Minimal value for gradient descent step, when the step size is reduced under this value, the corresponding axis is not updated any more.
         * @param drag Flag determining if drag modeling is used
         * @param debug 
         */
        PMM_MG_Trajectory3D(std::vector<Vector<3>> waypoints, const Vector<3> start_vel, const Vector<3> end_vel,
                                            const Scalar max_acc_norm, const Scalar max_vel_norm, const Scalar dT_precision = 1e-3, const int GD_max_iter = 30, 
                                            const Scalar GD_alpha = 30, const Scalar GD_alpha_red_fctr = 0.4, const Scalar GD_alpha_min_thr = 1e-2, 
                                            const Scalar TD_max_iter = 10, const Scalar TD_acc_precision = 1e-3, const bool scnd_round_optim = true,
                                            const int GD_max_iter_2nd_rnd = 10, const Scalar GD_alpha_2nd_rnd = 30, const Scalar GD_alpha_red_fctr_2nd_rnd = 0.4,
                                            const Scalar GD_alpha_min_thr_2nd_rnd = 1e-2, const bool drag=false, const bool debug = false);

        /**
         * @brief Method for optimizing velocities in via trajectory points using gradient descent to obtain a trajectory with smaller duration.
         *        Per-axis acceleration limits or limited acceleration norm are considered depending on the initially constructed PMM trajectory;
         * 
         * @param gd_step Gradient descent step
         * @param gd_step_reduction_factor Factor from <0,1> by which the gradient descent step should be reduced when close to singularity.
         * @param gd_step_min_threshold Minimal value for gradient descent step, when the step size is reduced under this value, the corresponding axis is not updated any more.
         * @param max_iter Maximum number of iterations of the gradient descent.
         * @param dT_threshold Stopping criterion, when the duration change of the updated trajectory is below this threshold, the optimization is terminated.
         * @param drag Flag determining if drag modeling is used
         * @param single_axis_optim_idx Indicates which trajectory axis should be optimized: (-1) -> all axis, (0) -> x-axis, (1) -> y-axis, (2) -> z-axis.
         * @param debug Flag for printing debugging info in terminal
         * @param export_trajectory_to_file Flag for exporting the trajectory data in each iteration of the optimization in to a csv file file_name
         * @param file_name Location of the file the optimization data should be saved to.
         * @param export_sampled_trajectories Flag that determines is if sampled trajectories are exported for every iteration of the optimization process
         * @return int 
         */
        int optimize_velocities_at_positions(Scalar gd_step=0.5, Scalar gd_step_reduction_factor=0.5, Scalar gd_step_min_threshold=1e-3,
                                            int max_iter=100, Scalar dT_threshold=1e-5, const bool drag=false, int single_axis_optim_idx=-1, bool debug=false,
                                            bool export_trajectory_to_file=false, std::string file_name="out.csv", bool export_sampled_trajectories=false);
        
        /**
         * @brief Check whether current trajectory is feasible.
         * 
         * @param drag Flag
         * @param verbose Flag determining if additional information regarding the thrust acceleration usage should be written out to the console
         * @return true 
         * @return false 
         */
        bool checkTrajectory(bool drag, bool verbose=false);

        /**
         * @brief Method for extracting computed trajectory values from PMMTrajectory
         *  class and storing them in one multi-point trajectory. It is assumed that the end
         *  configuration of a trajectory segment is the same as the start configuration of the
         *  following trajectory segment.
         * 
         * @param tr pointer to PMMTrajectory trajectory
         * @param segment_idx Index of trajectory segment starting with 0.
         * @param override Flag that disables error message when already filled trajectory segment is filled again.
         * @return int 
         */
        int fillTrajectorySegment(PointMassTrajectory3D *tr, int segment_idx, const bool override=false);

        /**
         * @brief Recompute trajectory for all segments given new values of velocities at each time stamp
         * 
         * @param v_new New velocity values for each time stamp and axis
         * @param a_max Symmetric maximal acceleration bound
         * @param equalize_times select if the times in all three axis should be equalized
         * @return int Success
         */
        int updateTrajectoryGivenVelocities(std::vector<Vector<3>> v_new, const bool equalize_times, const bool sync_traj=false, const bool drag=false);

        /**
         * @brief Recompute trajectory given new values of starting velocities for a given trajectory segment
         * 
         * @param v_new New velocity values of starting velocities of the given segment
         * @param segment_idx index of the given segment
         * @param sync_traj Flag determining if sync segment should be computed
         * @return int Success
         */
        int updateTrajectoryGivenVelocities(Vector<3> v_new, const int segment_idx, const bool equalize_times, const bool sync_traj=false, const bool drag=false);

        /**
         * @brief Compute gradient dt/dv12 given two consecutive MIN and SYNC trajectory segments
         * 
         * @param ax_idx index of the trajectory axis
         * @param sync_idx SYNC trajectory segment index
         * @return tuple (Gradient dt/dv12, grad lower bound, grad upper bound)
         */
        std::tuple<Scalar, Scalar, Scalar> get_MIN_SYNC_grad_and_bounds(const int ax_idx, const int min_idx);

        /**
         * @brief Compute gradient dt/dv12 given two consecutive SYNC and MIN trajectory segments
         * 
         * @param ax_idx index of the trajectory axis
         * @param min_idx MIN trajectory segment index
         * @return tuple (Gradient dt/dv12, grad lower bound, grad upper bound)
         */
        std::tuple<Scalar, Scalar, Scalar> get_SYNC_MIN_grad_and_bounds(const int ax_idx, const int min_idx);

        /**
         * @brief Compute gradient dt/dv12 given two consecutive MIN and MIN trajectory segments
         * 
         * @param ax_idx index of the trajectory axis
         * @param second_min_idx second MIN trajectory segment index
         * @return tuple (Gradient dt/dv12, grad lower bound, grad upper bound)
         */
        std::tuple<Scalar, Scalar, Scalar> get_MIN_MIN_grad_and_bounds(const int ax_idx, const int second_min_idx);
        
        /**
         * @brief Function that returns velocity bounds on the boundary velocity of the single-axis SYNC-MIN trajectory.
         *        The bounds are given by the limited acceleration scale value of the SYNC trajectory
         * 
         * @param seg_idx index of the SYNC segment
         * @param ax_idx trajectory axis index
         * @return std::tuple<Scalar, Scalar>  min and max velocity bounds, respectively.
         */
        std::tuple<Scalar, Scalar> get_SYNC_MIN_vel_bounds(const int seg_idx, const int ax_idx);
        
        /**
         * @brief Function that returns velocity bounds on the boundary velocity of the single-axis MIN-SYNC trajectory.
         *        The bounds are given by the limited acceleration scale value of the SYNC trajectory
         * 
         * @param seg_idx index of the SYNC segment
         * @param ax_idx trajectory axis index
         * @return std::tuple<Scalar, Scalar>  min and max velocity bounds, respectively.
         */
        std::tuple<Scalar, Scalar> get_MIN_SYNC_vel_bounds(const int seg_idx, const int ax_idx);
        
        /**
         * @brief Get a sampled trajectory as an tuple of time, position, velocity and acceleration sampled with sampling period sampling_period. 
         * 
         * @param sampling_period Sampling period with which the trajectory should be sampled
         * @return std::tuple<std::vector<Scalar>, std::vector<Vector<3>>, std::vector<Vector<3>>, std::vector<Vector<3>>> Tuple with vectors of sampled times,
         *         positions, velocities and acceleration.
         */
        std::tuple<std::vector<Scalar>, std::vector<Vector<3>>, std::vector<Vector<3>>, std::vector<Vector<3>>> get_sampled_trajectory(Scalar sampling_period);

        /**
         * @brief Get a sampled single-axis trajectory as an tuple of time, position, velocity and acceleration sampled with sampling period sampling_period.
         * 
         * @param sampling_period Sampling period with which the trajectory should be sampled
         * @param ax_idx Axis index
         * @return std::tuple<std::vector<Scalar>, std::vector<Scalar>, std::vector<Scalar>, std::vector<Scalar>> 
         */
        std::tuple<std::vector<Scalar>, std::vector<Scalar>, std::vector<Scalar>, std::vector<Scalar>> get_sampled_single_axis_trajectory(Scalar sampling_period, int ax_idx);

        /**
         * @brief Get the sampled trajectory positions which are sampled according to the sampling_distance, where sampling_distance = max(dx, dy, dz), 
         *        WARNING: Valid only for BANG-BANG trajectory!!!
         * 
         * @param sampling_distance per-axis maximal distance between two waypoints
         * @return std::vector<Vector<3>> Trajectory path waypoints
         */
        std::tuple<std::vector<Scalar>, std::vector<Vector<3>>> get_sampled_trajectory_path(const Scalar sampling_distance); 

        /**
         * @brief Sample current trajectory and export it to a CSV file
         * 
         * @param sampling_period  Sampling period
         * @param output_file Output file location
         */
        void sample_and_export_trajectory(Scalar sampling_period, std::string output_file);

        /**
         * @brief Sample current trajectory and export it to a CSV file row-wise
         * 
         * @param sampling_period  Sampling period
         * @param output_file Output file location
         */
        void sample_and_export_trajectory_column(Scalar sampling_period, std::string output_file);

        /**
         * @brief Sample current unsynchronized trajectory and export it to a CSV file
         * 
         * @param sampling_period  Sampling period
         * @param output_file Output file location
         */
        void sample_and_export_unsynchronized_trajectory(Scalar sampling_period, std::string output_file);

        /**
         * @brief Sample current trajectory and export it to a CSV file
         * 
         * @param sampling_distance  per-axis maximal distance between two waypoints
         * @param output_file Output file location
         */
        void sample_and_export_trajectory_path(Scalar sampling_distance, std::string output_file);

        /**
         * @brief Reset all variables for new trajectory values
         * 
         */
        void reset();

        /**
         * @brief Export all trajectory data to a CSV file
         * 
         * @param file_name Path of the file where the values should be stored
         * @param optim_data Flag for exporting also optimization data
         */
        void exportTrajectoryToCSV(std::string file_name, bool optim_data = false);

        /**
         * @brief Write trajectory data to an opened CSV file
         * 
         * @param outStream Pointer to the ofstream object
         * @param optim_data Flag for exporting also optimization data
         */
        void exportTrajectoryToCSV(std::ofstream *outStream, bool optim_data = false);

        /**
         * @brief Write out elements of a vector in to a given file
         * 
         * @param v Vector of Scalars
         * @param outFile Pointer to the ofstream object
         */
        void writeVectorToFile(std::vector<Scalar> v, std::ofstream *outFile);

        /**
         * @brief Write out elements of a vector in to a given file
         * 
         * @param v 3D vector
         * @param outFile Pointer to the ofstream object
         */
        void writeVectorToFile(Vector<3> v, std::ofstream *outFile);

        /**
         * @brief Write out elements of a vector in to a given file
         * 
         * @param v 4D vector
         * @param outFile Pointer to the ofstream object
         */
        void writeVectorToFile(Vector<4> v, std::ofstream *outFile);

        /**
         * @brief Write out elements of a 3D vector in to a given file, separate line 
         *        for each axis
         * 
         * @param v std vector of 3D vectors
         * @param outFile Pointer to the ofstream object
         */
        void writeVectorToFile(std::vector<Vector<3>> v, std::ofstream *outFile);

        /**
         * @brief Write out elements of a 4D vector in to a given file, separate line 
         *        for each axis
         * 
         * @param v std vector of 3D vectors
         * @param outFile Pointer to the ofstream object
         */
        void writeVectorToFile(std::vector<Vector<4>> v, std::ofstream *outFile);

        /**
         * @brief Get the total duration of the trajectory assuming equalized axis times
         * 
         * @return Scalar 
         */
        Scalar duration();

        /**
         * @brief Get the duration of the trajectory segment for a given axis
         * 
         * @param segment_idx
         * @param axis_idx
         * @return Scalar 
         */
        Scalar segment_duration(int segment_idx, int axis_idx);

        /**
         * @brief Get the current values of start and end velocities of all trajectory segments
         * 
         * @return std::vector<Vector<3>> Vector of all start and end velocities
         */
        std::vector<Vector<3>> get_current_velocities();

        /**
         * @brief Get the current value of start and end positions of all trajectory segments
         * 
         * @return std::vector<Vector<3>> Vector of all start and end positions
         */
        std::vector<Vector<3>> get_current_positions();

        /**
         * @brief Copy trajectory values from the given trajectory
         * 
         * @param tr Trajectory from which values are copied
         */
        void copy_trajectory(const PMM_MG_Trajectory3D &tr);
        
        /**
         * @brief Operator definition
         * 
         * @param tr
         * @return PMM_MG_Trajectory3D 
         */
        PMM_MG_Trajectory3D operator=(const PMM_MG_Trajectory3D &tr);

        bool _exists = false;
        bool _thrust_limit_flag = false; // Flag that determines if the thrust acceleration is limited
        int n_segments;
        int n_segments_filled;
        Vector<3> _a_max;   // Vector of per-axis initial max acceleration
        Vector<3> _a_min;   // Vector of per-axis initial min acceleration
        Vector<3> _v_max;   // Vector of per-axis initial max velocity

        std::vector<Vector<3>> _p;
        std::vector<Vector<3>> _v;
        std::vector<Vector<3>> _a;
        std::vector<Vector<3>> _dt;
        std::vector<Vector<3>> _dtdv;

        std::vector<Eigen::Matrix<int, 3, 1>> _traj_type_vec; // 0: MIN, 1: SYNC, 2:ZERO
        std::vector<Eigen::Matrix<int, 3, 1>> _sol_idx; // Index of the solution (scale for SYNC)
        std::vector<Eigen::Matrix<bool, 3, 1>> _update_vel_flag;    // flag determining if the corresponding axis in the corresponding segment should be updated
        std::vector<Vector<3>> _ax_min_traj_duration;   // min per-axis trajectory duration before equalization
        std::vector<Vector<3>> _used_grad;  // gradients used for velocity optimization
        
        Scalar _max_acc_norm; // maximal acceleration norm
        Scalar _max_vel_norm; // maximal velocity norm
        int _TD_max_iter = 100; // Newtons method max iteration
        Scalar _TD_alpha_factor = 1.0;  // NM step scaling factor
        Scalar _TD_norm_precision = 1e-3;   // Stopping criterion tolerance

    private:
};

/**
 * @brief Decompose limited thrust-acceleration vector in to per-axis (xyz) acceleration limits while taking the gravitation into account;
 * 
 * @param max_acc_norm  Limited thrust-acceleration vector norm 
 * @return std::tuple<Vector<3>, Vector<3>> Min per-axis acceleration, Max per-axis acceleration
 */
std::tuple<Vector<3>, Vector<3>> compute_per_axis_acc_limits_from_limited_thrust(const Scalar max_acc_norm);

/**
 * @brief Compute initial velocity values for trajectory generation
 * 
 * @param waypoints Waypoints which should be connected by the trajectory
 * @param start_velocity Initial velocity vector
 * @param end_velocity End velocity vector
 * @param max_per_axis_acc Maximum per-axis acceleration (absolute value). Used for velocity norm estimation.
 * @param max_per_axis_vel Maximum per-axis velocity (absolute value).
 * @param const_vel_norm OPTIONAL, if non-negative, the const_vel_norm value is used for velocity norm in all waypoints, where
 *                       the angle between two consecutive headings is below (3/4)*PI, instead of velocity norm estimation using the max_per_axis_acc value
 * @return std::vector<Vector<3>> Vector of initial velocities in all waypoints
 */
std::vector<Vector<3>> compute_initial_velocities(std::vector<Vector<3>> waypoints, Vector<3> start_velocity,
                     Vector <3> end_velocity, const Scalar max_per_axis_acc, const Scalar max_per_axis_vel, const double const_vel_norm = -1);

} // namespace pmm