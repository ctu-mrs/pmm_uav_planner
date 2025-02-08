/**
 * @file pmm_mg_trajectory3d.cpp
 * @author Krystof Teissing (teisskry@gmail.com), Matej Novosad (novosma2@fel.cvut.cz)
 * @version 0.1
 * @date 2024-09-20
 * 
 */

#include "pmm_mg_trajectory3d.hpp"

namespace pmm {

    PMM_MG_Trajectory3D::PMM_MG_Trajectory3D() : _exists(false), n_segments(0) {}


    PMM_MG_Trajectory3D::PMM_MG_Trajectory3D(const int n_path_segments) : n_segments(n_path_segments) {
        _dt.resize(3 * n_segments);
        _dtdv.resize(2 * n_segments);
        _p.resize(3 *n_segments +1);
        _v.resize(3 *n_segments +1);
        _a.resize(3 * n_segments);
        _traj_type_vec.resize(n_segments);
        _sol_idx.resize(n_segments);
        _used_grad.resize(n_segments);
        _update_vel_flag.resize(n_segments);
        _ax_min_traj_duration.resize(n_segments);
        n_segments_filled = 0;
    }

    PMM_MG_Trajectory3D::PMM_MG_Trajectory3D(std::vector<Vector<3>> p, std::vector<Vector<3>> v,
                                             Vector<3> a_max, Vector<3> v_max) : n_segments(p.size()-1){
        _dt.resize(3 * n_segments);
        _dtdv.resize(2 * n_segments);
        _p.resize(3 *n_segments +1);
        _v.resize(3 *n_segments +1);
        _a.resize(3 * n_segments);
        _traj_type_vec.resize(n_segments);
        _sol_idx.resize(n_segments);
        _update_vel_flag.resize(n_segments);
        _ax_min_traj_duration.resize(n_segments);
        _used_grad.resize(n_segments);
        n_segments_filled = 0;
        _a_max = a_max;
        _a_min = -a_max;
        _v_max = v_max;
        int succ = 1;

        // fill all trajectory segments
        for(int i=0; i<n_segments; i++){

            QuadState start_conf;
            start_conf.setZero();
            start_conf.p = p[i];
            start_conf.v = v[i];
            QuadState end_conf;
            end_conf.setZero();
            end_conf.p = p[i+1];
            end_conf.v = v[i+1];

            PointMassTrajectory3D segment_tr(start_conf, end_conf, a_max, v_max, true, true);    // equalize time

            // fill the computed trajectory in to the trajectory segment
            succ *= fillTrajectorySegment(&segment_tr, i);
        }
        if (succ == 0){
            std::cout << "ERROR: Unsuccessful PMM_MG_Trajectory3D initialization" << std::endl;
        } else {
            _exists = true;
        }
    }

    PMM_MG_Trajectory3D::PMM_MG_Trajectory3D(std::vector<Vector<3>> p, std::vector<Vector<3>> v,
                                             Vector<3> a_min, Vector<3> a_max, Vector<3> v_max) : n_segments(p.size()-1){
        _dt.resize(3 * n_segments);
        _dtdv.resize(2 * n_segments);
        _p.resize(3 *n_segments +1);
        _v.resize(3 *n_segments +1);
        _a.resize(3 * n_segments);
        _traj_type_vec.resize(n_segments);
        _sol_idx.resize(n_segments);
        _update_vel_flag.resize(n_segments);
        _ax_min_traj_duration.resize(n_segments);
        _used_grad.resize(n_segments);
        n_segments_filled = 0;
        _a_max = a_max;
        _a_min = a_min;
        _v_max = v_max;
        int succ = 1;

        bool equalize_times = true;
        bool calculate_gradients = true;

        // fill all trajectory segments
        for(int i=0; i<n_segments; i++){
            QuadState start_conf;
            start_conf.setZero();
            start_conf.p = p[i];
            start_conf.v = v[i];
            QuadState end_conf;
            end_conf.setZero();
            end_conf.p = p[i+1];
            end_conf.v = v[i+1];

            PointMassTrajectory3D segment_tr(start_conf, end_conf, a_max, a_min, v_max, equalize_times, calculate_gradients);

            // fill the computed trajectory in to the trajectory segment
            succ *= fillTrajectorySegment(&segment_tr, i);
        }
        if (succ == 0){
            std::cout << "ERROR: Unsuccessful PMM_MG_Trajectory3D initialization" << std::endl;
        } else {
            _exists = true;
        }
    }

    PMM_MG_Trajectory3D::PMM_MG_Trajectory3D(std::vector<Vector<3>> p, std::vector<Vector<3>> v,
                                             Vector<3> a_min, Vector<3> a_max, Vector<3> v_max, bool equalize_axes) : n_segments(p.size()-1){
        _dt.resize(3 * n_segments);
        _dtdv.resize(2 * n_segments);
        _p.resize(3 *n_segments +1);
        _v.resize(3 *n_segments +1);
        _a.resize(3 * n_segments);
        _traj_type_vec.resize(n_segments);
        _sol_idx.resize(n_segments);
        _update_vel_flag.resize(n_segments);
        _ax_min_traj_duration.resize(n_segments);
        _used_grad.resize(n_segments);
        n_segments_filled = 0;
        _a_max = a_max;
        _a_min = a_min;
        _v_max = v_max;
        int succ = 1;

        bool calculate_gradients = true;

        // fill all trajectory segments
        for(int i=0; i<n_segments; i++){
            QuadState start_conf;
            start_conf.setZero();
            start_conf.p = p[i];
            start_conf.v = v[i];
            QuadState end_conf;
            end_conf.setZero();
            end_conf.p = p[i+1];
            end_conf.v = v[i+1];

            PointMassTrajectory3D segment_tr(start_conf, end_conf, a_max, a_min, v_max, equalize_axes, calculate_gradients);

            // fill the computed trajectory in to the trajectory segment
            succ *= fillTrajectorySegment(&segment_tr, i);
        }
        if (succ == 0){
            std::cout << "ERROR: Unsuccessful PMM_MG_Trajectory3D initialization" << std::endl;
        } else {
            _exists = true;
        }
    }

    PMM_MG_Trajectory3D::PMM_MG_Trajectory3D(std::vector<Vector<3>> p, std::vector<Vector<3>> v, const Scalar max_acc_norm, const Scalar max_vel_norm,
                                             const int td_max_iter, const Scalar td_acc_norm_precision_limit, const bool drag, const bool debug) : n_segments(p.size()-1){
        _dt.resize(3 * n_segments);
        _dtdv.resize(2 * n_segments);
        _p.resize(3 *n_segments +1);
        _v.resize(3 *n_segments +1);
        _a.resize(3 * n_segments);
        _traj_type_vec.resize(n_segments);
        _sol_idx.resize(n_segments);
        _update_vel_flag.resize(n_segments);
        _ax_min_traj_duration.resize(n_segments);
        _used_grad.resize(n_segments);
        n_segments_filled = 0;
        _thrust_limit_flag = true;
        int succ = 1;

        _max_acc_norm = max_acc_norm;
        _max_vel_norm = max_vel_norm;
        _v_max = Vector<3>::Constant(max_vel_norm / sqrt(3));
        _TD_max_iter = td_max_iter;
        _TD_norm_precision = td_acc_norm_precision_limit;

        bool equalize_times = true;
        bool calculate_gradients = true;
        bool bang_singular_bang = false;

        // fill all trajectory segments
        for(int i=0; i<n_segments; i++){
            QuadState start_conf;
            start_conf.setZero();
            start_conf.p = p[i];
            start_conf.v = v[i];
            QuadState end_conf;
            end_conf.setZero();
            end_conf.p = p[i+1];
            end_conf.v = v[i+1];

            PointMassTrajectory3D segment_tr(start_conf, end_conf, max_acc_norm, max_vel_norm, drag, td_max_iter, td_acc_norm_precision_limit, debug);

            // fill the computed trajectory in to the trajectory segment
            succ *= fillTrajectorySegment(&segment_tr, i); 
        }
        if (succ == 0){
            std::cout << "ERROR: Unsuccessful PMM_MG_Trajectory3D initialization" << std::endl;
        } else {
            _exists = true;
        }
    }

    PMM_MG_Trajectory3D::PMM_MG_Trajectory3D(std::vector<Vector<3>> waypoints, const Vector<3> start_vel, const Vector<3> end_vel,
                                            const Scalar max_acc_norm, const Scalar max_vel_norm, const Scalar dT_precision, const int GD_max_iter, 
                                            const Scalar GD_alpha, const Scalar GD_alpha_red_fctr, const Scalar GD_alpha_min_thr,
                                            const Scalar TD_max_iter, const Scalar TD_acc_precision, const bool scnd_round_optim,
                                            const int GD_max_iter_2nd_rnd, const Scalar GD_alpha_2nd_rnd, const Scalar GD_alpha_red_fctr_2nd_rnd,
                                            const Scalar GD_alpha_min_thr_2nd_rnd, const bool drag, const bool debug) : n_segments(waypoints.size()-1){
        const std::vector<Vector<3>> wpts = waypoints;

        // initial values
        Vector<3> max_per_axis_acc_vec;
        Vector<3> min_per_axis_acc_vec;
        std::tie(min_per_axis_acc_vec, max_per_axis_acc_vec) = compute_per_axis_acc_limits_from_limited_thrust(max_acc_norm);

        Scalar max_vel_per_axis = max_vel_norm / sqrt(3);
        Vector<3> max_per_axis_vel_vec = Vector<3>::Constant(max_vel_per_axis);
        std::vector<Vector<3>> init_vel = compute_initial_velocities(waypoints, start_vel, end_vel, max_per_axis_acc_vec[0], max_vel_per_axis);

        // compute and optimize trajectory
        PMM_MG_Trajectory3D mp_tr(waypoints, init_vel, min_per_axis_acc_vec, max_per_axis_acc_vec, max_per_axis_vel_vec);
        if (debug){std::cout << "Original trajectory duration: " << mp_tr.duration() << " s." << std::endl;}
        mp_tr.optimize_velocities_at_positions(GD_alpha, GD_alpha_red_fctr, GD_alpha_min_thr, GD_max_iter, dT_precision, drag);

        // recompute trajectory using limited thrust decomposition
        std::vector<Vector<3>> v = mp_tr.get_current_velocities();
        PMM_MG_Trajectory3D mp_tr_o(wpts, v, max_acc_norm, max_vel_norm, TD_max_iter, TD_acc_precision, drag, debug);

        if (debug){std::cout << "Optimized trajectory duration after first optimization run with thrust decomposition: " << mp_tr_o.duration()  << " s." << std::endl;}

        // Second optimization run
        if (scnd_round_optim){
            mp_tr_o.optimize_velocities_at_positions(GD_alpha_2nd_rnd, GD_alpha_red_fctr_2nd_rnd, GD_alpha_min_thr_2nd_rnd,
                                                    GD_max_iter_2nd_rnd, dT_precision, drag);
            if (debug){std::cout << "Optimized trajectory duration after second optimization run with thrust decomposition: " << mp_tr_o.duration()  << " s." << std::endl<< std::endl;}
        }
        
        this->copy_trajectory(mp_tr_o);
    }

    int PMM_MG_Trajectory3D::optimize_velocities_at_positions(Scalar gd_step, Scalar gd_step_reduction_factor, Scalar gd_step_min_threshold,
                                            int max_iter, Scalar dT_threshold, const bool drag, int single_axis_optim_idx, bool debug,
                                            bool export_trajectory_to_file, std::string file_name, bool export_sampled_trajectories){
        
               
        // Check inputs
        if (gd_step_reduction_factor < 0. or gd_step_reduction_factor >= 1.){
            std::cout << "ERROR: Invalid gradient method reduction factor value." << std::endl;
        }

        if (gd_step_min_threshold >= gd_step){
            std::cout << "ERROR: Invalid gradient method minimum step threshold." << std::endl;
        }
        // Apply GD to adjust velocities to find a trajectory with minimal duration ------------------------------------------/
        std::ofstream outFile;

        if (export_trajectory_to_file){
            // initialize used_grad
            for (int wp=0; wp<_used_grad.size();wp++){
                _used_grad[wp] = Vector<3>(0.,0.,0.);
            }
            // export
            outFile = std::ofstream(file_name);
            this->exportTrajectoryToCSV(&outFile, true);
        }
        if(export_sampled_trajectories){
            std::stringstream ss;
            ss << "scripts/trajectory_data/dp_vo_sampled_trajectory_iter_"<<0<<".csv";
            this->sample_and_export_trajectory(1e-3, ss.str());
        }
        
        auto T_last = this->duration();
        auto T_new = T_last;
        double dT = 0;

        if(debug){
            std::cout << "Original trajectory time: " << T_last << std::endl;
        }

        // GD variables
        auto v_last = this->_v;;
        auto v_new = v_last;
        Vector<3> gd_step_vec = Vector<3>(gd_step, gd_step, gd_step);

        // Select if the time should be equalized throughout the GD optimization
        bool equalize_time = true;

        int iter = 0;
        for (iter=0; iter < max_iter; iter++) {
            if (debug){
                std::cout << std::endl << "Iteration: " << iter << std::endl;
            }

            // check preconditions
            if (!this->_exists){
                std::cout << "ERROR: PMM_MG_TRAJECTORY is infeasible, terminating." << iter << std::endl << std::endl;
                break;
            }

            // fill for first segment
            if (debug || export_trajectory_to_file){
                this->_used_grad[0] = Vector<3>(0.,0.,0.);
            }

            // update the velocities at the trajectory segment switch points using GD
            v_new = v_last;

            // reset update velocity flag for each iteration
            for (int seg = 1;seg < this->n_segments;seg++){
                this->_update_vel_flag[seg] = Eigen::Matrix<bool, 3, 1>(true, true, true);
            }

            int seg_idx = 0;
            for (int seg = 1; seg < (this->n_segments); seg++){
                if (iter % 2 == 1){
                    seg_idx = (this->n_segments)-seg;
                } else {
                    seg_idx = seg;
                }

                if (debug) {
                    std::cout << "Segment " << seg_idx << std::endl;
                }

                Eigen::Matrix<int, 3, 1> traj_type_1 = this->_traj_type_vec[seg_idx-1];
                Eigen::Matrix<int, 3, 1> traj_type_2 = this->_traj_type_vec[seg_idx];

                Eigen::Matrix<bool, 3, 1> axis_vel_updated_flag(false, false, false);

                // handle each axis separately
                for (int ax_idx = 0; ax_idx<3; ax_idx++){
                    // single axis optimization
                    if (single_axis_optim_idx >= 0 and ax_idx != single_axis_optim_idx){
                        continue;
                    }
                    if (!this->_update_vel_flag[seg_idx][ax_idx] or traj_type_1[ax_idx] == ZERO_TRAJ or traj_type_2[ax_idx] == ZERO_TRAJ){
                        // Do not update velocity
                        continue;
                    }                     
                    // UPDATE VELOCITIES BASED ON TWO SEGMENT SCOPE -------------------------------------------------------------------/
                    if (traj_type_1[ax_idx] == SYNC_TRAJ and traj_type_2[ax_idx] == SYNC_TRAJ){
                        // SYNC-SYNC -> no need for velocity update
                        if (debug){
                            std::cout << ax_idx << ": SYNC-SYNC" << std::endl;
                            std::cout << "v_new: " << v_last[(3*seg_idx)][ax_idx] << std::endl;
                        }
                        if (debug || export_trajectory_to_file){
                            this->_used_grad[seg_idx][ax_idx] = 0.;
                        }
                        continue;

                    } else if (traj_type_1[ax_idx] == MIN_TRAJ and traj_type_2[ax_idx] == MIN_TRAJ) {
                        // MIN-MIN -> update based on gradients calculated when computing trajectory (dtdve, dtdvs)                        
                        // Scalar dtdve_1 = this->_dtdv[(2*seg_idx)-1][ax_idx];
                        // Scalar dtdvs_2 = this->_dtdv[(2*seg_idx)][ax_idx];
                        // Scalar dtdv = dtdve_1 + dtdvs_2;
                        Scalar dtdv; Scalar v12_min; Scalar v12_max;
                        std::tie(dtdv, v12_min, v12_max) = this->get_MIN_MIN_grad_and_bounds(ax_idx, seg_idx);

                        Scalar v_last_tmp = v_last[(3*seg_idx)][ax_idx];
                        Scalar v_new_tmp;

                        v_new_tmp = v_last_tmp - (gd_step_vec[ax_idx] * dtdv);
                        
                        // Check velocity bounds and clip
                        if (v_new_tmp > v12_max){
                            v_new_tmp = v12_max;
                        } else if (v_new_tmp < v12_min) {
                            v_new_tmp = v12_min;
                        }

                        v_new[(3*seg_idx)][ax_idx] = v_new_tmp;
                        axis_vel_updated_flag[ax_idx] = true;

                        if (debug){
                            std::cout << ax_idx << ": MIN-MIN" << std::endl;
                            std::cout << "v_new: " << v_new[(3*seg_idx)][ax_idx] << std::endl;
                        }
                        if (debug || export_trajectory_to_file){
                            this->_used_grad[seg_idx][ax_idx] =  dtdv;
                        }
                    } else if (traj_type_1[ax_idx] == MIN_TRAJ and traj_type_2[ax_idx] == SYNC_TRAJ){
                        // MIN-SYNC -> update based on gradient of sum of absolute values of time for two consecutive segments
                        // Get velocity update bounds
                        Scalar v12_min_1, v12_max_1;
                        Scalar v_last_tmp = v_last[(3*seg_idx)][ax_idx];
                        std::tie(v12_min_1, v12_max_1) = get_MIN_SYNC_vel_bounds(seg_idx, ax_idx);

                        // Check bounds if an update is needed
                        if ((fabs(v_last_tmp-v12_min_1) < PRECISION_PMM_VALUES) or (fabs(v_last_tmp-v12_max_1) < PRECISION_PMM_VALUES)){
                            // No need for update
                            if (debug){
                                std::cout << "Bounds, no update!!" << std::endl;
                            }
                            continue;
                        }

                        Scalar dtdv; Scalar v12_min_2; Scalar v12_max_2;
                        std::tie(dtdv, v12_min_2, v12_max_2) = this->get_MIN_SYNC_grad_and_bounds(ax_idx, seg_idx-1);

                        // Check bounds given by the gradient
                        if ((fabs(v_last_tmp-v12_min_2) < PRECISION_PMM_VALUES) or (fabs(v_last_tmp - v12_max_2) < PRECISION_PMM_VALUES)){
                            // No need for update
                            if (debug){
                                std::cout << "Bounds, no update!!" << std::endl;
                            }
                            continue;
                        }
                        Scalar v_new_tmp = v_last_tmp - (gd_step_vec[ax_idx] * dtdv);

                        // clip new velocity based on SYNC segment
                        if (v_new_tmp > std::min(v12_max_1, v12_max_2)){
                            v_new_tmp = std::min(v12_max_1, v12_max_2);
                        } else if (v_new_tmp < std::max(v12_min_1, v12_min_2)) {
                            v_new_tmp = std::max(v12_min_1, v12_min_2);
                        }
                        // update
                        v_new[(3*seg_idx)][ax_idx] = v_new_tmp;
                        axis_vel_updated_flag[ax_idx] = true;

                        if (debug){
                            std::cout << ax_idx << ": MIN-SYNC" << std::endl;
                            std::cout << "v_new: " << v_new_tmp << std::endl;
                        }
                        if (debug || export_trajectory_to_file){
                            this->_used_grad[seg_idx][ax_idx] =  dtdv;
                        }
                    } else if (traj_type_1[ax_idx] == SYNC_TRAJ and traj_type_2[ax_idx] == MIN_TRAJ){
                        // SYNC-MIN-> update based on gradient of sum of absolute values of time for two consecutive segments

                        Scalar v12_min_1 ,v12_max_1;
                        Scalar v_last_tmp = v_last[(3*seg_idx)][ax_idx];
                        std::tie(v12_min_1, v12_max_1) = get_SYNC_MIN_vel_bounds(seg_idx-1, ax_idx);

                        if ((fabs(v_last_tmp-v12_min_1) < PRECISION_PMM_VALUES) or (fabs(v_last_tmp-v12_max_1) < PRECISION_PMM_VALUES)){
                            // No need for update
                            if (debug){
                                std::cout << "Bounds, no update!!" << std::endl;
                            }
                            continue;
                        }
                        
                        Scalar dtdv; Scalar v12_min_2; Scalar v12_max_2;
                        std::tie(dtdv, v12_min_2, v12_max_2) = this->get_SYNC_MIN_grad_and_bounds(ax_idx, seg_idx);

                        // Check bounds given by the gradient
                        if ((fabs(v_last_tmp-v12_min_2) < PRECISION_PMM_VALUES) or (fabs(v_last_tmp-v12_max_2) < PRECISION_PMM_VALUES)){
                            // No need for update
                            if (debug){
                                std::cout << "Bounds, no update!!" << std::endl;
                            }
                            continue;
                        }
                        // dtdv = this->_dtdv[(2*seg_idx)-1][ax_idx];
                        Scalar v_new_tmp = v_last_tmp - (gd_step_vec[ax_idx] * dtdv);

                        // clip new velocity based on SYNC segment
                        if (v_new_tmp > std::min(v12_max_1, v12_max_2)){
                            v_new_tmp = std::min(v12_max_1, v12_max_2);
                        } else if (v_new_tmp < std::max(v12_min_1, v12_min_2)) {
                            v_new_tmp = std::max(v12_min_1, v12_min_2);
                        }
                        // update
                        v_new[(3*seg_idx)][ax_idx] = v_new_tmp;
                        axis_vel_updated_flag[ax_idx] = true;

                        if (debug){
                            std::cout << ax_idx << ": SYNC-MIN" << std::endl;
                            std::cout <<  " v_new: " << v_new_tmp << std::endl;
                        }
                        if (debug || export_trajectory_to_file){
                            this->_used_grad[seg_idx][ax_idx] =  dtdv;
                        }
                    } else {
                        std::cout << "ERROR: SYNC flag error." << iter << std::endl << std::endl;
                        return -1;
                    }
                }
                
                // Keep track of segment duration for comparison with updated values
                Scalar seg1_time_old = this->segment_duration(seg_idx-1, 0);
                Scalar seg2_time_old = this->segment_duration(seg_idx, 0);

                for (int a=0;a<3;a++) {
                    if (fabs(v_new[3*seg_idx](a)) > _v_max(a)) {
                        v_new[3*seg_idx](a) = std::copysign( _v_max(a), v_new[3*seg_idx](a) );
                    }
                }
                // if (v_new[3*seg_idx].norm() > _v_max.norm()) {
                //     v_new[3*seg_idx] = v_new[3*seg_idx].normalized() * _v_max.norm();
                // }         

                // Update only the segments that depend on the updated velocities --------------------------------------------------------/
                int succ = this->updateTrajectoryGivenVelocities(v_new[3*seg_idx], seg_idx, equalize_time, true, drag);

                // Check if the trajectory time has not been increased due to the update
                if ((this->segment_duration(seg_idx-1, 0) + this->segment_duration(seg_idx, 0)) > (seg1_time_old + seg2_time_old + PRECISION_PMM_VALUES) or (succ==-1)){
                    // two segment time increased => one or both segment times increased 
                    // => the slowest axis for both of the segments must be checked
                    if (debug) {
                        std::cout << "succ " << succ << " time increase " << (this->segment_duration(seg_idx-1, 0) + this->segment_duration(seg_idx, 0)) << " > " << (seg1_time_old + seg2_time_old + PRECISION_PMM_VALUES) << std::endl;
                    }
                    bool change_occurred = false;
                    if (succ==1){
                        for (int ax=0;ax<3;ax++){
                            // firs segment
                            if (this->_ax_min_traj_duration[seg_idx-1][ax] > seg1_time_old and axis_vel_updated_flag[ax]){    // !this->_sync_traj_flag[seg_idx-1][ax] and 
                                gd_step_vec[ax] = gd_step_reduction_factor * gd_step_vec[ax];
                                change_occurred = true;

                                if (gd_step_vec[ax] < gd_step_min_threshold) {
                                    // stop velocity update for this axis
                                    gd_step_vec[ax] = gd_step;
                                    this->_update_vel_flag[seg_idx][ax] = false;
                                }

                                if (debug){
                                    std::cout << "INFO: Reduced gd_step for axis " << ax << " to gd_step=" << gd_step_vec[ax] << " for segment " << seg_idx << std::endl;
                                    if (!this->_update_vel_flag[seg_idx][ax]){
                                        std::cout << "INFO: Stopping velocity update for axis" << ax << ", gd_step below threshold " << gd_step_min_threshold << std::endl;
                                    }
                                }
                            } else if (this->_ax_min_traj_duration[seg_idx][ax] > seg2_time_old and axis_vel_updated_flag[ax]) {      // !this->_sync_traj_flag[seg_idx][ax] and 
                                gd_step_vec[ax] = gd_step_reduction_factor * gd_step_vec[ax];
                                change_occurred = true;
                                
                                if (gd_step_vec[ax] < gd_step_min_threshold) {
                                    // stop velocity update for this axis
                                    gd_step_vec[ax] = gd_step_min_threshold;
                                    this->_update_vel_flag[seg_idx][ax] = false;
                                }

                                if (debug){
                                    std::cout << "INFO: Reduced gd_step for axis " << ax << " to gd_step=" << gd_step_vec[ax] << " for segment " << seg_idx << std::endl;
                                    if (!this->_update_vel_flag[seg_idx][ax]){
                                        std::cout << "INFO: Stopping velocity update for axis" << ax << ", gd_step below threshold " << gd_step_min_threshold << std::endl;
                                    }
                                }
                            }
                        }
                    }

                    if (!change_occurred){
                        // Update of other axis caused time increase in the rest of the axis -> reduce factor in all axis that were updated
                        if (debug){
                            std::cout << "WARNING: Velocity update in one axis caused time increase in other axis" << std::endl;

                            std::cout << seg1_time_old << " s -> " << this->segment_duration(seg_idx-1, 0) << " s" << std::endl;
                            std::cout << seg2_time_old << " s -> " << this->segment_duration(seg_idx, 0) << " s" << std::endl;
                        }
                        for (int ax=0;ax<3;ax++){
                            if (axis_vel_updated_flag[ax]){
                                gd_step_vec[ax] = gd_step_reduction_factor * gd_step_vec[ax];

                                if (gd_step_vec[ax] < gd_step_min_threshold) {
                                    // stop velocity update for this axis
                                    gd_step_vec[ax] = gd_step_min_threshold;
                                    this->_update_vel_flag[seg_idx][ax] = false;
                                }

                                 if (debug){
                                    std::cout << "INFO: Reduced gd_step for axis " << ax << " to gd_step=" << gd_step_vec[ax] << " for segment " << seg_idx << std::endl;
                                 }
                            }
                        }
                    }
                    // start from previous values
                    v_new[3*seg_idx] = v_last[3*seg_idx];
                    this->updateTrajectoryGivenVelocities(v_new[3*seg_idx], seg_idx, equalize_time, true, drag);
                    seg--;
                } else {
                    gd_step_vec = Vector<3>(gd_step, gd_step, gd_step);
                }
            }

            // write out convergence data to a csv file
            if(export_trajectory_to_file){
                this->exportTrajectoryToCSV(&outFile, true);
            }
            if(export_sampled_trajectories){
                std::stringstream ss;
                ss << "scripts/trajectory_data/dp_vo_sampled_trajectory_iter_"<<iter+1<<".csv";
                this->sample_and_export_trajectory(1e-3, ss.str());
            }

            // check terminating conditions
            T_new = this->duration();
            dT = fabs(T_new-T_last);

            if (debug){
                std::cout << "dT: " << dT << std::endl;
            }

            if (dT < dT_threshold){
                if(debug){
                    std::cout << "Termination, dT_traj smaller then threshold" << std::endl;
                }
            break;
            }

            if (T_new > T_last){
                this->updateTrajectoryGivenVelocities(v_last, equalize_time, true, drag);
                if (debug){
                    std::cout << "Termination, T_traj would increase" << std::endl;
                }
                if (export_trajectory_to_file){
                    this->exportTrajectoryToCSV(&outFile, true);
                }
                break;
            }

            T_last = T_new;
            v_last = v_new;
        }

        if(debug){
            std::cout << "Optimized trajectory time: " << T_new << std::endl << "Iteration count: " << iter+1 << std::endl;
        }
        return iter;
    }

    bool PMM_MG_Trajectory3D::checkTrajectory(bool drag, bool verbose){
        if (n_segments_filled != n_segments){
            std::cout << "ERROR: Cannot check incomplete trajectory" << std::endl;
            return false;
        }

        Scalar precision = 1e-3;

        // -------------------------------------------/
        Scalar avg_thrust_usage = 0;
        Scalar max_thrust = 0;
        Scalar min_thrust = _max_acc_norm;

        for (int seg_idx=0; seg_idx< n_segments; seg_idx++){

            if (_thrust_limit_flag){
                QuadState start_conf;
                QuadState end_conf;
                start_conf.setZero();
                start_conf.p = _p[3*(seg_idx)];
                start_conf.v = _v[(3*seg_idx)];
                end_conf.setZero();
                end_conf.p = _p[3*(seg_idx)+3];
                end_conf.v = _v[3*(seg_idx)+3];

                auto start_acc = _a[(3*seg_idx)];
                auto end_acc = _a[3*(seg_idx)+2];
                Vector<3> max_per_axis_acc_vec;
                Vector<3> min_per_axis_acc_vec;

                for (int ax=0; ax<3;ax++){
                    if (start_acc[ax] < end_acc[ax]){
                        min_per_axis_acc_vec[ax] = start_acc[ax];
                        max_per_axis_acc_vec[ax] = end_acc[ax];
                    } else {
                        min_per_axis_acc_vec[ax] = end_acc[ax];
                        max_per_axis_acc_vec[ax] = start_acc[ax];
                    }
                }

                PointMassTrajectory3D pmm3d(start_conf, end_conf, max_per_axis_acc_vec, min_per_axis_acc_vec, _v_max, true, true);

                Vector<3> switch_times = pmm3d.switch_times();
                // std::sort(switch_times.begin(), switch_times.end());
                std::sort(switch_times.data(), switch_times.data() + switch_times.size());


                Scalar duration = pmm3d.time();


                Vector<3> acc1 = pmm3d.acc_in_time(switch_times(0));
                Vector<3> acc2 = pmm3d.acc_in_time(switch_times(1));
                Vector<3> acc3 = pmm3d.acc_in_time(switch_times(2));
                Vector<3> acc4 = pmm3d.end_acc();

                Scalar thrust1;
                Scalar thrust2;
                Scalar thrust3;
                Scalar thrust4;
                if (drag){
                    Vector<3> drag1 = pmm3d.drag(acc1 - GVEC, switch_times(0));
                    Vector<3> drag2 = pmm3d.drag(acc2 - GVEC, switch_times(1));
                    Vector<3> drag3 = pmm3d.drag(acc3 - GVEC, switch_times(2));
                    Vector<3> drag4 = pmm3d.drag(acc4 - GVEC, pmm3d.time());

                    thrust1 = (acc1 - drag1 - GVEC).norm();
                    thrust2 = (acc2 - drag2 - GVEC).norm();
                    thrust3 = (acc3 - drag3 - GVEC).norm();
                    thrust4 = (acc4 - drag4 - GVEC).norm();
                } else {
                    thrust1 = (acc1 - GVEC).norm();
                    thrust2 = (acc2 - GVEC).norm();
                    thrust3 = (acc3 - GVEC).norm();
                    thrust4 = (acc4 - GVEC).norm(); 
                }

                Scalar largest_thrust = std::max(std::max(thrust1, thrust2), std::max(thrust3, thrust4));
                Scalar smallest_thrust = std::min(std::min(thrust1, thrust2), std::min(thrust3, thrust4));

                if (largest_thrust > (_max_acc_norm + _TD_norm_precision)){
                    std::cout << "ERROR: Invalid trajectory, acceleration norm above maximal values." << std::endl;
                    return false;
                }

                // Thrust usage stats
                Scalar k1 = switch_times(0)/duration;
                Scalar k2 = (switch_times(1)-switch_times(0))/duration;
                Scalar k3 = (switch_times(2)-switch_times(1))/duration;
                Scalar k4 = (duration-switch_times(2))/duration;
                // std::cout << k1+k2+k3+k4 << std::endl;

                // averaging
                Scalar tmp_avg_thrust = k1*thrust1 + k2*thrust2 + k3*thrust3 + k4*thrust4; // averaging over segment
                avg_thrust_usage += (tmp_avg_thrust/_max_acc_norm)/n_segments; // averaging over whole trajectory
                if (largest_thrust > max_thrust){
                    max_thrust = largest_thrust;
                }
                if (smallest_thrust < min_thrust){
                    min_thrust = smallest_thrust;
                }

            }


            for (int ax=0; ax<3;ax++){
                Scalar p0 = _p[3*seg_idx][ax];
                Scalar p11 = _p[(3*seg_idx)+1][ax];
                Scalar p12 = _p[(3*seg_idx)+2][ax];
                Scalar p13 = _p[(3*seg_idx)+3][ax];
                Scalar v0 = _v[3*seg_idx][ax];
                Scalar v11 = _v[(3*seg_idx)+1][ax];
                Scalar v12 = _v[(3*seg_idx)+2][ax];
                Scalar v13 = _v[(3*seg_idx)+3][ax];
                Scalar a11 = _a[(3*seg_idx)][ax];
                Scalar a12 = _a[(3*seg_idx)+1][ax];
                Scalar a13 = _a[(3*seg_idx)+2][ax];
                Scalar t11 = _dt[(3*seg_idx)][ax];
                Scalar t12 = _dt[(3*seg_idx)+1][ax];
                Scalar t13 = _dt[(3*seg_idx)+2][ax];

                // check acceleration
                if (!_thrust_limit_flag){
                    if ((a11 > (_a_max[ax]+precision) or a11 < (_a_min[ax]-precision)) or
                        (a12 > (_a_max[ax]+precision) or a12 < (_a_min[ax]-precision)) or
                        (a13 > (_a_max[ax]+precision) or a13 < (_a_min[ax]-precision))){
                        std::cout << "ERROR: Invalid trajectory, acceleration above maximal values." << std::endl;
                        return false;
                    }
                }
                

                // check times
                if (t11 < 0 or t12 < 0 or t13 < 0){
                    std::cout << "ERROR: Invalid trajectory, negative time values." << std::endl;
                    return false;
                }

                // integrate and check results
                Scalar p11_s = p0 + v0*t11+0.5*a11*pow(t11, 2);
                Scalar p12_s = p11 + v11*t12+0.5*a12*pow(t12, 2);
                Scalar p13_s = p12 + v12*t13+0.5*a13*pow(t13, 2);
                Scalar v11_s = v0 + a11*t11;
                Scalar v12_s = v11 + a12*t12;
                Scalar v13_s = v12 + a13*t13;

                if (fabs(p11-p11_s) > precision or fabs(p12-p12_s) > precision or
                    fabs(p13-p13_s) > precision or fabs(v11-v11_s) > precision or
                    fabs(v12-v12_s) > precision or fabs(v13-v13_s) > precision) {
                        std::cout << "ERROR: Invalid trajectory, positions or velocities do not match the rest of the trajectory data." << std::endl;
                        return false;
                    }
            }
            // check if all axis have same segment duration
            if (fabs(segment_duration(seg_idx, 0) - segment_duration(seg_idx, 1)) > precision or 
                fabs(segment_duration(seg_idx, 0) - segment_duration(seg_idx, 2)) > precision or 
                fabs(segment_duration(seg_idx, 1) - segment_duration(seg_idx, 2)) > precision){
                    std::cout << "ERROR: Invalid trajectory, segment durations are not equal across all axis." << std::endl;
                    return false;
                }

        }

        if (verbose) {
            std::cout << "INFO: Trajectory is feasible and valid." << std::endl;
            std::cout << "INFO: Average thrust usage: "  << avg_thrust_usage << std::endl;
            std::cout << "INFO: Minimum thrust usage: "  << min_thrust/_max_acc_norm << std::endl;
            std::cout << "INFO: Maximum thrust usage: "  << max_thrust/_max_acc_norm << std::endl;
        }
        
        return true;
    }

    void PMM_MG_Trajectory3D::reset(){
        n_segments_filled = 0;
        _exists = false;
        _thrust_limit_flag = false;
    }

    int PMM_MG_Trajectory3D::fillTrajectorySegment(PointMassTrajectory3D *tr, int segment_idx, const bool override){
        if ((segment_idx < this->n_segments_filled) and (!override)){
            // this trajectory segment has been already defined
            std::cout << "ERROR:fillTrajectorySegment method failed, trajectory segment has already been defined." << std::endl;
            return 0;
        }

        if (!tr->exists()){
            // filling trajectory that is unfeasible
            std::cout << "ERROR:fillTrajectorySegment method failed due to invalid input trajectory." << std::endl;
            return 0;
        }

        if (segment_idx == 0){
            // first segment, fill also starting positions
            _p[0] = Vector<3>(tr->x_.p_(0), tr->y_.p_(0), tr->z_.p_(0));
            _v[0] = Vector<3>(tr->x_.v_(0), tr->y_.v_(0), tr->z_.v_(0));
        }

        // for all others it is assumed that the starting configuration is the same as
        // the previous one

        // p_(0) is the end configuration of the previous trajectory segment
        _p[(3*segment_idx)+1] = Vector<3>(tr->x_.p_(1), tr->y_.p_(1), tr->z_.p_(1));
        _p[(3*segment_idx)+2] = Vector<3>(tr->x_.p_(2), tr->y_.p_(2), tr->z_.p_(2));
        _p[(3*segment_idx)+3] = Vector<3>(tr->x_.p_(3), tr->y_.p_(3), tr->z_.p_(3));

        // v_(0) is the end velocity of the previous trajectory segment
        _v[(3*segment_idx)+1] = Vector<3>(tr->x_.v_(1), tr->y_.v_(1), tr->z_.v_(1));
        _v[(3*segment_idx)+2] = Vector<3>(tr->x_.v_(1), tr->y_.v_(1), tr->z_.v_(1));    // same as v_(1) even for bang-singular-bang
        _v[(3*segment_idx)+3] = Vector<3>(tr->x_.v_(2), tr->y_.v_(2), tr->z_.v_(2));

        _dt[(3*segment_idx)] = Vector<3>(tr->x_.t_(0), tr->y_.t_(0), tr->z_.t_(0));
        _dt[(3*segment_idx)+1] = Vector<3>(tr->x_.t_(1), tr->y_.t_(1), tr->z_.t_(1));
        _dt[(3*segment_idx)+2] = Vector<3>(tr->x_.t_(2), tr->y_.t_(2), tr->z_.t_(2));

        _dtdv[2*segment_idx] = Vector<3>(tr->x_.dt_dvs_, tr->y_.dt_dvs_, tr->z_.dt_dvs_);
        _dtdv[(2*segment_idx)+1] = Vector<3>(tr->x_.dt_dve_, tr->y_.dt_dve_, tr->z_.dt_dve_);

        _a[3*segment_idx] = Vector<3>(tr->x_.a_(0), tr->y_.a_(0), tr->z_.a_(0));
        _a[(3*segment_idx)+1] = Vector<3>(0, 0, 0);
        _a[(3*segment_idx)+2] = Vector<3>(tr->x_.a_(1), tr->y_.a_(1), tr->z_.a_(1));

        // variables for optimization
        _traj_type_vec[segment_idx] = Eigen::Matrix<int, 3, 1>(tr->x_.trajectory_type_, tr->y_.trajectory_type_, tr->z_.trajectory_type_);
        _sol_idx[segment_idx] = Eigen::Matrix<int, 3, 1>(tr->x_.sol_idx_, tr->y_.sol_idx_, tr->z_.sol_idx_);
        _ax_min_traj_duration[segment_idx] = tr->min_axis_trajectory_duration_;

        if (!override){
            this->n_segments_filled++;
            _update_vel_flag[segment_idx] = Eigen::Matrix<bool, 3, 1>(true, true, true);
        }
        return 1;
    }

    int PMM_MG_Trajectory3D::updateTrajectoryGivenVelocities(std::vector<Vector<3>> v_new, const bool equalize_times, const bool sync_traj, const bool drag){
        if (n_segments_filled != n_segments){
            // The trajectory is not complete, some data is missing
            std::cout << "ERROR:updateTrajectoryGivenVelocities method failed due to incomplete trajectory." << std::endl;
            return 0;
        }
        
        int succ = 1;
        this->reset();
        for (int i = 0; i < n_segments; i++) {
            // recompute trajectory segment and update trajectory values
            QuadState start_conf;
            start_conf.setZero();
            start_conf.p = _p[3*i];
            start_conf.v = v_new[3*i];
            QuadState end_conf;
            end_conf.setZero();
            end_conf.p = _p[(3*i)+3];
            end_conf.v = v_new[(3*i)+3];

            if (!sync_traj){
                PointMassTrajectory3D tr(start_conf, end_conf, _a_max, _v_max, equalize_times, true);    // equalize time
                if (!tr.exists()){
                    return -1;
                }
                succ *= fillTrajectorySegment(&tr, i);
            } else {
                if (_thrust_limit_flag){
                    PointMassTrajectory3D tr(start_conf, end_conf, _max_acc_norm, _max_vel_norm, drag, _TD_max_iter, _TD_norm_precision, false);
                    if (!tr.exists()){
                        return -1;
                    }
                    succ *= fillTrajectorySegment(&tr, i);
                } else {
                    PointMassTrajectory3D tr(start_conf, end_conf, _a_max, _v_max, equalize_times, true);
                    if (!tr.exists()){
                        return -1;
                    }
                    succ *= fillTrajectorySegment(&tr, i);
                }
            }
        }
        if (succ == 1){
            _exists = true;
        }

        return succ;
    }

    int PMM_MG_Trajectory3D::updateTrajectoryGivenVelocities(Vector<3> v_new, const int segment_idx, const bool equalize_times, const bool sync_traj, const bool drag){
        if (n_segments_filled != n_segments){
            // The trajectory is not complete, some data is missing
            std::cout << "ERROR:updateTrajectoryGivenVelocities method failed due to incomplete trajectory." << std::endl;
            return 0;
        }

        int succ = 1;
        _exists = false;
        // update previous and selected segment (end velocity of previous segment == start velocity of current segment == v_new)
        for (int i = 0; i < 2; i++) {
            // recompute trajectory segment and update trajectory values
            QuadState start_conf;
            QuadState end_conf;

            if (i==0){
                start_conf.setZero();
                start_conf.p = _p[3*(segment_idx-1)];
                start_conf.v = _v[3*(segment_idx-1)];
                end_conf.setZero();
                end_conf.p = _p[3*(segment_idx-1)+3];
                end_conf.v = v_new;
            } else {
                start_conf.setZero();
                start_conf.p = _p[3*(segment_idx)];
                start_conf.v = v_new;
                end_conf.setZero();
                end_conf.p = _p[(3*segment_idx)+3];
                end_conf.v = _v[(3*segment_idx)+3];
            }

            if (!sync_traj){
                PointMassTrajectory3D tr(start_conf, end_conf, _a_max, _v_max, equalize_times, true);    // equalize time
                if (!tr.exists()){
                    return -1;
                }
                succ *= fillTrajectorySegment(&tr, (segment_idx-1)+i, true);
            } else {
                if (_thrust_limit_flag){
                    PointMassTrajectory3D tr(start_conf, end_conf, _max_acc_norm, _max_vel_norm, drag, _TD_max_iter, _TD_norm_precision, false);
                    if (!tr.exists()){
                        return -1;
                    }
                    succ *= fillTrajectorySegment(&tr, (segment_idx-1)+i, true);
                } else {
                    PointMassTrajectory3D tr(start_conf, end_conf, _a_max, _a_min, _v_max, equalize_times, true);
                    if (!tr.exists()){
                        return -1;
                    }
                    succ *= fillTrajectorySegment(&tr, (segment_idx-1)+i, true);
                }
                
            }
        }
        if (succ == 1){
            _exists = true;
        }

        return succ;
    }

    std::tuple<Scalar, Scalar, Scalar>  PMM_MG_Trajectory3D::get_MIN_SYNC_grad_and_bounds(const int ax_idx, const int min_idx){
        // get solution index of the MIN segment

        // // get current parameters
        Scalar a11 = _a[3*min_idx][ax_idx];
        Scalar a12 = _a[3*min_idx+2][ax_idx];
        Scalar v0 = _v[3*min_idx][ax_idx];
        Scalar v11 = _v[3*min_idx+1][ax_idx];
        Scalar v12 = _v[3*min_idx+3][ax_idx];
        Scalar p0 = _p[3*min_idx][ax_idx];
        Scalar p12 = _p[3*min_idx+3][ax_idx];


        Scalar v21 = _v[3*(min_idx+1)+1][ax_idx];

        
        Scalar dtdv = this->_dtdv[(2*min_idx)+1][ax_idx];

        // Gradient bounds/ singularities 
        Scalar b_max = std::max(std::max(fabs(v11), fabs(v21)), _v_max(ax_idx));
        Scalar b_min = -b_max;
        std::vector<Scalar> b_vect;

        b_vect.push_back((sqrt(a12)*sqrt(pow(v0,2) - 2*a11*p0 + 2*a11*p12))/sqrt(a11));
        b_vect.push_back(-(sqrt(a12)*sqrt(pow(v0,2) - 2*a11*p0 + 2*a11*p12))/sqrt(a11));
        b_vect.push_back(-sqrt(pow(v0,2) - 2*a11*p0 + 2*a11*p12));
        b_vect.push_back(sqrt(pow(v0,2) - 2*a11*p0 + 2*a11*p12));
        b_vect.push_back(-sqrt(pow(v0,2) - 2*a12*p0 + 2*a12*p12));
        b_vect.push_back(sqrt(pow(v0,2) - 2*a12*p0 + 2*a12*p12));

        for (int i=0; i<6 ; i++){
            if (std::isfinite(b_vect[i])){
                if ((b_vect[i] < v12) and (b_vect[i] > b_min)){
                    b_min = b_vect[i];
                } else if ((b_vect[i] >= v12) and (b_vect[i] < b_max)){
                    b_max = b_vect[i];
                }
            }
        }

        return std::tuple(dtdv, b_min, b_max);
    }

    std::tuple<Scalar, Scalar, Scalar> PMM_MG_Trajectory3D::get_SYNC_MIN_grad_and_bounds(const int ax_idx, const int min_idx){
        // get solution index of the MIN segment

        // get current parameters
        Scalar a21 = _a[3*min_idx][ax_idx];
        Scalar a22 = _a[3*min_idx+2][ax_idx];
        Scalar v12 = _v[3*min_idx][ax_idx];
        Scalar v11 = _v[3*min_idx+1][ax_idx];
        Scalar v22 = _v[3*min_idx+3][ax_idx];
        Scalar p12 = _p[3*min_idx][ax_idx];
        Scalar p22 = _p[3*min_idx+3][ax_idx];


        Scalar v21 = _v[3*(min_idx-1)+1][ax_idx];


        Scalar dtdv = this->_dtdv[(2*min_idx)][ax_idx];
        
        // Gradient bounds/ singularities
        Scalar b_max = std::max(std::max(fabs(v11), fabs(v21)), _v_max(ax_idx));
        Scalar b_min = -b_max;
        std::vector<Scalar> b_vect;

        b_vect.push_back((sqrt(a21)*sqrt(pow(v22,2) + 2*a22*p12 - 2*a22*p22))/sqrt(a22));
        b_vect.push_back(-(sqrt(a21)*sqrt(pow(v22,2) + 2*a22*p12 - 2*a22*p22))/sqrt(a22));
        b_vect.push_back(-sqrt(pow(v22,2) + 2*a22*p12 - 2*a22*p22));
        b_vect.push_back(sqrt(pow(v22,2) + 2*a22*p12 - 2*a22*p22));
        b_vect.push_back(-sqrt(pow(v22,2) + 2*a21*p12 - 2*a21*p22));
        b_vect.push_back(sqrt(pow(v22,2) + 2*a21*p12 - 2*a21*p22));

        for (int i=0; i<6 ; i++){
            if (std::isfinite(b_vect[i])){
                if ((b_vect[i] < v12) and (b_vect[i] > b_min)){
                    b_min = b_vect[i];
                } else if ((b_vect[i] >= v12) and (b_vect[i] < b_max)){
                    b_max = b_vect[i];
                }
            }
        }
        
        return std::tuple(dtdv, b_min, b_max);
    }

     std::tuple<Scalar, Scalar, Scalar> PMM_MG_Trajectory3D::get_MIN_MIN_grad_and_bounds(const int ax_idx, const int second_min_idx){
        // get solution indexes of the MIN segments
        const int min_idx1 = second_min_idx-1;
        const int min_idx2 = second_min_idx;

        // get current parameters
        Scalar a11 = _a[3*min_idx1][ax_idx];
        Scalar a12 = _a[3*min_idx1+2][ax_idx];
        Scalar v0 = _v[3*min_idx1][ax_idx];
        Scalar v12 = _v[3*min_idx1+3][ax_idx];
        Scalar p0 = _p[3*min_idx1][ax_idx];
        Scalar p12 = _p[3*min_idx1+3][ax_idx];
        Scalar a21 = _a[3*min_idx2][ax_idx];
        Scalar a22 = _a[3*min_idx2+2][ax_idx];
        Scalar v22 = _v[3*min_idx2+3][ax_idx];
        Scalar p22 = _p[3*min_idx2+3][ax_idx];


        Scalar v11 = _v[3*min_idx1+1][ax_idx];
        Scalar v21 = _v[3*min_idx2+1][ax_idx];

        // Gradient
        Scalar dtdv = this->_dtdv[(2*min_idx2)-1][ax_idx] + this->_dtdv[(2*min_idx2)][ax_idx];
        
        // Gradient bounds/singularities
        Scalar b_max = std::max(std::max(fabs(v11), fabs(v21)), _v_max(ax_idx));
        Scalar b_min = -b_max;
        std::vector<Scalar> b_vect;

        b_vect.push_back((sqrt(a12)*sqrt(pow(v0,2) - 2*a11*p0 + 2*a11*p12))/sqrt(a11));
        b_vect.push_back(-(sqrt(a12)*sqrt(pow(v0,2) - 2*a11*p0 + 2*a11*p12))/sqrt(a11));
        b_vect.push_back(-sqrt(pow(v0,2) - 2*a11*p0 + 2*a11*p12));
        b_vect.push_back(sqrt(pow(v0,2) - 2*a11*p0 + 2*a11*p12));
        b_vect.push_back(-sqrt(pow(v0,2) - 2*a12*p0 + 2*a12*p12));
        b_vect.push_back(sqrt(pow(v0,2) - 2*a12*p0 + 2*a12*p12));

        b_vect.push_back((sqrt(a21)*sqrt(pow(v22,2) + 2*a22*p12 - 2*a22*p22))/sqrt(a22));
        b_vect.push_back(-(sqrt(a21)*sqrt(pow(v22,2) + 2*a22*p12 - 2*a22*p22))/sqrt(a22));
        b_vect.push_back(-sqrt(pow(v22,2) + 2*a22*p12 - 2*a22*p22));
        b_vect.push_back(sqrt(pow(v22,2) + 2*a22*p12 - 2*a22*p22));
        b_vect.push_back(-sqrt(pow(v22,2) + 2*a21*p12 - 2*a21*p22));
        b_vect.push_back(sqrt(pow(v22,2) + 2*a21*p12 - 2*a21*p22));

        for (int i=0; i<12 ; i++){
            if (std::isfinite(b_vect[i])){
                if ((b_vect[i] < v12) and (b_vect[i] > b_min)){
                    b_min = b_vect[i];
                } else if ((b_vect[i] >= v12) and (b_vect[i] < b_max)){
                    b_max = b_vect[i];
                }
            }
        }
        
        return std::tuple(dtdv, b_min, b_max);
    }

    std::tuple<Scalar, Scalar> PMM_MG_Trajectory3D::get_SYNC_MIN_vel_bounds(const int seg_idx, const int ax_idx){
        Scalar p0 = _p[3*seg_idx][ax_idx];
        Scalar p12 = _p[(3*seg_idx)+3][ax_idx];
        Scalar v0 = _v[(3*seg_idx)][ax_idx];
        Scalar v11 = _v[(3*seg_idx+1)][ax_idx];
        Scalar v12 = _v[(3*seg_idx+3)][ax_idx];
        Scalar T = segment_duration(seg_idx, ax_idx);
        Scalar sc = 1.0 - 1e-5;

        Scalar v21 = _v[(3*(seg_idx+1)+1)][ax_idx];

        Scalar v12_max = std::max(std::max(fabs(v11), fabs(v21)), _v_max(ax_idx));
        Scalar v12_min = -v12_max;

        if (!_thrust_limit_flag){
            std::vector<Vector<2>> test_acc_vec = {Vector<2>(_a_min[ax_idx], _a_max[ax_idx]),Vector<2>(_a_max[ax_idx], _a_min[ax_idx])};
            for (Vector<2> a : test_acc_vec) {
                Scalar a11 = a(0);
                Scalar a12 = a(1);

                Scalar sigma_1 = sqrt(sc*(a11 - a12)*(a11*sc*pow(T,2) + 2*v0*T + 2*p0 - 2*p12));
                Scalar t1_1 = (v0 - sigma_1 + T*a11*sc)/(sc*(a11 - a12)) - (v0 + T*a12*sc)/(sc*(a11 - a12));
                Scalar t1_2 = (v0 + sigma_1 + T*a11*sc)/(sc*(a11 - a12)) - (v0 + T*a12*sc)/(sc*(a11 - a12));
                Scalar t2_1 = (v0 + T*a11*sc)/(sc*(a11 - a12)) - (v0 - sigma_1 + T*a11*sc)/(sc*(a11 - a12));
                Scalar t2_2 = (v0 + T*a11*sc)/(sc*(a11 - a12)) - (v0 + sigma_1 + T*a11*sc)/(sc*(a11 - a12));

                if (std::isfinite(t1_1) and std::isfinite(t2_1) and t1_1 > -PRECISION_PMM_VALUES and t2_1 > -PRECISION_PMM_VALUES) {
                    Scalar v12 = v0 - sigma_1 + T*a11*sc;
                    if (v12 > 0. and v12 < v12_max){
                        v12_max = v12;
                    } else if (v12 < 0. and v12 > v12_min){
                        v12_min = v12;
                    }
                }
                if (std::isfinite(t1_2) and std::isfinite(t2_2) and t1_2 > -PRECISION_PMM_VALUES and t2_2 > -PRECISION_PMM_VALUES) {
                    Scalar v12 = v0 + sigma_1 + T*a11*sc;
                    if (v12 > 0. and v12 < v12_max){
                        v12_max = v12;
                    } else if (v12 < 0. and v12 > v12_min){
                        v12_min = v12;
                    }
                }
            }
        }

        return std::tuple(v12_min, v12_max);

    }

    std::tuple<Scalar, Scalar> PMM_MG_Trajectory3D::get_MIN_SYNC_vel_bounds(const int seg_idx, const int ax_idx){
        Scalar p0 = _p[3*seg_idx][ax_idx];
        Scalar p12 = _p[(3*seg_idx)+3][ax_idx];
        Scalar v0 = _v[(3*seg_idx)][ax_idx];
        Scalar v11 = _v[(3*seg_idx+1)][ax_idx];
        Scalar v12 = _v[(3*seg_idx)+3][ax_idx];
        Scalar T = segment_duration(seg_idx, ax_idx);
        const Scalar sc = 1.0 - 1e-5;

        Scalar v21 = _v[(3*(seg_idx-1)+1)][ax_idx];

        Scalar v0_max = std::max(std::max(fabs(v11), fabs(v21)), _v_max(ax_idx));
        Scalar v0_min = -v0_max;
        
        if (!_thrust_limit_flag){
            std::vector<Vector<2>> test_acc_vec = {Vector<2>(_a_min[ax_idx], _a_max[ax_idx]),Vector<2>(_a_max[ax_idx], _a_min[ax_idx])};
            for (Vector<2> a : test_acc_vec) {
                Scalar a11 = a(0);
                Scalar a12 = a(1);

                Scalar sigma_1 = sqrt(sc*(a11 - a12)*(- a12*sc*pow(T,2) + 2*v12*T + 2*p0 - 2*p12));
                Scalar t1_1 = (v12 - T*a12*sc)/(sc*(a11 - a12)) + (sigma_1 - v12 + T*a12*sc)/(sc*(a11 - a12));
                Scalar t1_2 = (v12 - T*a12*sc)/(sc*(a11 - a12)) - (v12 + sigma_1 - T*a12*sc)/(sc*(a11 - a12));
                Scalar t2_1 = - (v12 - T*a11*sc)/(sc*(a11 - a12)) - (sigma_1 - v12 + T*a12*sc)/(sc*(a11 - a12));
                Scalar t2_2 = (v12 + sigma_1 - T*a12*sc)/(sc*(a11 - a12)) - (v12 - T*a11*sc)/(sc*(a11 - a12));

                if (std::isfinite(t1_1) and std::isfinite(t2_1) and t1_1 > -PRECISION_PMM_VALUES and t2_1 > -PRECISION_PMM_VALUES) {
                    Scalar v0 = v12 - sigma_1 - T*a12*sc;
                    if (v0 > 0. and v0 < v0_max){
                        v0_max = v0;
                    } else if (v0 < 0. and v0 > v0_min){
                        v0_min = v0;
                    }
                }
                if (std::isfinite(t1_2) and std::isfinite(t2_2) and t1_2 > -PRECISION_PMM_VALUES and t2_2 > -PRECISION_PMM_VALUES) {
                    Scalar v0 = v12 + sigma_1 - T*a12*sc;
                    if (v0 > 0. and v0 < v0_max){
                        v0_max = v0;
                    } else if (v12 < 0. and v0 > v0_min){
                        v0_min = v0;
                    }
                }
            }
        }

        return std::tuple(v0_min, v0_max);

    }

    std::tuple<std::vector<Scalar>, std::vector<Vector<3>>, std::vector<Vector<3>>, std::vector<Vector<3>>> PMM_MG_Trajectory3D::get_sampled_trajectory(Scalar sampling_period){
        if (!_exists or n_segments_filled != n_segments){
            std::cout << "ERROR: Can not sample incomplete or unfeasible trajectory." << std::endl;
            std::cout << n_segments_filled << " <  " << n_segments << std::endl;
            return std::tuple(std::vector<Scalar>(), std::vector<Vector<3>>(), std::vector<Vector<3>>(), std::vector<Vector<3>>());
        }

        Scalar t_integrator_traj = 0.;
        Scalar t_traj = 0.;
        Scalar t_carry = 0.;
        std::vector<Scalar> t;
        std::vector<Vector<3>> p;
        std::vector<Vector<3>> v;
        std::vector<Vector<3>> a;

        for (int seg_idx=0; seg_idx<n_segments; seg_idx++){
            // Assuming equalized axis
            Scalar segment_duration = _dt[(3*seg_idx)][0] + _dt[(3*seg_idx)+1][0] + _dt[(3*seg_idx)+2][0];
            // incorporate t_carry 
            Scalar n_wp_per_segment = (segment_duration-t_carry)/sampling_period;
            int n_wp_to_sample = std::floor(n_wp_per_segment);

            // time left for the next sample in the following segment
            t_carry = (n_wp_per_segment - n_wp_to_sample) * sampling_period;

            // add initial state
            if (seg_idx == 0){
                n_wp_to_sample += 1;
            }

            Scalar t_integrator_seg = 0.;

            for (int wp_idx=0;wp_idx<n_wp_to_sample; wp_idx++){
                Vector<3> p_wp;
                Vector<3> v_wp;
                Vector<3> a_wp;
                // sample
                for (int ax_idx=0; ax_idx<3;ax_idx++){
                    if (t_integrator_seg < _dt[(3*seg_idx)][ax_idx]) {
                        p_wp[ax_idx] = _p[3*seg_idx][ax_idx] + _v[3*seg_idx][ax_idx] * t_integrator_seg + 0.5 * _a[(3*seg_idx)][ax_idx] * pow(t_integrator_seg,2);
                        v_wp[ax_idx] = _v[3*seg_idx][ax_idx] + _a[(3*seg_idx)][ax_idx] * t_integrator_seg;
                        a_wp[ax_idx] = _a[(3*seg_idx)][ax_idx];
                    } else if (t_integrator_seg < _dt[(3*seg_idx)][ax_idx] + _dt[(3*seg_idx)+1][ax_idx]) {
                        Scalar t1 = _dt[(3*seg_idx)][ax_idx];
                        p_wp[ax_idx] = _p[(3*seg_idx)+1][ax_idx] + _v[(3*seg_idx)+1][ax_idx] * (t_integrator_seg-t1) + 0.5 * _a[(3*seg_idx)+1][ax_idx] * pow(t_integrator_seg-t1,2);
                        v_wp[ax_idx] = _v[(3*seg_idx)+1][ax_idx] + _a[(3*seg_idx)+1][ax_idx] * (t_integrator_seg-t1);
                        a_wp[ax_idx] = _a[(3*seg_idx)+1][ax_idx];
                    } else {
                        Scalar t2 = _dt[(3*seg_idx)][ax_idx] + _dt[(3*seg_idx)+1][ax_idx];
                        p_wp[ax_idx] = _p[(3*seg_idx)+2][ax_idx] + _v[(3*seg_idx)+2][ax_idx] * (t_integrator_seg-t2) + 0.5 * _a[(3*seg_idx)+2][ax_idx] * pow(t_integrator_seg-t2,2);
                        v_wp[ax_idx] = _v[(3*seg_idx)+2][ax_idx] + _a[(3*seg_idx)+2][ax_idx] * (t_integrator_seg-t2);
                        a_wp[ax_idx] = _a[(3*seg_idx)+2][ax_idx];
                    }
                }
                t.push_back(t_integrator_traj);
                p.push_back(p_wp);
                v.push_back(v_wp);
                a.push_back(a_wp);

                // update
                t_integrator_seg += sampling_period;
                t_integrator_traj += sampling_period;
            }
            t_traj += segment_duration;
            t_integrator_traj = t_traj;
        }

        // add end configuration
        if (t_carry != 0.){
            t.push_back(t_integrator_traj);
            p.push_back(_p[3*(n_segments)]);
            v.push_back(_v[3*(n_segments)]);
            a.push_back(_a[3*(n_segments)-1]);
        }

        return std::tuple(t, p, v, a);
    }

    std::tuple<std::vector<Scalar>, std::vector<Scalar>, std::vector<Scalar>, std::vector<Scalar>> PMM_MG_Trajectory3D::get_sampled_single_axis_trajectory(Scalar sampling_period, int ax_idx){
        if (!_exists or n_segments_filled != n_segments){
            std::cout << "ERROR: Can not sample incomplete or unfeasible trajectory." << std::endl;
            return std::tuple(std::vector<Scalar>(), std::vector<Scalar>(), std::vector<Scalar>(), std::vector<Scalar>());
        }

        Scalar t_integrator_traj = 0.;
        Scalar t_carry = 0.;
        std::vector<Scalar> t;
        std::vector<Scalar> p;
        std::vector<Scalar> v;
        std::vector<Scalar> a;

        for (int seg_idx=0; seg_idx<n_segments; seg_idx++){
            // Assuming equalized axis
            Scalar segment_duration = _dt[(3*seg_idx)][ax_idx] + _dt[(3*seg_idx)+1][ax_idx] + _dt[(3*seg_idx)+2][ax_idx];
            // incorporate t_carry 
            Scalar n_wp_per_segment = (segment_duration-t_carry)/sampling_period;
            int n_wp_to_sample = std::floor(n_wp_per_segment);

            // time left for the next sample in the following segment
            t_carry = (n_wp_per_segment - n_wp_to_sample) * sampling_period;

            // add initial state
            if (seg_idx == 0){
                n_wp_to_sample += 1;
            }

            Scalar t_integrator_seg = 0.;

            for (int wp_idx=0;wp_idx<n_wp_to_sample; wp_idx++){
                Scalar p_wp;
                Scalar v_wp;
                Scalar a_wp;
                // sample
                if (t_integrator_seg < _dt[(3*seg_idx)][ax_idx]) {
                    p_wp = _p[3*seg_idx][ax_idx] + _v[3*seg_idx][ax_idx] * t_integrator_seg + 0.5 * _a[(3*seg_idx)][ax_idx] * pow(t_integrator_seg,2);
                    v_wp = _v[3*seg_idx][ax_idx] + _a[(3*seg_idx)][ax_idx] * t_integrator_seg;
                    a_wp = _a[(3*seg_idx)][ax_idx];
                } else if (t_integrator_seg < _dt[(3*seg_idx)][ax_idx] + _dt[(3*seg_idx)+1][ax_idx]) {
                    Scalar t1 = _dt[(3*seg_idx)][ax_idx];
                    p_wp = _p[(3*seg_idx)+1][ax_idx] + _v[(3*seg_idx)+1][ax_idx] * (t_integrator_seg-t1) + 0.5 * _a[(3*seg_idx)+1][ax_idx] * pow(t_integrator_seg-t1,2);
                    v_wp = _v[(3*seg_idx)+1][ax_idx] + _a[(3*seg_idx)+1][ax_idx] * (t_integrator_seg-t1);
                    a_wp = _a[(3*seg_idx)+1][ax_idx];
                } else {
                    Scalar t2 = _dt[(3*seg_idx)][ax_idx] + _dt[(3*seg_idx)+1][ax_idx];
                    p_wp = _p[(3*seg_idx)+2][ax_idx] + _v[(3*seg_idx)+2][ax_idx] * (t_integrator_seg-t2) + 0.5 * _a[(3*seg_idx)+2][ax_idx] * pow(t_integrator_seg-t2,2);
                    v_wp = _v[(3*seg_idx)+2][ax_idx] + _a[(3*seg_idx)+2][ax_idx] * (t_integrator_seg-t2);
                    a_wp = _a[(3*seg_idx)+2][ax_idx];
                }
                t.push_back(t_integrator_traj);
                p.push_back(p_wp);
                v.push_back(v_wp);
                a.push_back(a_wp);

                // update
                t_integrator_seg += sampling_period;
                t_integrator_traj += sampling_period;
            }
        }

        // add end configuration
        if (t_carry != 0.){
            t.push_back(t_integrator_traj);
            p.push_back(_p[3*(n_segments)][ax_idx]);
            v.push_back(_v[3*(n_segments)][ax_idx]);
            a.push_back(_a[3*(n_segments)-1][ax_idx]);
        }

        return std::tuple(t, p, v, a);
    }

    std::tuple<std::vector<Scalar>, std::vector<Vector<3>>> PMM_MG_Trajectory3D::get_sampled_trajectory_path(const Scalar sampling_distance){
        if (!_exists or n_segments_filled != n_segments){
            std::cout << "ERROR: Can not sample incomplete or unfeasible trajectory." << std::endl;
            return std::tuple(std::vector<Scalar>(), std::vector<Vector<3>>());
        }

        Scalar t_integrator_traj = 0.;
        Vector<3> ds_carry = Vector<3>(0.,0.,0.);
        int ds_carry_ax_idx = 0;
        std::vector<Scalar> t;
        std::vector<Vector<3>> path;

        path.push_back(_p[0]);
        t.push_back(0.0);

        Scalar step = sampling_distance;
        bool force_carry = false;

        int cnt = 0;

        for (int seg_idx=0; seg_idx<n_segments; seg_idx++){
            // Assuming equalized axis
            Scalar segment_duration = _dt[(3*seg_idx)][0] + _dt[(3*seg_idx)+1][0] + _dt[(3*seg_idx)+2][0];
            
            Scalar t_integrator_seg = 0.;

            bool terminate = false;
            Eigen::Matrix<int, 3, 1> subseg_idx;

            while (!terminate){
                Vector<3> t_tmp;
                Scalar min_ax_t = MAX_SCALAR;

                for (int ax_idx=0; ax_idx<3; ax_idx++){
                    Scalar a_ax;
                    Scalar v_ax;
                    bool tmp_terminate;

                    // select in which subsegment are we working
                    Scalar t1 = _dt[(3*seg_idx)][ax_idx];
                    if (t_integrator_seg < t1){
                        subseg_idx[ax_idx] = 0;
                        v_ax = _v[(3*seg_idx)+subseg_idx[ax_idx]][ax_idx] + _a[(3*seg_idx)+subseg_idx[ax_idx]][ax_idx] * t_integrator_seg;
                    } else {
                        subseg_idx[ax_idx] = 2;
                        v_ax = _v[(3*seg_idx)+subseg_idx[ax_idx]][ax_idx] + _a[(3*seg_idx)+subseg_idx[ax_idx]][ax_idx] * (t_integrator_seg-t1);
                    }
                    a_ax = _a[(3*seg_idx)+subseg_idx[ax_idx]][ax_idx];

                    if (fabs(a_ax)< PRECISION_PMM_VALUES){
                        continue;
                    }
                    
                    // solve ds == v*t + 0.5*a*t^2
                    if (force_carry){
                        step = ds_carry[ax_idx];
                    }
                    Scalar sigma_1 = pow(v_ax,2) + 2*a_ax*step;
                    Scalar sigma_2 = pow(v_ax,2) - 2*a_ax*step;
                    Scalar t_1 =  - (v_ax - sqrt(sigma_1))/a_ax;
                    Scalar t_2 =  - (v_ax + sqrt(sigma_1))/a_ax;
                    Scalar t_3 =  - (v_ax - sqrt(sigma_2))/a_ax;
                    Scalar t_4 =  - (v_ax + sqrt(sigma_2))/a_ax;

                    t_tmp[ax_idx] = MAX_SCALAR;
                    if (std::isfinite(t_1) and t_1 > 0. and t_1 < t_tmp[ax_idx]) {
                        t_tmp[ax_idx] = t_1;
                    } if (std::isfinite(t_2) and t_2 > 0. and t_2 < t_tmp[ax_idx]){
                        t_tmp[ax_idx] = t_2;
                    } if (std::isfinite(t_3) and t_3 > 0. and t_3 < t_tmp[ax_idx]){
                        t_tmp[ax_idx] = t_3;
                    } if (std::isfinite(t_4) and t_4 > 0. and t_4 < t_tmp[ax_idx]){
                        t_tmp[ax_idx] = t_4;
                    }

                    // check if time did not overflow the current segment
                    if (t_integrator_seg + t_tmp[ax_idx] > segment_duration){
                        Scalar t_left = segment_duration-t_integrator_seg;
                        t_tmp[ax_idx] = t_left;
                        Scalar ds_seg = abs(v_ax*t_left + 0.5 * a_ax * pow(t_left,2));
                        ds_carry[ax_idx] = sampling_distance - ds_seg;
                        tmp_terminate = true;
                    } else {
                        tmp_terminate = false;
                    }

                    // check if acc did not change throughout the segment
                    if ((subseg_idx[ax_idx] == 0) and (t_integrator_seg+t_tmp[ax_idx] > t1)){
                        Scalar t_subseg1 = t1-t_integrator_seg;
                        Scalar ds_subseg1 = abs(v_ax*t_subseg1 + 0.5 * a_ax * pow(t_subseg1,2));
                        Scalar ds_left = step - ds_subseg1;

                        // compute t for subseg2
                        subseg_idx[ax_idx] = 2;
                        a_ax = _a[(3*seg_idx)+subseg_idx[ax_idx]][ax_idx];
                        v_ax = _v[(3*seg_idx)+subseg_idx[ax_idx]][ax_idx];

                        sigma_1 = pow(v_ax,2) + 2*a_ax*ds_left;
                        sigma_2 = pow(v_ax,2) - 2*a_ax*ds_left;
                        t_1 =  - (v_ax - sqrt(sigma_1))/a_ax;
                        t_2 =  - (v_ax + sqrt(sigma_1))/a_ax;
                        t_3 =  - (v_ax - sqrt(sigma_2))/a_ax;
                        t_4 =  - (v_ax + sqrt(sigma_2))/a_ax;

                        t_tmp[ax_idx] = MAX_SCALAR;
                        if (std::isfinite(t_1) and t_1 > PRECISION_PMM_VALUES and t_1 < t_tmp[ax_idx]) {
                            t_tmp[ax_idx] = t_1;
                        } if (std::isfinite(t_2) and t_2 > PRECISION_PMM_VALUES and t_2 < t_tmp[ax_idx]){
                            t_tmp[ax_idx] = t_2;
                        } if (std::isfinite(t_3) and t_3 > PRECISION_PMM_VALUES and t_3 < t_tmp[ax_idx]){
                            t_tmp[ax_idx] = t_3;
                        } if (std::isfinite(t_4) and t_4 > PRECISION_PMM_VALUES and t_4 < t_tmp[ax_idx]){
                            t_tmp[ax_idx] = t_4;
                        }

                        // check again if time did not overflow the current segment
                        if (t_integrator_seg + t_tmp[ax_idx] > segment_duration){
                            Scalar t_left = segment_duration-t_integrator_seg;
                            t_tmp[ax_idx] = t_left;
                            Scalar ds_seg = abs(v_ax*t_left + 0.5 * a_ax * pow(t_left,2));
                            ds_carry[ax_idx] = sampling_distance - ds_seg;
                            tmp_terminate = true;
                        } else {
                            tmp_terminate = false;
                        }
                    }

                    if (t_tmp[ax_idx] < min_ax_t){
                        min_ax_t = t_tmp[ax_idx];
                        terminate = tmp_terminate;
                    }
                }

                if (fabs(min_ax_t - MAX_SCALAR) < PRECISION_PMM_VALUES){
                    // All three axis were zero trajectories
                    terminate = true;
                    force_carry = true;
                    continue;
                }

                if (terminate) {
                    force_carry = true;
                    continue;
                }
                
                cnt += 1 ;
                t_integrator_seg += min_ax_t;
                t.push_back(t_integrator_traj+t_integrator_seg);

                Vector<3> p_wp;
                for (int ax_idx=0; ax_idx<3;ax_idx++){
                    if (t_integrator_seg < _dt[(3*seg_idx)][ax_idx]) {
                        p_wp[ax_idx] = _p[3*seg_idx][ax_idx] + _v[3*seg_idx][ax_idx] * t_integrator_seg + 0.5 * _a[(3*seg_idx)][ax_idx] * pow(t_integrator_seg,2);
                    } else if (t_integrator_seg < _dt[(3*seg_idx)][ax_idx] + _dt[(3*seg_idx)+1][ax_idx]) {
                        Scalar t1 = _dt[(3*seg_idx)][ax_idx];
                        p_wp[ax_idx] = _p[(3*seg_idx)+1][ax_idx] + _v[(3*seg_idx)+1][ax_idx] * (t_integrator_seg-t1) + 0.5 * _a[(3*seg_idx)+1][ax_idx] * pow(t_integrator_seg-t1,2);
                    } else {
                        Scalar t2 = _dt[(3*seg_idx)][ax_idx] + _dt[(3*seg_idx)+1][ax_idx];
                        p_wp[ax_idx] = _p[(3*seg_idx)+2][ax_idx] + _v[(3*seg_idx)+2][ax_idx] * (t_integrator_seg-t2) + 0.5 * _a[(3*seg_idx)+2][ax_idx] * pow(t_integrator_seg-t2,2);
                    }
                }
                path.push_back(p_wp);

                if (force_carry){
                    // reset step
                    step = sampling_distance;
                    force_carry = false;
                }
            }

            t_integrator_traj += segment_duration;
        }

        // add end configuration
        if (fabs(ds_carry[0]) > PRECISION_PMM_VALUES or fabs(ds_carry[1]) > PRECISION_PMM_VALUES or fabs(ds_carry[2]) > PRECISION_PMM_VALUES){
            path.push_back(_p[_p.size()-1]);
            t.push_back(t_integrator_traj);
        }
        return std::tuple(t, path);
    }

    void PMM_MG_Trajectory3D::sample_and_export_trajectory(Scalar sampling_period, std::string output_file){
        // Get sampled trajectory
        std::vector<Scalar> t_s;
        std::vector<Vector<3>> p_s;
        std::vector<Vector<3>> v_s;
        std::vector<Vector<3>> a_s;

        std::tie(t_s, p_s, v_s, a_s) = this->get_sampled_trajectory(sampling_period);

        // export sampled trajectory
        std::ofstream outFile(output_file);
        for (int i=0;i<t_s.size();i++) {
            outFile << t_s[i] << "," << p_s[i](0) << "," << p_s[i](1) << "," << p_s[i](2) << "," << v_s[i](0) << "," << v_s[i](1) << "," << v_s[i](2) << "," << a_s[i](0) << "," << a_s[i](1) << "," << a_s[i](2) << std::endl;
        }
        
        outFile.close();
    }

    void PMM_MG_Trajectory3D::sample_and_export_trajectory_column(Scalar sampling_period, std::string output_file){
        // Get sampled trajectory
        std::vector<Scalar> t_s;
        std::vector<Vector<3>> p_s;
        std::vector<Vector<3>> v_s;
        std::vector<Vector<3>> a_s;

        std::tie(t_s, p_s, v_s, a_s) = this->get_sampled_trajectory(sampling_period);

        // export sampled trajectory
        std::ofstream outFile(output_file);
        writeVectorToFile(t_s, &outFile);
        writeVectorToFile(p_s, &outFile);
        writeVectorToFile(v_s, &outFile);
        writeVectorToFile(a_s, &outFile);
        outFile.close();
    }

    void PMM_MG_Trajectory3D::sample_and_export_unsynchronized_trajectory(Scalar sampling_period, std::string output_file){
        // Get sampled trajectory
        std::vector<Scalar> t_x;
        std::vector<Scalar> p_x;
        std::vector<Scalar> v_x;
        std::vector<Scalar> a_x;
        std::tie(t_x, p_x, v_x, a_x) = this->get_sampled_single_axis_trajectory(sampling_period, 0);

        std::vector<Scalar> t_y;
        std::vector<Scalar> p_y;
        std::vector<Scalar> v_y;
        std::vector<Scalar> a_y;
        std::tie(t_y, p_y, v_y, a_y) = this->get_sampled_single_axis_trajectory(sampling_period, 1);

        std::vector<Scalar> t_z;
        std::vector<Scalar> p_z;
        std::vector<Scalar> v_z;
        std::vector<Scalar> a_z;
        std::tie(t_z, p_z, v_z, a_z) = this->get_sampled_single_axis_trajectory(sampling_period, 2);

        // export sampled trajectory
        std::ofstream outFile(output_file);
        writeVectorToFile(t_x, &outFile);
        writeVectorToFile(t_x, &outFile);
        writeVectorToFile(t_y, &outFile);
        writeVectorToFile(t_z, &outFile);
        writeVectorToFile(p_x, &outFile);
        writeVectorToFile(p_y, &outFile);
        writeVectorToFile(p_z, &outFile);
        writeVectorToFile(v_x, &outFile);
        writeVectorToFile(v_y, &outFile);
        writeVectorToFile(v_z, &outFile);
        writeVectorToFile(a_x, &outFile);
        writeVectorToFile(a_y, &outFile);
        writeVectorToFile(a_z, &outFile);
        outFile.close();
    }

    void PMM_MG_Trajectory3D::sample_and_export_trajectory_path(Scalar sampling_distance, std::string output_file){
        // Get sampled trajectory
        std::vector<Scalar> t_s;
        std::vector<Vector<3>> p_s;

        std::tie(t_s, p_s) = this->get_sampled_trajectory_path(sampling_distance);

        // export sampled trajectory
        std::ofstream outFile(output_file);
        writeVectorToFile(t_s, &outFile);
        writeVectorToFile(p_s, &outFile);
        outFile.close();
    }
    
    Scalar PMM_MG_Trajectory3D::duration(){
        Scalar t_integrator = 0.;
        for (int seg_idx=0;seg_idx<n_segments;seg_idx++){
            t_integrator += segment_duration(seg_idx, 0);
        }
        return t_integrator;
    }

    Scalar PMM_MG_Trajectory3D::segment_duration(int segment_idx, int ax_idx){
        // assuming synchronized axis
        return _dt[(3*segment_idx)][ax_idx] + _dt[(3*segment_idx)+1][ax_idx] + _dt[(3*segment_idx)+2][ax_idx];
    }

    std::vector<Vector<3>> PMM_MG_Trajectory3D::get_current_velocities(){
        std::vector<Vector<3>> v_out;
        if (n_segments == 0){
            return v_out;
        }
        for (int i=0; i<n_segments+1;i++){
            v_out.push_back(_v[i*3]);
        }
        return v_out;
    }

    std::vector<Vector<3>> PMM_MG_Trajectory3D::get_current_positions(){
        std::vector<Vector<3>> p_out;
        if (n_segments == 0){
            return p_out;
        }
        for (int i=0; i<n_segments+1;i++){
            p_out.push_back(_p[i*3]);
        }
        return p_out;
    }

    void PMM_MG_Trajectory3D::exportTrajectoryToCSV(std::string file_name, bool optim_data){
        std::ofstream outFile(file_name);

        // fill the trajectory data
        outFile << this->n_segments_filled << std::endl;

        // dt
        writeVectorToFile(_p, &outFile);
        writeVectorToFile(_v, &outFile);
        writeVectorToFile(_a, &outFile);
        writeVectorToFile(_dtdv, &outFile);
        writeVectorToFile(_dt, &outFile);
        if(optim_data){
            writeVectorToFile(_used_grad, &outFile);
        }
        outFile.close();
    }

    void PMM_MG_Trajectory3D::exportTrajectoryToCSV(std::ofstream *outStream, bool optim_data){
        // fill the trajectory data
        *outStream << n_segments_filled << std::endl;
        // dt
        writeVectorToFile(_p, outStream);
        writeVectorToFile(_v, outStream);
        writeVectorToFile(_a, outStream);
        writeVectorToFile(_dtdv, outStream);
        writeVectorToFile(_dt, outStream);
        if(optim_data){
            writeVectorToFile(_used_grad, outStream);
        }
    }

    void PMM_MG_Trajectory3D::writeVectorToFile(std::vector<Scalar> v, std::ofstream *outFile){
        for (auto it = v.begin(); it != v.end(); it++){
            if (it != v.begin()) {
                *outFile << ",";
            }
            *outFile << *it;
        }
        *outFile << std::endl;
    }

    void PMM_MG_Trajectory3D::writeVectorToFile(Vector<3> v, std::ofstream *outFile){
        for (int i = 0; i < 3; i++){
            if (i != 0) {
                *outFile << ",";
            }
            *outFile << v[i];
        }
        *outFile << std::endl;
    }

    void PMM_MG_Trajectory3D::writeVectorToFile(Vector<4> v, std::ofstream *outFile){
        for (int i = 0; i < 4; i++){
            if (i != 0) {
                *outFile << ",";
            }
            *outFile << v[i];
        }
        *outFile << std::endl;
    }

    void PMM_MG_Trajectory3D::writeVectorToFile(std::vector<Vector<3>> v, std::ofstream *outFile){
        for (int i = 0; i < 3; i++){
            for (int j=0; j<v.size(); j++){
                if (j != 0) {
                    *outFile << ",";
                }
                *outFile << v[j][i];
            }
            *outFile << std::endl;
        }
    }

    void PMM_MG_Trajectory3D:: writeVectorToFile(std::vector<Vector<4>> v, std::ofstream *outFile){
        for (int i = 0; i < 4; i++){
            for (int j=0; j<v.size(); j++){
                if (j != 0) {
                    *outFile << ",";
                }
                *outFile << v[j][i];
            }
            *outFile << std::endl;
        }
    }

    PMM_MG_Trajectory3D PMM_MG_Trajectory3D::operator=(const PMM_MG_Trajectory3D &tr){
        this->copy_trajectory(tr);
        return *this;
    }

    void PMM_MG_Trajectory3D::copy_trajectory(const PMM_MG_Trajectory3D &tr){
        this->n_segments = tr.n_segments;
        _dt.resize(3 * n_segments);
        _dtdv.resize(2 * n_segments);
        _p.resize(3 *n_segments +1);
        _v.resize(3 *n_segments +1);
        _a.resize(3 * n_segments);
        _traj_type_vec.resize(n_segments);
        _sol_idx.resize(n_segments);
        _update_vel_flag.resize(n_segments);
        _ax_min_traj_duration.resize(n_segments);
        _used_grad.resize(n_segments);

        this->_exists = tr._exists;
        this->_thrust_limit_flag = tr._thrust_limit_flag;
        this->n_segments_filled = tr.n_segments_filled;
        this->_a_max = tr._a_max;
        this->_a_min = tr._a_min;
        this->_v_max = tr._v_max;

        this->_p = tr._p;
        this->_v = tr._v;
        this->_a = tr._a;
        this->_dt = tr._dt;
        this->_dtdv = tr._dtdv;
        this->_traj_type_vec = tr._traj_type_vec;
        this->_ax_min_traj_duration = tr._ax_min_traj_duration;

        this->_max_acc_norm = tr._max_acc_norm;
        this->_max_vel_norm = tr._max_vel_norm;
        this->_TD_max_iter = tr._TD_max_iter;
        this->_TD_alpha_factor = tr._TD_alpha_factor;
        this->_TD_norm_precision = tr._TD_norm_precision;
    }

    std::tuple<Vector<3>, Vector<3>> compute_per_axis_acc_limits_from_limited_thrust(const Scalar max_acc_norm){
        // Compute per axis limits from limited thrust acceleration
        // B = T + GVEC , |T|=max_acc_norm
        // initially B equal per axis with b_x=b_y=b_z -> |B-GVEC|^2 = |T|^2
        // -> 3*bs_x^2 + 2*g*a_x + g^2 - |T|^2 = 0 --> roots are the possible acc

        const Scalar max_acc_norm_pow2 = max_acc_norm * max_acc_norm;
        const Scalar a_equal_acc = 3;
        const Scalar b_equal_acc = 2 * G;
        const Scalar c_equal_acc = G * G - max_acc_norm_pow2;
        const Scalar equal_acc_1 = (-b_equal_acc + sqrt(b_equal_acc * b_equal_acc - 4 * a_equal_acc * c_equal_acc)) / (2 * a_equal_acc);

        Vector<3> per_axis_acc_vec;
        if ((Vector<3>::Constant(-equal_acc_1) - GVEC).norm() - max_acc_norm < PRECISION_PMM_VALUES) {
            per_axis_acc_vec = Vector<3>::Constant(equal_acc_1);
        } else {
            const Scalar equal_acc_2 =
            (-b_equal_acc - sqrt(b_equal_acc * b_equal_acc - 4 * a_equal_acc * c_equal_acc)) / (2 * a_equal_acc);
            per_axis_acc_vec = Vector<3>::Constant(equal_acc_2);
            std::cout << "WARNING: This does no usually happen! " << (Vector<3>::Constant(-equal_acc_1) - GVEC).norm() << std::endl;
        }

        Vector<3> max_per_axis_acc_vec = per_axis_acc_vec;
        Vector<3> min_per_axis_acc_vec = -per_axis_acc_vec + 2*GVEC;

        return std::tuple(min_per_axis_acc_vec, max_per_axis_acc_vec);
    }

    std::vector<Vector<3>> compute_initial_velocities(std::vector<Vector<3>> waypoints, Vector<3> start_velocity,
                     Vector <3> end_velocity, const Scalar max_per_axis_acc, const Scalar max_per_axis_vel, const Scalar const_vel_norm){

        int n_waypoints = waypoints.size();

        std::vector<Vector<3>> initial_velocities;
        initial_velocities.push_back(start_velocity);
        
        // compute normalized headings vectors towards the next waypoint and path segment lengths
        std::vector<Vector<3>> heading_vectors;
        std::vector<Scalar> path_segment_len;

        for (int i = 0; i < n_waypoints-1; i++){
            Vector<3> tmp_heading = waypoints[i+1] - waypoints[i];
            path_segment_len.push_back(tmp_heading.norm());
            heading_vectors.push_back(tmp_heading);

            if (tmp_heading.norm() < PRECISION_PMM_VALUES){
                std::cout << "WARNING: Trajectory contains two identical waypoints." << std::endl;
            }
        }

        for (int i = 1; i<n_waypoints-1; i++){
            Vector<3> h_1 = heading_vectors[i-1].normalized();
            Vector<3> h_2 = heading_vectors[i].normalized();
            Scalar l_1 = path_segment_len[i-1];
            Scalar l_2 = path_segment_len[i];

            if (l_1 < PRECISION_PMM_VALUES or l_2 < PRECISION_PMM_VALUES){
                initial_velocities.push_back(Vector<3>(0,0,0));
                continue;
            }

            // angle between vectors
            Scalar theta = std::acos((h_1.dot(h_2))/(h_1.norm()*h_2.norm()));
            Vector<3> v = h_1.cross(h_2);

            // compute heading - resulting angle as ratio of segment lengths
            Scalar theta_new = (0.5 + 0.6*((l_1/(l_1 + l_2))-0.5)) * theta;

            // compute resulting velocity norm according to the selected method
            Scalar v_norm_new;
            Scalar angl_sc_factor;
            if (const_vel_norm >= 0.0) {
                // default vel norm for all waypoints
                v_norm_new = const_vel_norm;

                if (theta > (3*M_PI)/4){
                    v_norm_new = 0;
                }
            } else {
                // consider worst case scenario, i.e. zero velocity at the previous and next wpt
                Vector<3> hv_1 = heading_vectors[i-1].cwiseAbs();
                Vector<3> hv_2 = heading_vectors[i].cwiseAbs();

                Scalar a_max = fabs(max_per_axis_acc);

                Vector<3> v_tmp;

                for (int ax=0;ax<3;ax++){
                    Scalar v_1 = sqrt(2*a_max*hv_1[ax]);
                    Scalar v_2 = sqrt(2*a_max*hv_2[ax]);
                    // v_tmp[ax] = std::min(v_1, v_2);
                    v_tmp[ax] = (v_1 + v_2) / 2;
                }

                // utilize prior knowledge about the sharpness of the angle
                if (theta < M_PI_2){
                    angl_sc_factor = 1;
                } else if (theta < ((3*M_PI)/4)){
                    angl_sc_factor = 0.8;
                } else {
                    angl_sc_factor = 0;
                }
                
                v_norm_new = 0.5 * angl_sc_factor * v_tmp.norm();
            }

            // rotate vector using Rodriguez formula
            Vector<3> v_new_tmp = h_1 * std::cos(theta_new) + v.cross(h_1)* std::sin(theta_new) + v * (v.dot(h_1))*(1-cos(theta_new));

            Vector<3> v_new = v_norm_new*v_new_tmp;

            for (int i=0;i<3;i++) {
                if (std::isnan(v_new(i))) {
                    v_new(i) = 0;
                } else if (fabs(v_new(i)) > max_per_axis_vel) {
                    v_new(i) = std::copysign( max_per_axis_vel, v_new(i) );
                } 
            }
            
            initial_velocities.push_back(v_new);
        }

        initial_velocities.push_back(end_velocity);
        return initial_velocities;

    }
} // namespace pmm