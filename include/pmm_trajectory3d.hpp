/**
 * @file pmm_trajectory3d.hpp
 * @author Matej Novosad (novosma2@fel.cvut.cz), Robert Penicka (penicrob@fel.cvut.cz), Krystof Teissing (teisskry@gmail.com)
 * @version 0.1
 * @date 2024-09-20
 * 
 */

#pragma once
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>

#include "common.hpp"
#include "pmm_trajectory.hpp"

#define CONVERGENCE_PRECISION (1E-3)
#define TIME_PRECISION (1E-8)

namespace pmm {

class PointMassTrajectory3D {
 public:
  PointMassTrajectory3D();

  /*
  version of per axis acceleration limits with sync segments
  */
  PointMassTrajectory3D(const QuadState &from, const QuadState &to,
                        const Vector<3> max_acc, const Vector<3> max_vel,
                        const bool equalize_time,
                        const bool calc_gradient);
  PointMassTrajectory3D(const QuadState &from, const QuadState &to,
                        const Vector<3> max_acc1, const Vector<3> max_acc2,
                        const Vector<3> max_vel,
                        const bool equalize_time,
                        const bool calc_gradient);

  /*
  version that limit the thrust by norm using iterative scaling
  */
  PointMassTrajectory3D(const QuadState &from, const QuadState &to,
                        const Scalar max_acc_norm, const Scalar max_vel_norm, const int max_iter,
                        const double precision_acc_limit, const bool debug = false);

  /*
  version that limit the thrust by norm using iterative scaling with drag
  */
  PointMassTrajectory3D(const QuadState &from, const QuadState &to,
                        const Scalar max_acc_norm, const Scalar max_vel_norm, const bool drag, const int max_iter,
                        const double precision_acc_limit, const bool debug = false);

  bool exists() const { return x_.exists_ && y_.exists_ && z_.exists_; };
  Scalar time() const { return std::max({x_.time(), y_.time(), z_.time()}); };
  Scalar time_min() const {
    return std::min({x_.time(), y_.time(), z_.time()});
  }
  Vector<3> acc_in_time(const Scalar time_in_tr) const;
  Vector<3> vel_in_time(const Scalar time_in_tr) const;

  QuadState state_in_time(const Scalar time_in_tr) const;
  QuadState get_start_state() const;
  QuadState get_end_state() const;
  
  Vector<3> start_acc() const;
  Vector<3> end_acc() const;
  Vector<3> switch_times() const;
  Vector<3> zero_time() const;
  Vector<3> max_velocity() const;
  Vector<3> drag(Vector<3> acc, const Scalar time_in_tr) const;
  Vector<3> dt_dvs_() const {
    return Vector<3>(x_.dt_dvs_, y_.dt_dvs_, z_.dt_dvs_);
  };
  Vector<3> dt_dve_() const {
    return Vector<3>(x_.dt_dve_, y_.dt_dve_, z_.dt_dve_);
  };

  Vector<3> max_end_velocity_abs() const {
    return Vector<3>(x_.max_end_velocity_abs(), y_.max_end_velocity_abs(),
                     z_.max_end_velocity_abs());
  };
  void set_axis_trajectory(const int i, const PMMTrajectory tr);
  PMMTrajectory &get_axis_trajectory(const int i);
  Scalar get_axis_switch_time(const int i) const;

  PointMassTrajectory3D operator=(const PointMassTrajectory3D &tr);
  void copy_trajectory(const PointMassTrajectory3D &tr);

  PMMTrajectory x_;
  PMMTrajectory y_;
  PMMTrajectory z_;

  Vector<3> min_axis_trajectory_duration_;

 private:
};

std::ostream &operator<<(std::ostream &o, const PointMassTrajectory3D &f);
}  // namespace pmm