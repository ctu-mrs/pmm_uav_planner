/**
 * @file pmm_trajectory.hpp
 * @author Robert Penicka (penicrob@fel.cvut.cz), Krystof Teissing (teisskry@gmail.com), Matej Novosad (novosma2@fel.cvut.cz)
 * @version 0.1
 * @date 2024-09-20
 * 
 */


#pragma once
#include <cmath>
#include "common.hpp"

#define PRECISION_PMM_VALUES (1.0e-8)
#define PRECISION_DTDV_DEF (1e-4)
#define PRECISION_TRANS3D (1E-4)
#define MIN_ACC_REQ (1.0)
#define INITIAL_MIN_ACC (0.5)

namespace pmm {

// Enum to store trajectory type
enum trajectory_type {
    MIN_TRAJ,
    SYNC_TRAJ,
    ZERO_TRAJ
};

static constexpr Scalar MAX_SCALAR = std::numeric_limits<Scalar>::max();

class PMMTrajectory {
 public:  
  PMMTrajectory();
  PMMTrajectory(const PMMTrajectory &in);
  /*
  this solves the equations
  v1 == vs + t1*(a1);
  ve == v1 + t2*(a2);
  p1 == ps + t1*vs + 1/2*t1^2*(a1);
  pe == p1 + t2*v1 + 1/2*t2^2*(a2);
  */
  PMMTrajectory(const Scalar ps, const Scalar vs, const Scalar pe,
                const Scalar ve, const Scalar a1_in, const Scalar a2_in,
                const Scalar max_v, const int i, 
                const double desired_time = 0,
                const bool keep_acc_sign = false,
                const bool calc_gradient = false,
                const bool check_result = false);

  /**
   * @brief Computes synchronization trajectory using acceleration scaling for a given trajectory duration
   * 
   * @param total_time Synchronization time
   * @param a_max1 acceleration bounds
   * @param a_max2 acceleration bounds
   * @param max_v maximal velocity norm
   * @return Scalar 
   */
  Scalar computeSyncTrajectory(const Scalar total_time, const Scalar a_max1, const Scalar a_max2, const Scalar max_v);

  static Scalar minRequiredAcc(const Scalar ps, const Scalar vs,
                               const Scalar pe, const Scalar ve);

  Scalar time() const { return t_.sum(); };
  Scalar max_end_velocity_abs() const;
  Vector<3> state_in_time(const Scalar time_in_tr) const;
  Scalar acc_in_time(const Scalar time_in_tr) const;
  Scalar vel_in_time(const Scalar time_in_tr) const;
  std::tuple<PMMTrajectory, PMMTrajectory> split_in_time(
    const Scalar time_in_tr);
  void copy_trajectory(const PMMTrajectory &in);
  void save_to_file(std::string filename);

  bool exists_;
  int trajectory_type_ = MIN_TRAJ;
  Vector<3> t_;
  Vector<4> p_;
  Vector<3> v_;
  Vector<2> a_;
  Scalar dt_dvs_;
  Scalar dt_dve_;
  int sol_idx_;
  int i_;

  Scalar v_max;
  Scalar v_min;
};

std::ostream &operator<<(std::ostream &o, const PMMTrajectory &f);
}  // namespace pmm