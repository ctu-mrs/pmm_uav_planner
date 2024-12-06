/**
 * @file pmm_trajectory.cpp
 * @author Robert Penicka (penicrob@fel.cvut.cz), Krystof Teissing (teisskry@gmail.com), Matej Novosad (novosma2@fel.cvut.cz)
 * @version 0.1
 * @date 2024-09-20
 * 
 */

#include "pmm_trajectory.hpp"

#include <cmath>
#include <fstream>
#include <iostream>


namespace pmm {

PMMTrajectory::PMMTrajectory()
  : exists_(false),
    t_(Vector<3>(0, 0, 0)),
    p_(Vector<4>(0, 0, 0, 0)),
    v_(Vector<3>(0, 0, 0)),
    a_(Vector<2>(0, 0)),
    i_(0) {}

PMMTrajectory::PMMTrajectory(const PMMTrajectory &in) { copy_trajectory(in); }

PMMTrajectory::PMMTrajectory(const Scalar ps, const Scalar vs, const Scalar pe,
                             const Scalar ve, const Scalar a1_in,
                             const Scalar a2_in, const Scalar max_v, const int i,
                             const double desired_time,
                             const bool keep_acc_sign, const bool calc_gradient,
                             const bool check_result) {

  v_max = max_v;
  v_min = -max_v;

  i_ = i;
  p_(0) = ps;
  p_(3) = pe;
  v_(0) = vs;
  v_(2) = ve;

  // already there
  if (fabs(pe - ps) < PRECISION_PMM_VALUES &&
      fabs(ve - vs) < PRECISION_PMM_VALUES && fabs(ve) < PRECISION_PMM_VALUES) {
    a_(0) = 0;
    a_(1) = 0;
    p_(0) = p_(1) = p_(2) = p_(3) = ps;
    v_(0) = v_(1) = v_(2) = vs;
    trajectory_type_ = ZERO_TRAJ;
    t_ = {0, 0, 0};
    // std::cout << "already there" << std::endl;
    exists_ = true;
    return;
  }

  // can not move without acc
  if (fabs(a1_in) <= PRECISION_PMM_VALUES &&
      fabs(a2_in) <= PRECISION_PMM_VALUES) {
    exists_ = false;
    return;
  }

  const Scalar pow_ve2 = ve * ve;
  const Scalar pow_vs2 = vs * vs;
  const Scalar pow_vm2 = v_max * v_max;

  Scalar t1 = MAX_SCALAR;
  Scalar t2 = MAX_SCALAR;
  Scalar t3 = MAX_SCALAR;
  Scalar dt_dvs = MAX_SCALAR;
  Scalar dt_dve = MAX_SCALAR;

  Scalar used_acc1 = a1_in;
  Scalar used_acc2 = a2_in;
  std::vector<Vector<2>> test_acc_vec = {Vector<2>(a1_in, a2_in),
                                         Vector<2>(a2_in, a1_in)};
  if (keep_acc_sign) {
    test_acc_vec = {Vector<2>(a1_in, a2_in)};
  }
  for (Vector<2> a : test_acc_vec) {
    const Scalar a1 = a(0);
    const Scalar a2 = a(1);
    const Scalar pow_a2_2 = a2 * a2;
    const Scalar pow_a1_2 = a1 * a1;
    const Scalar pow_a1_3 = pow_a1_2 * a1;


    //  t1 and t2 defined for:
    //      ((vs == ve) and (((a1 - a2)*(a1*ve*ve - a2*ve*ve - 2*a1*a2*pe + 2*a1*a2*ps)) >= 0))
    //  or  ((vs != ve) and ((-(a1 - a2)*(a1*ve*ve - a2*vs*vs - 2*a1*a2*pe + 2*a1*a2*ps)) <= 0))
    //
    //  dtdvs and dtdve defined for
    //     ((vs == ve) and (((a1 - a2)*(a1*ve*ve - a2*ve*ve - 2*a1*a2*pe + 2*a1*a2*ps)) > 0))
    //  or   ((vs != ve) and ((-(a1 - a2)*(a1*ve*ve - a2*vs*vs - 2*a1*a2*pe + 2*a1*a2*ps)) < 0))
 
    if (v_max < MAX_SCALAR) {

      // case 1 bang-singular-bang
      const Scalar t1_3 = (v_max - vs) / a1;
      const Scalar t3_3 = (ve - v_max) / a2;

      const Scalar t2_3 = (a1 * (2 * a2 * (pe - ps) + pow_vm2 - pow_ve2) + a2 * (pow_vs2 - pow_vm2)) / 2 / a1 / a2 / v_max;

      // case 2 bang-singular-bang
      const Scalar t1_4 = (v_min - vs) / a1;
      const Scalar t2_4 = (a1 * (2 * a2 * (pe - ps) + pow_vm2 - pow_ve2) + a2 * (pow_vs2 - pow_vm2)) / 2 / a1 / a2 / v_min;
      const Scalar t3_4 = (ve - v_min) / a2;

      if (std::isfinite(t1_3) and std::isfinite(t2_3) and std::isfinite(t3_3) and
          t1_3 > -PRECISION_PMM_VALUES and t2_3 > -PRECISION_PMM_VALUES and t3_3 > -PRECISION_PMM_VALUES and
          fabs(desired_time - (t1_3 + t2_3 + t3_3)) < fabs(desired_time - (t1 + t2 + t3))) {
        // clip time values
        // std ::cout << "t1 " << t1_3 << " t2 " << t2_3 << " t3 " << t3_3 << std::endl;
        t1 = std::max(t1_3, 0.0);
        t2 = std::max(t2_3, 0.0);
        t3 = std::max(t3_3, 0.0);

        sol_idx_ = 2;

        if (calc_gradient) {
          const Scalar d_t_dvs = (-v_max + vs) / a1 / v_max;
          const Scalar d_t_dve = (v_max - ve) / a2 / v_max;;

          dt_dvs = d_t_dvs;
          dt_dve = d_t_dve;
        }

        used_acc1 = a1;
        used_acc2 = a2;
      }
      else if (std::isfinite(t1_4) and std::isfinite(t2_4) and std::isfinite(t3_4) and
          t1_4 > -PRECISION_PMM_VALUES and t2_4 > -PRECISION_PMM_VALUES and t3_4 > -PRECISION_PMM_VALUES and
          fabs(desired_time - (t1_4 + t2_4 + t3_4)) < fabs(desired_time - (t1 + t2 + t3))) {
        // clip time values
        t1 = std::max(t1_4, 0.0);
        t2 = std::max(t2_4, 0.0);
        t3 = std::max(t3_4, 0.0);
 
        sol_idx_ = 3;

        if (calc_gradient) {
          const Scalar d_t_dvs = (-v_min + vs) / a1 / v_min;
          const Scalar d_t_dve = (v_min - ve) / a2 / v_min;;

          dt_dvs = d_t_dvs;
          dt_dve = d_t_dve;
        }

        used_acc1 = a1;
        used_acc2 = a2;
      }
    }

    if ((vs == ve) and (((a1 - a2)*(a1*ve*ve - a2*ve*ve - 2*a1*a2*pe + 2*a1*a2*ps)) < -PRECISION_DTDV_DEF)){
      // std::cout << "ERROR: PMM equations solution does not exist" << ((a1 - a2)*(a1*ve*ve - a2*ve*ve - 2*a1*a2*pe + 2*a1*a2*ps)) << std::endl;
      continue;
    } else if ((vs != ve) and (((-a1 + a2)*(-a1*ve*ve + a2*vs*vs + 2*a1*a2*pe - 2*a1*a2*ps)) < -PRECISION_DTDV_DEF)){  
      // std::cout << "ERROR: PMM equations solution does not exist, discriminant: " << ((-a1 + a2)*(-a1*ve*ve + a2*vs*vs + 2*a1*a2*pe - 2*a1*a2*ps)) << std::endl;
      continue;
    } 

    Scalar tst1, tst2;
    if (fabs((-a1 + a2) * (2 * a1 * a2 * (pe - ps) - a1 * pow_ve2 + a2 * pow_vs2)) < PRECISION_DTDV_DEF) {
      tst1 = 0;
    } else {
      tst1 = sqrt(
        (-a1 + a2) * (2 * a1 * a2 * (pe - ps) - a1 * pow_ve2 + a2 * pow_vs2));
    }
    if ((a1 - a2) * (a1 * (-2 * a2 * pe + 2 * a2 * ps + pow_ve2) - a2 * pow_vs2) < PRECISION_DTDV_DEF) {
      tst2 = 0;
    } else {
      tst2 = sqrt(
        (a1 - a2) * (a1 * (-2 * a2 * pe + 2 * a2 * ps + pow_ve2) - a2 * pow_vs2));  
    }
      
    // case 1
    const Scalar t1_1 = (-(a1 * vs) + a2 * vs + tst1) / (a1 * (a1 - a2));
    const Scalar t2_1 = -((-(a1 * ve) + a2 * ve + tst2) / ((a1 - a2) * a2));

    // case 2
    const Scalar t1_2 = -((a1 * vs - a2 * vs + tst1) / (a1 * (a1 - a2)));
    const Scalar t2_2 = (a1 * ve - a2 * ve + tst2) / ((a1 - a2) * a2);

    if (std::isfinite(t1_1) and std::isfinite(t2_1) and vs + a1 * t1_1 >= v_min and vs + a1 * t1_1 <= v_max and
        t1_1 > -PRECISION_PMM_VALUES and t2_1 > -PRECISION_PMM_VALUES and
        fabs(desired_time - (t1_1 + t2_1)) < fabs(desired_time - (t1 + t2 + t3))) {
      // clip time values
      t1 = std::max(t1_1, 0.0);
      t2 = 0.0;
      t3 = std::max(t2_1, 0.0);

      sol_idx_ = 0;

      if (calc_gradient) {
        const Scalar d_t_dvs = (a1 * vs - a2 * vs - tst1) / (a1 * tst1);
        const Scalar d_t_dve = (-(a1 * ve) + a2 * ve + tst2) / (a2 * tst1);

        dt_dvs = d_t_dvs;
        dt_dve = d_t_dve;
      }
      used_acc1 = a1;
      used_acc2 = a2;
    }
    else if (std::isfinite(t1_2) and std::isfinite(t2_2) and vs + a1 * t1_2 >= v_min and vs + a1 * t1_2 <= v_max and
        t1_2 > -PRECISION_PMM_VALUES and t2_2 > -PRECISION_PMM_VALUES and
        fabs(desired_time - (t1_2 + t2_2)) < fabs(desired_time - (t1 + t2 + t3))) {
      t1 = std::max(t1_2, 0.0);
      t2 = 0.0;
      t3 = std::max(t2_2, 0.0);

      sol_idx_ = 1;

      if (calc_gradient) {
        // dt/da == gradient we can use to optimize time
        const Scalar d_t_da1 =
          (a1 * (-2 * a2 * pe + 2 * a2 * ps + pow_ve2 + pow_vs2) +
          2 * vs * (-(a2 * vs) + tst1)) /
          (2 * pow_a1_2 * tst1);
        const Scalar d_t_da2 = (2 * a1 * a2 * (pe - ps) - 2 * a1 * pow_ve2 +
                                a2 * (pow_ve2 + pow_vs2) - 2 * ve * tst2) /
                              (2 * pow_a2_2 * tst1);


        const Scalar d_t_dvs = (-(a1 * vs) + a2 * vs - tst1) / (a1 * tst1);
        const Scalar d_t_dve = (a1 * ve - a2 * ve + tst2) / (a2 * tst1);

        dt_dvs = d_t_dvs;
        dt_dve = d_t_dve;
      }
      used_acc1 = a1;
      used_acc2 = a2;
    }
  }

  if (t1 != MAX_SCALAR and t2 != MAX_SCALAR and t3 != MAX_SCALAR and t1 >= 0 and t2 >= 0 and t3 >= 0) {
    exists_ = true;
    t_ = Vector<3>(t1, t2, t3);
    p_(1) = ps + t1 * vs + 0.5 * used_acc1 * t1 * t1;
    v_(1) = vs + used_acc1 * t1;
    p_(2) = p_(1) + t2 * v_(1);
    a_(0) = used_acc1;
    a_(1) = used_acc2;
    dt_dvs_ = dt_dvs;
    dt_dve_ = dt_dve;

    // if (fabs(v_(1)) > v_max+PRECISION_PMM_VALUES) {
    //   std::cout << "DANGER " << v_(1) << " " << sol_idx_ << std::endl;
    // } 

    if (check_result) {
      const Scalar ve_tst = v_(1) + a_(1) * t3;
      const Scalar pe_tst = p_(2) + t3 * v_(1) + 0.5 * a_(1) * t3 * t3;
      if (fabs(ve_tst - ve) > 1e-3 ||
          fabs(pe_tst - pe) > 1e-3) {
        std::cout << "wrong ve or pe oddi two acc" << std::endl;

        std::cout << "start ps " << ps << " vs " << vs << std::endl;
        
        std::cout << "ve_tst " << ve_tst << " ve " << ve << " diff " << fabs(ve_tst - ve) << std::endl;
        std::cout << "pe_tst " << pe_tst << " pe " << pe << " diff " << fabs(pe_tst - pe) << std::endl;

        std::cout << "t1 " << t1 << " t2 " << t2 << " t3 " << t3 << " a1 " << used_acc1
                  << " a2 " << used_acc2 << std::endl;
        std::cout << "sol_idx " << sol_idx_ << std::endl;
        exists_ = false;
      }
    }

  } else {
    exists_ = false;
    // std::cout << "failed " << ps << " to " << pe << " with " << a1_in << " " << a2_in << " and |v|  < " << v_max <<  std::endl;
  }
}

Scalar PMMTrajectory::computeSyncTrajectory(const Scalar total_time, const Scalar a_max1, const Scalar a_max2, const Scalar max_v) {

  v_max = max_v;
  v_min = -max_v;
  
  if ((!exists_) or (this->time() < 0.0)){
    std::cout << "ERROR: Invalid trajectory to be rescaled, exists: "<< exists_ << ", time: " << time() << ", axis: " << i_ << std::endl;
    exists_ = false;
    return -1;
  }

  if (total_time - this->time() < -PRECISION_PMM_VALUES) {
    std::cout << "ERROR: Desired time: " << total_time << ", is smaller than current time: " << this->time() << std::endl;
    exists_ = false;
    return -1;
  }

  if (fabs(p_(0)-p_(3)) < PRECISION_PMM_VALUES and fabs(v_(0)-v_(2)) < PRECISION_PMM_VALUES and fabs(v_(0)) < PRECISION_PMM_VALUES){
    exists_ = true;
    trajectory_type_ = ZERO_TRAJ;
    t_ = Vector<3>(total_time, 0, 0);
    v_(1) = v_(0);
    p_(1) = p_(2) = p_(0);
    a_(0) = a_(1) = 0.;
    return 0;
  }

  // given trajectory values
  const Scalar &ps = p_(0);
  const Scalar &pe = p_(3);
  const Scalar &vs = v_(0);
  const Scalar &ve = v_(2);
  Scalar used_a1 = a_max1;
  Scalar used_a2 = a_max2;

  std::vector<Vector<2>> test_acc_vec = {Vector<2>(a_max1, a_max2),
                                         Vector<2>(a_max2, a_max1)};

  for (Vector<2> a : test_acc_vec) {
    Scalar a1 = a(0);
    Scalar a2 = a(1);

    // recompute current trajectory to the input time using bang-singular-bang approach
    Scalar t1_3;
    Scalar t2_3;
    Scalar t3_3;
    Scalar t1_4;
    Scalar t2_4;
    Scalar t3_4;
    Scalar sc_3;
    Scalar sc_4;

    // BANG-SINGULAR-BANG

      t1_3 = - (2 * a2 * (pe - ps - total_time * v_max) * (v_max - vs)) / (-a1 * pow(v_max - ve, 2) + a2 * pow(v_max - vs, 2));
      t2_3 = (-a1 * (v_max - ve) * (-2 * pe + 2* ps + total_time * (v_max + ve)) + a2 * (v_max - vs) * (-2 * pe + 2 * ps + total_time * (v_max + vs))) / (a1 * pow(v_max - ve, 2) - a2 * pow(v_max - vs, 2));
      t3_3 = - (2 * a1 * (pe - ps - total_time * v_max) * (v_max - ve)) / (a1 * pow(v_max - ve, 2) - a2 * pow(v_max - vs, 2));   

      t1_4 = - (2 * a2 * (pe - ps - total_time * v_min) * (v_min - vs)) / (-a1 * pow(v_min - ve, 2) + a2 * pow(v_min - vs, 2));
      t2_4 = (-a1 * (v_min - ve) * (-2 * pe + 2* ps + total_time * (v_min + ve)) + a2 * (v_min - vs) * (-2 * pe + 2 * ps + total_time * (v_min + vs))) / (a1 * pow(v_min - ve, 2) - a2 * pow(v_min - vs, 2));
      t3_4 = - (2 * a1 * (pe - ps - total_time * v_min) * (v_min - ve)) / (a1 * pow(v_min - ve, 2) - a2 * pow(v_min - vs, 2));

      sc_3 = (a1 * pow(v_max - ve, 2) - a2 * pow(v_max - vs, 2)) / (2 * a1 * a2 * (pe - ps - total_time * v_max));
      
      sc_4 = (a1 * pow(v_min - ve, 2) - a2 * pow(v_min - vs, 2)) / (2 * a1 * a2 * (pe - ps - total_time * v_min));
 
    
  
    // select feasible solution
    if (std::isfinite(t1_3) and std::isfinite(t2_3) and std::isfinite(t3_3) and std::isfinite(sc_3) and sc_3 > 0 and
            sc_3 <= (1+PRECISION_PMM_VALUES) and t1_3 > -PRECISION_PMM_VALUES and t2_3 > -PRECISION_PMM_VALUES and t3_3 > -PRECISION_PMM_VALUES) {
      // clip time values
      Scalar t1 = std::max(0.0, t1_3);
      Scalar t2 = std::max(0.0, t2_3);
      Scalar t3 = std::max(0.0, t3_3);
      Scalar sc = sc_3;
      
      // calculate remaining trajectory values
      used_a1 = sc*a1;
      used_a2 = sc*a2;
      exists_ = true;
      trajectory_type_ = SYNC_TRAJ;
      t_ = Vector<3>(t1, t2, t3);
      p_(1) = ps + t1 * vs + 0.5 * used_a1 * t1 * t1;
      v_(1) = vs + used_a1 * t1;
      p_(2) = p_(1) + t2 * v_(1);
      a_(0) = used_a1;
      a_(1) = used_a2;
      sol_idx_ = 2;

      return 0;
    }

    if (std::isfinite(t1_4) and std::isfinite(t2_4) and std::isfinite(t3_4) and std::isfinite(sc_4) and sc_4 > 0 and
            sc_4 <= (1+PRECISION_PMM_VALUES) and t1_4 > -PRECISION_PMM_VALUES and t2_4 > -PRECISION_PMM_VALUES and t3_4 > -PRECISION_PMM_VALUES) {
      // clip time values
      Scalar t1 = std::max(0.0, t1_4);
      Scalar t2 = std::max(0.0, t2_4);
      Scalar t3 = std::max(0.0, t3_4);
      Scalar sc = sc_4;

      // calculate remaining trajectory values
      used_a1 = sc*a1;
      used_a2 = sc*a2;
      exists_ = true;
      trajectory_type_ = SYNC_TRAJ;
      t_ = Vector<3>(t1, t2, t3);
      p_(1) = ps + t1 * vs + 0.5 * used_a1 * t1 * t1;
      v_(1) = vs + used_a1 * t1;
      p_(2) = p_(1) + t2 * v_(1);
      a_(0) = used_a1;
      a_(1) = used_a2;
      sol_idx_ = 3;

      return 0;
    }

    Scalar t1_1, t1_2, t2_1, t2_2, sc_1, sc_2;

    // BANG-BANG
    if (fabs(vs-ve) < PRECISION_PMM_VALUES){
      // simplified solution preventing division by zero
      t1_1 = -(a2*total_time)/(a1 - a2);
      t1_2 = t1_1;
      t2_1 = (a1*total_time)/(a1 - a2);
      t2_2 = t2_1;
      sc_1 = ((2*a1 - 2*a2)*(ps - pe + total_time*ve))/(a1*a2*pow(total_time,2));
      sc_2 = sc_1;

    } else {
      const Scalar sigma_1 = a1*a2*sqrt(((a1 - a2)*(a1*pow(pe,2) - a2*pow(pe,2) + a1*pow(ps,2) - a2*pow(ps,2) + a1*pow(total_time,2)*pow(ve,2) - a2*pow(total_time,2)*pow(vs,2) - 2*a1*pe*ps + 2*a2*pe*ps - 2*a1*pe*total_time*ve + 2*a1*ps*total_time*ve + 2*a2*pe*total_time*vs - 2*a2*ps*total_time*vs))/(pow(a1,2)*pow(a2,2)));

      t1_1 = -(a1*pe - a2*pe - a1*ps + a2*ps + sigma_1 - a1*total_time*ve + a2*total_time*ve)/((a1 - a2)*(ve - vs));
      t1_2 = (a2*pe - a1*pe + a1*ps - a2*ps + sigma_1 + a1*total_time*ve - a2*total_time*ve)/((a1 - a2)*(ve - vs));

      t2_1 = (a1*pe - a2*pe - a1*ps + a2*ps + sigma_1 - a1*total_time*vs + a2*total_time*vs)/((a1 - a2)*(ve - vs));
      t2_2 = -(a2*pe - a1*pe + a1*ps - a2*ps + sigma_1 + a1*total_time*vs - a2*total_time*vs)/((a1 - a2)*(ve - vs));

      sc_1 = (a2*pe - a1*pe + a1*ps - a2*ps + sigma_1 + a1*total_time*ve - a2*total_time*vs)/(a1*a2*pow(total_time,2));
      sc_2 = -(a1*pe - a2*pe - a1*ps + a2*ps + sigma_1 - a1*total_time*ve + a2*total_time*vs)/(a1*a2*pow(total_time,2));
    }

    // std::cout << "sync " << t1_1 << " " << t2_1 << " " << sc_1 << std::endl;
    // std::cout << "sync " << t1_2 << " " << t2_2 << " " << sc_2 << std::endl;

    if (std::isfinite(t1_1) and std::isfinite(t2_1) and std::isfinite(sc_1) and sc_1 > 0 and
            /*sc_1 <= (1+PRECISION_PMM_VALUES) and*/ t1_1 > -PRECISION_PMM_VALUES and t2_1 > -PRECISION_PMM_VALUES and
            vs + sc_1*a1 * t1_1 <= v_max and vs + sc_1*a1 * t1_1 >= v_min) {
      // clip time values
      Scalar t1 = std::max(0.0, t1_1);
      Scalar t2 = std::max(0.0, t2_1);
      Scalar sc = sc_1;

      // calculate remaining trajectory values
      used_a1 = sc*a1;
      used_a2 = sc*a2;
      exists_ = true;
      trajectory_type_ = SYNC_TRAJ;
      t_ = Vector<3>(t1, 0, t2);
      v_(1) = vs + used_a1 * t1;
      p_(1) = ps + t1 * vs + 0.5 * used_a1 * t1 * t1;
      p_(2) = p_(1);
      a_(0) = used_a1;
      a_(1) = used_a2;
      sol_idx_ = 0;

      return 0;
    }

    if (std::isfinite(t1_2) and std::isfinite(t2_2) and std::isfinite(sc_2) and sc_2 > 0 and
            /*sc_2 <= (1+PRECISION_PMM_VALUES) and*/ t1_2 > -PRECISION_PMM_VALUES and t2_2 > -PRECISION_PMM_VALUES and 
            vs + sc_2*a1 * t1_2 <= v_max and vs + sc_2*a1 * t1_2 >= v_min) {
      // clip time values
      Scalar t1 = std::max(0.0, t1_2);
      Scalar t2 = std::max(0.0, t2_2);
      Scalar sc = sc_2;

      // calculate remaining trajectory values
      used_a1 = sc*a1;
      used_a2 = sc*a2;
      exists_ = true;
      trajectory_type_ = SYNC_TRAJ;
      t_ = Vector<3>(t1, 0, t2);
      v_(1) = vs + used_a1 * t1;
      p_(1) = ps + t1 * vs + 0.5 * used_a1 * t1 * t1;
      p_(2) = p_(1);
      a_(0) = used_a1;
      a_(1) = used_a2;
      sol_idx_ = 1;
      return 0;
    }

  }
  

  // compute the smallest possible time the axis is able to sync
  for (Vector<2> a : test_acc_vec) {
    Scalar a11 = a(1);
    Scalar a12 = a(0);

    Scalar t1_3 = (v_max - vs) / a11;
    Scalar t2_3 = (a11 * (2 * a12 * (pe - ps) + pow(v_max, 2) - pow(ve, 2)) + a12 * (pow(vs, 2) - pow(v_max, 2))) / 2 / a11 / a12 / v_max;
    Scalar t3_3 = (ve - v_max) / a12;

    Scalar t1_4 = (v_min - vs) / a11;
    Scalar t2_4 = (a11 * (2 * a12 * (pe - ps) + pow(v_min, 2) - pow(ve, 2)) + a12 * (pow(vs, 2) - pow(v_min, 2))) / 2 / a11 / a12 / v_min;
    Scalar t3_4 = (ve - v_min) / a12;

    Scalar feasible_sync_time = -1;

    if (std::isfinite(t1_3) and std::isfinite(t2_3) and std::isfinite(t3_3) and
    t1_3 > -PRECISION_PMM_VALUES and t2_3 > -PRECISION_PMM_VALUES and t3_3 > -PRECISION_PMM_VALUES and (t1_3 + t2_3 + t3_3) >= total_time){
        // clip time values
        Scalar t1 = std::max(0.0, t1_3);
        Scalar t2 = std::max(0.0, t2_3);
        Scalar t3 = std::max(0.0, t3_3);

        // calculate remaining trajectory values
        exists_ = true;
        trajectory_type_ = MIN_TRAJ;
        t_ = Vector<3>(t1, t2, t3);
        p_(1) = ps + t1 * vs + 0.5 * a11 * t1 * t1;
        v_(1) = vs + a11 * t1;
        p_(2) = p_(1) + t2 * v_(1);
        a_(0) = a11;
        a_(1) = a12;
        sol_idx_ = 2;
        // return new slowest traj time
        feasible_sync_time = t1 + t2 + t3;
        return feasible_sync_time;
        
    } else if (std::isfinite(t1_4) and std::isfinite(t2_4) and std::isfinite(t3_4) and
    t1_4 > -PRECISION_PMM_VALUES and t2_4 > -PRECISION_PMM_VALUES and t3_4 > -PRECISION_PMM_VALUES and (t1_4 + t2_4 + t3_4) >= total_time){        
        // clip time values
        Scalar t1 = std::max(0.0, t1_4);
        Scalar t2 = std::max(0.0, t2_4);
        Scalar t3 = std::max(0.0, t3_4);

        // calculate remaining trajectory values
        exists_ = true;
        trajectory_type_ = MIN_TRAJ;
        t_ = Vector<3>(t1, t2, t3);
        p_(1) = ps + t1 * vs + 0.5 * a11 * t1 * t1;
        v_(1) = vs + a11 * t1;
        p_(2) = p_(1) + t2 * v_(1);
        a_(0) = a11;
        a_(1) = a12;
        sol_idx_ = 2;
        // return new slowest traj time
        feasible_sync_time = t1 + t2 + t3;
        return feasible_sync_time;      

    }

    Scalar sigma_1;
    if (-(a12*pow(vs,2) - a11*pow(ve,2) - 2*a11*a12*ps + 2*a11*a12*pe)/(a11 - a12) < -PRECISION_PMM_VALUES) {
      continue;
    } else if (fabs(-(a12*pow(vs,2) - a11*pow(ve,2) - 2*a11*a12*ps + 2*a11*a12*pe)/(a11 - a12)) < PRECISION_PMM_VALUES) {
      sigma_1 = 0;
    } else {
      sigma_1 = sqrt(-(a12*pow(vs,2) - a11*pow(ve,2) - 2*a11*a12*ps + 2*a11*a12*pe)/(a11 - a12));
    }

    Scalar t1_1 = -(vs + sigma_1)/a11;
    Scalar t2_1 = (ve + sigma_1)/a12;
    Scalar t1_2 = - (vs - sigma_1)/a11;
    Scalar t2_2 = (ve - sigma_1)/a12;

    // std::cout << t1_1 << " " << t2_1 << std::endl;
    // std::cout << t1_2 << " " << t2_2 << std::endl;

    if (std::isfinite(t1_1) and std::isfinite(t2_1) and vs + a11 * t1_1 <= v_max and vs + a11 * t1_1 >= v_min and
    t1_1 > -PRECISION_PMM_VALUES and t2_1 > -PRECISION_PMM_VALUES and (t1_1 + t2_1) >= total_time) {
        // clip time values
        Scalar t1 = std::max(0.0, t1_1);
        Scalar t2 = std::max(0.0, t2_1);

        // calculate remaining trajectory values
        exists_ = true;
        trajectory_type_ = MIN_TRAJ;
        t_ = Vector<3>(t1, 0, t2);
        v_(1) = vs + a11 * t1;
        p_(1) = ps + t1 * vs + 0.5 * a11 * t1 * t1;
        p_(2) = p_(1);
        a_(0) = a11;
        a_(1) = a12;
        sol_idx_ = 0;
        // return new slowest traj time
        feasible_sync_time = t1 + t2;
        return feasible_sync_time;
    } else if (std::isfinite(t1_2) and std::isfinite(t2_2) and vs + a11 * t1_2 <= v_max and vs + a11 * t1_2 >= v_min and
    t1_2 > -PRECISION_PMM_VALUES and t2_2 > -PRECISION_PMM_VALUES and (t1_2 + t2_2) >= total_time){
        // clip time values
        Scalar t1 = std::max(0.0, t1_2);
        Scalar t2 = std::max(0.0, t2_2);

        // calculate remaining trajectory values
        exists_ = true;
        trajectory_type_ = MIN_TRAJ;
        t_ = Vector<3>(t1, 0, t2);
        v_(1) = vs + a11 * t1;
        p_(1) = ps + t1 * vs + 0.5 * a11 * t1 * t1;
        p_(2) = p_(1);
        a_(0) = a11;
        a_(1) = a12;
        sol_idx_ = 1;
        // return new slowest traj time
        feasible_sync_time = t1 + t2;
        return feasible_sync_time;
    }
  }
  return -1;
}

Vector<3> PMMTrajectory::state_in_time(const Scalar time_in_tr) const {
  Scalar pos, vel, acc;

  if (time_in_tr < t_(0)) {  // const acc part with a1
    pos = p_(0) + v_(0) * time_in_tr + 0.5 * a_(0) * time_in_tr * time_in_tr;
    vel = v_(0) + a_(0) * time_in_tr;
    acc = a_(0);
  } else if (time_in_tr < t_(0) + t_(1)) {  // const vel part with a=0
    const Scalar time_part = (time_in_tr - t_(0));
    pos = p_(1) + v_(1) * time_part;
    vel = v_(1);
    acc = 0.0;
  } else if (time_in_tr < time()) {  // const vel part with a2
    const Scalar time_part = (time_in_tr - t_(0) - t_(1));
    pos = p_(2) + v_(1) * time_part + 0.5 * a_(1) * time_part * time_part;
    vel = v_(1) + a_(1) * time_part;
    acc = a_(1);
  } else {  // return the last state
    pos = p_(3);
    vel = v_(2);
    acc = a_(1);
  }

  if (t_(0) == 0 && t_(1) == 0 && t_(2) > 0) {
    acc = a_(1);
  } else if (t_(0) > 0 and t_(1) == 0 and t_(2) == 0) {
    acc = a_(0);
  }

  return Vector<3>(pos, vel, acc);
}

std::tuple<PMMTrajectory, PMMTrajectory> PMMTrajectory::split_in_time(
  const Scalar time_in_tr) {
  Vector<3> state = state_in_time(time_in_tr);

  PMMTrajectory bef;
  bef.exists_ = true;
  bef.i_ = i_;
  bef.a_ = a_;
  PMMTrajectory aft;
  aft.exists_ = true;
  aft.i_ = i_;
  aft.a_ = a_;

  if (time_in_tr <= t_(0)) {
    bef.t_(0) = time_in_tr;
    bef.t_(1) = bef.t_(2) = .0;
    bef.t_(2) = .0;
    bef.p_(0) = p_(0);
    bef.p_(1) = bef.p_(2) = bef.p_(3) = state(0);
    bef.v_(0) = v_(0);
    bef.v_(1) = bef.v_(2) = state(1);

    aft.t_(0) = t_(0) - time_in_tr;
    aft.t_(1) = t_(1);
    aft.t_(2) = t_(2);
    aft.p_(0) = state(0);
    aft.p_(1) = p_(1);
    aft.p_(2) = p_(2);
    aft.p_(3) = p_(3);
    aft.v_(0) = state(1);
    aft.v_(1) = v_(1);
    aft.v_(2) = v_(2);
  } else if (time_in_tr <= t_(0) + t_(1)) {
    bef.t_(0) = t_(0);
    bef.t_(1) = time_in_tr - t_(0);
    bef.t_(2) = .0;
    bef.p_(0) = p_(0);
    bef.p_(1) = p_(1);
    bef.p_(2) = bef.p_(3) = state(0);
    bef.v_(0) = v_(0);
    bef.v_(1) = v_(1);
    bef.v_(2) = state(1);

    aft.t_(0) = .0;
    aft.t_(1) = (t_(0) + t_(1)) - time_in_tr;
    aft.t_(2) = t_(2);
    aft.p_(0) = aft.p_(1) = state(0);
    aft.p_(2) = p_(2);
    aft.p_(3) = p_(3);
    aft.v_(0) = aft.v_(1) = state(1);
    aft.v_(2) = v_(2);
  } else if (time_in_tr <= t_(0) + t_(1) + t_(2)) {
    bef.t_ = t_;
    bef.t_(2) = time_in_tr - (t_(0) + t_(1));
    bef.p_ = p_;
    bef.p_(3) = state(0);
    bef.v_ = v_;
    bef.v_(2) = state(1);

    aft.t_(0) = aft.t_(1) = .0;
    aft.t_(2) = (t_(0) + t_(1) + t_(2)) - time_in_tr;
    aft.p_(0) = aft.p_(1) = aft.p_(2) = state(0);
    aft.p_(3) = p_(3);
    aft.v_(0) = aft.v_(1) = state(1);
    aft.v_(2) = v_(2);
  } else {
    bef.t_ = t_;
    bef.p_ = p_;
    bef.v_ = v_;

    aft.t_(0) = aft.t_(1) = aft.t_(2) = .0;
    aft.p_(0) = aft.p_(1) = aft.p_(3) = state(0);
    aft.v_(0) = aft.v_(1) = aft.v_(2) = state(1);
  }
  return {bef, aft};
}

Scalar PMMTrajectory::acc_in_time(const Scalar time_in_tr) const {
  Scalar acc;

  if (time_in_tr < t_(0)) {  // const acc part with a1
    acc = a_(0);
  } else if (time_in_tr < t_(0) + t_(1)) {  // const vel part with a=0
    // const Scalar time_part = (time_in_tr - t_(0));
    acc = 0.0;
  } else if (time_in_tr < time()) {  // const vel part with a2
    // const Scalar time_part = (time_in_tr - t_(0) - t_(1));
    acc = a_(1);
  } else {  // return the last state
    acc = a_(1);
  }

  // if (t_(0) == 0 && t_(1) == 0 && t_(2) > 0) {
  //   acc = a_(1);
  // } else if (t_(0) > 0 and t_(1) == 0 and t_(2) == 0) {
  //   acc = a_(0);
  // }

  return acc;
}

Scalar PMMTrajectory::vel_in_time(const Scalar time_in_tr) const {
  Scalar vel;

  if (time_in_tr < t_(0)) {  // a1
    vel = v_(0) + a_(0)*time_in_tr;
  } else if (time_in_tr < t_(0) + t_(1)) {  // const vel part with a=0
    // const Scalar time_part = (time_in_tr - t_(0));
    vel = v_(1);
  } else if (time_in_tr < time()) {  // a2
    const Scalar time_part = (time_in_tr - t_(0) - t_(1));
    vel = v_(1) + a_(1)*time_part;
  } else {  // return the last state
    vel = v_(2);
  }

  return vel;
}

void PMMTrajectory::copy_trajectory(const PMMTrajectory &in) {
  exists_ = in.exists_;
  t_ = in.t_;
  p_ = in.p_;
  v_ = in.v_;
  a_ = in.a_;
  i_ = in.i_;
  dt_dvs_ = in.dt_dvs_;
  dt_dve_ = in.dt_dve_;
  sol_idx_ = in.sol_idx_;
  trajectory_type_ = in.trajectory_type_;
}

void PMMTrajectory::save_to_file(std::string filename) {
  std::ofstream myfile;
  myfile.open(filename.c_str());
  if (myfile.is_open()) {
    myfile << "a1,a2,t1,t2,t3,p0,p1,p2,p3,v0,v1,"
              "v3,exists,axis"
           << std::endl;
    myfile << a_(0) << "," << a_(1) << "," << t_(0) << "," << t_(1) << ","
           << t_(2) << "," << p_(0) << "," << p_(1) << "," << p_(2) << ","
           << p_(3) << "," << v_(0) << "," << v_(1) << "," << v_(2) << ","
           << exists_ << "," << i_;
    myfile.close();
  }
}

Scalar PMMTrajectory::max_end_velocity_abs() const {
  /*what is the end velocity that is when only one acc part is used*/
  Scalar a = a_(0);
  if (v_(2) > v_(0)) {
    if (a < 0) a = a_(1);  // need positive acc
  } else {
    if (a > 0) a = a_(1);  // need negative acc
  }
  const Scalar vs_pow_2 = v_(0) * v_(0);
  Scalar max_v_a = 1.0 / (2.0 * a);
  Scalar max_v_c = p_(0) - p_(3) - vs_pow_2 / (2.0 * a);
  Scalar max_v_disc_sqrt = sqrt(-4 * max_v_a * max_v_c);
  Scalar max_v_abs = max_v_disc_sqrt / (2.0 * max_v_a);
  return max_v_abs;
}

Scalar PMMTrajectory::minRequiredAcc(const Scalar ps, const Scalar vs,
                                     const Scalar pe, const Scalar ve) {
  // required acc to be able to fulfill the task
  // i.e. accelerate between required velocities and in the same time not
  // overshoot the position derived from
  // t_change_vel = (ve-vs)/amax
  // and pe-ps = vs*t_change_vel + 0.5*amax*t_change_vel**2
  // from that what is the amax....
  if (fabs(pe - ps) < PRECISION_TRANS3D) {
    return MIN_ACC_REQ;
  } else {
    const Scalar pow_ve2 = ve * ve;
    const Scalar pow_vs2 = vs * vs;
    if (fabs(pow_ve2 - pow_vs2) < PRECISION_TRANS3D) {
      // does not need any acc basically
      return std::copysign(INITIAL_MIN_ACC, pe - ps);
    } else {
      const Scalar a_min_req = (-pow_ve2 + pow_vs2) / (2 * (-pe + ps));
      return a_min_req;
    }
  }
}

std::ostream &operator<<(std::ostream &o, const PMMTrajectory &f) {
  o << "maxacc: "
    << "t_tot:" << (f.time()) << ";t1:" << f.t_(0) << ";t2:" << f.t_(1)
    << ";t3:" << f.t_(2) << ";exists:" << f.exists_ << ";a1:" << f.a_(0)
    << ";a2:" << f.a_(1) << ";i:" << f.i_ << "\n";
  o << " \t: p0:" << f.p_(0) << ";p1:" << f.p_(1) << ";p2:" << f.p_(2)
    << ";p3:" << f.p_(3) << "\n";
  o << " \t: v0:" << f.v_(0) << ";v1:" << f.v_(1) << ";v3:" << f.v_(2) << "\n";
  o << " \t: dt/dvs:" << f.dt_dvs_ << ";dt/dve:" << f.dt_dve_ << "\n";
  return o;
}
} // namespace pmm