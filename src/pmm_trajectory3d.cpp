/**
 * @file pmm_trajectory3d.cpp
 * @author Matej Novosad (novosma2@fel.cvut.cz), Robert Penicka (penicrob@fel.cvut.cz), Krystof Teissing (teisskry@gmail.com)
 * @version 0.1
 * @date 2024-09-20
 * 
 */

#include "pmm_trajectory3d.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>

namespace pmm {

PointMassTrajectory3D::PointMassTrajectory3D() {}


PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from, const QuadState &to,
                        const Scalar max_acc_norm, const Scalar max_vel_norm, const int max_iter,
                        const double precision_acc_limit, const bool debug) {
  // B = T + GVEC , |T|=max_acc_norm
  // initially B equal per axis with b_x=b_y=b_z -> |B-GVEC|^2 = |T|^2
  // -> 3*bs_x^2 + 2*g*a_x + g^2 - |T|^2 = 0 --> roots are the possible acc

  // const Scalar precision_acc_limit = 0.1;

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
    std::cout << "that is happening ever?" << std::endl;
  }

  Vector<3> max_per_axis_acc_vec = per_axis_acc_vec;
  Vector<3> min_per_axis_acc_vec = -per_axis_acc_vec + 2*GVEC;
  Vector<3> max_per_axis_vel_vec = Vector<3>::Constant(max_vel_norm / sqrt(3));

  Vector<3> endpoint_velocity = from.v.cwiseAbs().cwiseMax(to.v.cwiseAbs());
  // if (abs(endpoint_velocity.maxCoeff()) > max_vel_norm / sqrt(3)) {
  //   std::cout << "WARNING: start or end velocity in one of the axis is greater than maximum velocity / sqrt(3)!" << std::endl;
  // }
  Vector<3> max_velocity = endpoint_velocity.cwiseAbs().cwiseMax(max_per_axis_vel_vec.cwiseAbs());

  PointMassTrajectory3D pmm3d(from, to, max_per_axis_acc_vec, min_per_axis_acc_vec, max_velocity, true, true);

  max_velocity = pmm3d.max_velocity().normalized() * max_vel_norm;
  // max_velocity = endpoint_velocity.cwiseAbs().cwiseMax(max_velocity.cwiseAbs());
  
  Vector<3> start_acc = pmm3d.start_acc();
  Vector<3> end_acc = pmm3d.end_acc();
  Scalar start_thrust = (start_acc - GVEC).norm();
  Scalar end_thrust = (end_acc - GVEC).norm();

  Vector<3> biggest_acc = start_thrust > end_thrust ? start_acc : end_acc;
  Scalar largest_thrust = std::max(start_thrust, end_thrust);

  Scalar duration = pmm3d.time();

  if (largest_thrust > max_acc_norm) {
    duration = MAX_SCALAR;
    largest_thrust = 0;
  } else {
    x_ = pmm3d.x_;
    y_ = pmm3d.y_;
    z_ = pmm3d.z_;

    min_axis_trajectory_duration_ = pmm3d.min_axis_trajectory_duration_;
  }

  int iter = 0;
  while ((fabs(largest_thrust - max_acc_norm) > precision_acc_limit or
         largest_thrust >= max_acc_norm) and iter < max_iter) {
    iter++;
    // B = T + GVEC , |T|=max_acc_norm, |B-GVEC|^2 = |T|^2
    // scale the T parts by same factor k ->
    // k^2*b_x^2 + k^2*b_y^2 + (k*b_z + g)^2 - |T|^2 = 0 -> find k
    // (b_x^2 + b_y^2 + b_z^2) * k^2 + (2*g*b_z)*k


    const Scalar a_dis = biggest_acc.squaredNorm();
    const Scalar b_dis = 2 * biggest_acc(2) * G;
    const Scalar c_dis = G * G - (max_acc_norm * max_acc_norm);
    const Scalar k1 = (-b_dis + sqrt(b_dis * b_dis - 4 * a_dis * c_dis)) / (2 * a_dis);
    const Scalar k2 = (-b_dis - sqrt(b_dis * b_dis - 4 * a_dis * c_dis)) / (2 * a_dis);
    Vector<3> thrust_acc_new_k1 = k1 * biggest_acc - GVEC;
    Vector<3> thrust_acc_new_k2 = k2 * biggest_acc - GVEC;
    Vector<3> max_acc_new1 = thrust_acc_new_k1 + GVEC;
    Vector<3> max_acc_new2 = -thrust_acc_new_k1 + GVEC;

    // for (int a=0;a<3;a++) {
    //   if (fabs(max_acc_new1[a]) < MIN_ACC_REQ) {
    //     max_acc_new1[a] = std::copysign(MIN_ACC_REQ, max_acc_new1[a]);
    //   }
    //   if (fabs(max_acc_new2[a]) < MIN_ACC_REQ) {
    //     max_acc_new2[a] = std::copysign(MIN_ACC_REQ, max_acc_new2[a]);
    //   }
    // }

    if (debug) {
      std::cout << std::endl << "Iter: " << iter << std::endl;
      std::cout << "Thrust: " << largest_thrust << std::endl;
      std::cout << "max_acc_new1: " << max_acc_new1[0] << "," << max_acc_new1[1] << ","<< max_acc_new1[2] << std::endl;
      std::cout << "max_acc_new2: " << max_acc_new2[0] << "," << max_acc_new2[1] << ","<< max_acc_new2[2] << std::endl;
    }

    pmm3d = PointMassTrajectory3D(from, to, max_acc_new1, max_acc_new2, max_velocity, true, true);

    if (!pmm3d.exists()) {
      max_velocity = endpoint_velocity.cwiseAbs().cwiseMax(max_velocity.cwiseAbs());
      // redistribute(biggest_acc);
      continue;
    }

    max_velocity = pmm3d.max_velocity().normalized() * max_vel_norm;
    max_velocity = endpoint_velocity.cwiseAbs().cwiseMax(max_velocity.cwiseAbs());

    if (debug){
      std::cout << "T new: " << pmm3d.time() << std::endl;
    }

    if (debug) {
      Vector<3> start_body_acc = (pmm3d.start_acc() - GVEC);
      Vector<3> end_body_acc = (pmm3d.end_acc() - GVEC);
      std::cout << std::endl << "Start thrust: " << (pmm3d.start_acc() - GVEC).norm() << std::endl;
      std::cout <<  "Start body acc: " << start_body_acc[0] << "," << start_body_acc[1] << ","<< start_body_acc[2] << std::endl;
      std::cout << "End thrust: " << (pmm3d.end_acc() - GVEC).norm() << std::endl;
      std::cout <<  "End body acc: " << end_body_acc[0] << "," << end_body_acc[1] << ","<< end_body_acc[2] << std::endl;
    }

    start_acc = pmm3d.start_acc();
    end_acc = pmm3d.end_acc();
    start_thrust = (start_acc - GVEC).norm();
    end_thrust = (end_acc - GVEC).norm();

    if (start_thrust > end_thrust) {
      biggest_acc = start_acc;
      largest_thrust = start_thrust;
    } else {
      biggest_acc = end_acc;
      largest_thrust = end_thrust;
    }

    // if (pmm3d.x_.trajectory_type_ == MIN_TRAJ){
    //   biggest_acc[0] += std::copysign(max_acc_norm*0.01, biggest_acc[0]);
    // }
    // else if (pmm3d.y_.trajectory_type_ == MIN_TRAJ){
    //   biggest_acc[1] += std::copysign(max_acc_norm*0.01, biggest_acc[1]);
    // }
    // else if (pmm3d.z_.trajectory_type_ == MIN_TRAJ){
    //   biggest_acc[2] += std::copysign(max_acc_norm*0.01, biggest_acc[2]);
    // }

    if (pmm3d.time() < duration && largest_thrust <= max_acc_norm && pmm3d.exists()) {
      x_ = pmm3d.x_;
      y_ = pmm3d.y_;
      z_ = pmm3d.z_;

      min_axis_trajectory_duration_ = pmm3d.min_axis_trajectory_duration_;

      duration = pmm3d.time();
    }
  }

  if (debug){
    std::cout << "Iter: " << iter << std::endl;
    std::cout << "Thrust: " << largest_thrust << std::endl;
  }

  // x_ = pmm3d.x_;
  // y_ = pmm3d.y_;
  // z_ = pmm3d.z_;

  // min_axis_trajectory_duration_ = pmm3d.min_axis_trajectory_duration_;
}

/*
Used by PMM-MG-Trajectory
version that converges to thrust limit by iterative increasing scaled (to match
time) to the acc norm WITH DRAG
this version scales time by default
*/
PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from, const QuadState &to,
                        const Scalar max_acc_norm, const Scalar max_vel_norm, const bool drag, const int max_iter,
                        const double precision_acc_limit, const bool debug) {
  // B = T + GVEC , |T|=max_acc_norm
  // initially B equal per axis with b_x=b_y=b_z -> |B-GVEC|^2 = |T|^2
  // -> 3*bs_x^2 + 2*g*a_x + g^2 - |T|^2 = 0 --> roots are the possible acc

  if (!drag) {
    PointMassTrajectory3D pmm3d(from, to, max_acc_norm, max_vel_norm, max_iter, precision_acc_limit, debug);
    x_ = pmm3d.x_;
    y_ = pmm3d.y_;
    z_ = pmm3d.z_;

    min_axis_trajectory_duration_ = pmm3d.min_axis_trajectory_duration_;
    return;
  }

  // const Scalar precision_acc_limit = 0.1;

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
    std::cout << "that is happening ever? " << std::endl;
  }

  Vector<3> max_per_axis_acc_vec = per_axis_acc_vec;
  Vector<3> min_per_axis_acc_vec = -per_axis_acc_vec + 2*GVEC;
  Vector<3> max_per_axis_vel_vec = Vector<3>::Constant(max_vel_norm / sqrt(3));

  Vector<3> endpoint_velocity = from.v.cwiseAbs().cwiseMax(to.v.cwiseAbs());
  Vector<3> max_velocity = endpoint_velocity.cwiseAbs().cwiseMax(max_per_axis_vel_vec.cwiseAbs());

  PointMassTrajectory3D pmm3d(from, to, max_per_axis_acc_vec, min_per_axis_acc_vec, max_velocity, true, true);

  max_velocity = pmm3d.max_velocity().normalized() * max_vel_norm;

  Vector<3> switch_times = pmm3d.switch_times();

  Vector<3> acc1 = pmm3d.acc_in_time(switch_times(0));
  Vector<3> acc2 = pmm3d.acc_in_time(switch_times(1));
  Vector<3> acc3 = pmm3d.acc_in_time(switch_times(2));
  Vector<3> acc4 = pmm3d.end_acc();

  Vector<3> drag1 = pmm3d.drag(acc1 - GVEC, switch_times(0));
  Vector<3> drag2 = pmm3d.drag(acc2 - GVEC, switch_times(1));
  Vector<3> drag3 = pmm3d.drag(acc3 - GVEC, switch_times(2));
  Vector<3> drag4 = pmm3d.drag(acc4 - GVEC, pmm3d.time());

  Scalar thrust1 = (acc1 - drag1 - GVEC).norm();
  Scalar thrust2 = (acc2 - drag2 - GVEC).norm();
  Scalar thrust3 = (acc3 - drag3 - GVEC).norm();
  Scalar thrust4 = (acc4 - drag4 - GVEC).norm();

  Scalar largest_thrust = std::max(std::max(thrust1, thrust2), std::max(thrust3, thrust4));

  Scalar duration = pmm3d.time();

  if (largest_thrust > max_acc_norm) {
    duration = MAX_SCALAR;
    largest_thrust = 0;
  } else {
    x_ = pmm3d.x_;
    y_ = pmm3d.y_;
    z_ = pmm3d.z_;

    min_axis_trajectory_duration_ = pmm3d.min_axis_trajectory_duration_;
  }
  
  int iter = 0;
  while (duration >= pmm3d.time() and (fabs(largest_thrust - max_acc_norm) > precision_acc_limit or largest_thrust > max_acc_norm) and
         (iter < max_iter)) {
    iter++;

    // B = T + GVEC , |T|=max_acc_norm, |B-GVEC|^2 = |T|^2
    // scale the T parts by same factor k ->
    // k^2*b_x^2 + k^2*b_y^2 + (k*b_z + g)^2 - |T|^2 = 0 -> find k
    // (b_x^2 + b_y^2 + b_z^2) * k^2 + (2*g*b_z)*k

    // (k*b_x + d_x)^2 + (k*b_y + d_y)^2 + (k*b_z + d_z + g)^2 - |T|^2 = 0 -> find k

    // no drag
    // const Scalar a_dis = biggest_acc.squaredNorm();
    // const Scalar b_dis = 2 * biggest_acc(2) * G;
    // const Scalar c_dis = G * G - (max_acc_norm * max_acc_norm);

    Vector<3> max_acc_new1 = {2*max_acc_norm, 2*max_acc_norm, 2*max_acc_norm};
    Vector<3> max_acc_new2 = -max_acc_new1;

    const Scalar a_dis1 = acc1.squaredNorm();
    const Scalar b_dis1 = - 2 * acc1(0) * drag1(0) - 2 * acc1(1) * drag1(1) - 2 * acc1(2) * drag1(2) + 2 * acc1(2) * G;
    const Scalar c_dis1 = drag1.squaredNorm() - 2 * drag1(2) * G + G * G - (max_acc_norm * max_acc_norm);
    const Scalar k1 = (-b_dis1 + sqrt(b_dis1 * b_dis1 - 4 * a_dis1 * c_dis1)) / (2 * a_dis1);
    Vector<3> thrust_acc_new_k1 = k1 * acc1 - drag1 - GVEC;
    Vector<3> max_acc_new11 = thrust_acc_new_k1 + drag1 + GVEC;

    for (int a=0;a<3;a++) {
      if (max_acc_new11(a) > 0 && max_acc_new11(a) < max_acc_new1(a)) {
        max_acc_new1(a) = max_acc_new11(a);
      } else if  (max_acc_new11(a) < 0 && max_acc_new11(a) > max_acc_new2(a)) {
        max_acc_new2(a) = max_acc_new11(a);
      }
    }

    const Scalar a_dis2 = acc2.squaredNorm();
    const Scalar b_dis2 = - 2 * acc2(0) * drag2(0) - 2 * acc2(1) * drag2(1) - 2 * acc2(2) * drag2(2) + 2 * acc2(2) * G;
    const Scalar c_dis2 = drag2.squaredNorm() - 2 * drag2(2) * G + G * G - (max_acc_norm * max_acc_norm);
    const Scalar k2 = (-b_dis2 + sqrt(b_dis2 * b_dis2 - 4 * a_dis2 * c_dis2)) / (2 * a_dis2);
    Vector<3> thrust_acc_new_k2 = k2 * acc2 - drag2 - GVEC;
    Vector<3> max_acc_new12 = thrust_acc_new_k2 + drag2 + GVEC;

    for (int a=0;a<3;a++) {
      if (max_acc_new12(a) > 0 && max_acc_new12(a) < max_acc_new1(a)) {
        max_acc_new1(a) = max_acc_new12(a);
      } else if  (max_acc_new12(a) < 0 && max_acc_new12(a) > max_acc_new2(a)) {
        max_acc_new2(a) = max_acc_new12(a);
      }
    }

    const Scalar a_dis3 = acc3.squaredNorm();
    const Scalar b_dis3 = - 2 * acc3(0) * drag3(0) - 2 * acc3(1) * drag3(1) - 2 * acc3(2) * drag3(2) + 2 * acc3(2) * G;
    const Scalar c_dis3 = drag3.squaredNorm() - 2 * drag3(2) * G + G * G - (max_acc_norm * max_acc_norm);
    const Scalar k3 = (-b_dis3 + sqrt(b_dis3 * b_dis3 - 4 * a_dis3 * c_dis3)) / (2 * a_dis3);
    Vector<3> thrust_acc_new_k3 = k3 * acc3 - drag3 - GVEC;
    Vector<3> max_acc_new13 = thrust_acc_new_k3 + drag3 + GVEC;

    for (int a=0;a<3;a++) {
      if (max_acc_new13(a) > 0 && max_acc_new13(a) < max_acc_new1(a)) {
        max_acc_new1(a) = max_acc_new13(a);
      } else if  (max_acc_new13(a) < 0 && max_acc_new13(a) > max_acc_new2(a)) {
        max_acc_new2(a) = max_acc_new13(a);
      }
    }

    const Scalar a_dis4 = acc4.squaredNorm();
    const Scalar b_dis4 = - 2 * acc4(0) * drag4(0) - 2 * acc4(1) * drag4(1) - 2 * acc4(2) * drag4(2) + 2 * acc4(2) * G;
    const Scalar c_dis4 = drag4.squaredNorm() - 2 * drag4(2) * G + G * G - (max_acc_norm * max_acc_norm);
    const Scalar k4 = (-b_dis4 + sqrt(b_dis4 * b_dis4 - 4 * a_dis4 * c_dis4)) / (2 * a_dis4);
    Vector<3> thrust_acc_new_k4 = k4 * acc4 - drag4 - GVEC;
    Vector<3> max_acc_new14 = thrust_acc_new_k4 + drag4 + GVEC;

    for (int a=0;a<3;a++) {
      if (max_acc_new14(a) > 0 && max_acc_new14(a) < max_acc_new1(a)) {
        max_acc_new1(a) = max_acc_new14(a);
      } else if  (max_acc_new14(a) < 0 && max_acc_new14(a) > max_acc_new2(a)) {
        max_acc_new2(a) = max_acc_new14(a);
      }  

      if (max_acc_new1(a) > max_acc_norm)
        max_acc_new1(a) = max_per_axis_acc_vec(a);
      if (max_acc_new2(a) < -max_acc_norm)
        max_acc_new2(a) = min_per_axis_acc_vec(a);

    }

    if (debug) {
      std::cout << std::endl << "Iter: " << iter << std::endl;
      std::cout << "Thrust: " << largest_thrust << std::endl;
      std::cout << "Duration: " << duration << std::endl;

      std::cout << "max_acc_new11: " << max_acc_new11.transpose() << std::endl;
      std::cout << "max_acc_new12: " << max_acc_new12.transpose() << std::endl;
      std::cout << "max_acc_new13: " << max_acc_new13.transpose() << std::endl;
      std::cout << "max_acc_new14: " << max_acc_new14.transpose() << std::endl;

      std::cout << "max_acc_new1: " << max_acc_new1[0] << "," << max_acc_new1[1] << ","<< max_acc_new1[2] << std::endl;
      std::cout << "max_acc_new2: " << max_acc_new2[0] << "," << max_acc_new2[1] << ","<< max_acc_new2[2] << std::endl << std::endl;
    }

    pmm3d = PointMassTrajectory3D(from, to, max_acc_new1, max_acc_new2, max_velocity, true, true);

    max_velocity = pmm3d.max_velocity().normalized() * max_vel_norm;

    if (debug){
      std::cout << "T new: " << pmm3d.time() << " trajectory valid " << pmm3d.exists() << std::endl;
    }

    Vector<3> switch_times = pmm3d.switch_times();

    acc1 = pmm3d.acc_in_time(switch_times(0));
    acc2 = pmm3d.acc_in_time(switch_times(1));
    acc3 = pmm3d.acc_in_time(switch_times(2));
    acc4 = pmm3d.end_acc();

    drag1 = pmm3d.drag(acc1 - drag1 - GVEC, switch_times(0));
    drag2 = pmm3d.drag(acc2 - drag2 - GVEC, switch_times(1));
    drag3 = pmm3d.drag(acc3 - drag3 - GVEC, switch_times(2));
    drag4 = pmm3d.drag(acc4 - drag4 - GVEC, pmm3d.time());

    Scalar thrust1 = (acc1 - drag1 - GVEC).norm();
    Scalar thrust2 = (acc2 - drag2 - GVEC).norm();
    Scalar thrust3 = (acc3 - drag3 - GVEC).norm();
    Scalar thrust4 = (acc4 - drag4 - GVEC).norm();


    Scalar biggest_thrust = std::max(std::max(thrust1, thrust2), std::max(thrust3, thrust4));

    if (debug) {
      std::cout << std::endl << "switch times: " << switch_times.transpose() << std::endl;
      std::cout << "Thrust1: " << thrust1 << std::endl;
      std::cout << "Acc1: " << acc1.transpose() << " = " << acc1.norm() <<  std::endl;
      std::cout << "Drag1: " << drag1.transpose() << " = " << drag1.norm() <<  std::endl;
      std::cout << "Thrust2: " << thrust2 << std::endl;
      std::cout << "Acc2: " << acc2.transpose() << " = " << acc2.norm() <<  std::endl;
      std::cout << "Drag2: " << drag2.transpose() << " = " << drag2.norm() <<  std::endl;
      std::cout << "Thrust3: " << thrust3 << std::endl;
      std::cout << "Acc3: " << acc3.transpose() << " = " << acc3.norm() <<  std::endl;
      std::cout << "Drag3: " << drag3.transpose() << " = " << drag3.norm() <<  std::endl;
      std::cout << "Thrust4: " << thrust4 << std::endl;
      std::cout << "Acc4: " << acc4.transpose() << " = " << acc4.norm() <<  std::endl;
      std::cout << "Drag4: " << drag4.transpose() << " = " << drag4.norm() <<  std::endl;
    }

    if ((pmm3d.x_.t_(0) < PRECISION_PMM_VALUES || pmm3d.y_.t_(0) < PRECISION_PMM_VALUES || pmm3d.z_.t_(0) < PRECISION_PMM_VALUES) && 
        (pmm3d.x_.t_(2) < PRECISION_PMM_VALUES || pmm3d.y_.t_(2) < PRECISION_PMM_VALUES || pmm3d.z_.t_(2) < PRECISION_PMM_VALUES)) {
      if (debug) {
        std::cout << "ENDING EARLY: acceleration scaling makes no difference!" << std::endl;
      }
      break;
    }

    if (pmm3d.time() < duration && biggest_thrust <= max_acc_norm + PRECISION_PMM_VALUES && pmm3d.exists()) {
      x_ = pmm3d.x_;
      y_ = pmm3d.y_;
      z_ = pmm3d.z_;

      min_axis_trajectory_duration_ = pmm3d.min_axis_trajectory_duration_;

      duration = pmm3d.time();

      largest_thrust = biggest_thrust;
    }
  }

  if (debug){
    std::cout << "Iter: " << iter << std::endl;
    std::cout << "Thrust: " << largest_thrust << std::endl;
  }
}

/*
basic version with symmetric acc limits in all axis using sync segments
*/
// Creates symmetrical acc limits and calls overladed function
PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from,
                                             const QuadState &to,
                                             const Vector<3> max_acc,
                                             const Vector<3> max_vel,
                                             const bool equalize_time,
                                             const bool calc_gradient)
  : PointMassTrajectory3D(from, to, max_acc, -max_acc, max_vel, equalize_time, calc_gradient) {}

PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from,
                                             const QuadState &to,
                                             const Vector<3> max_acc1,
                                             const Vector<3> max_acc2,
                                             const Vector<3> max_vel,
                                             const bool equalize_time,
                                             const bool calc_gradient) {
  
  x_ = PMMTrajectory(from.p(0), from.v(0), to.p(0), to.v(0), max_acc1(0),
                     max_acc2(0), max_vel(0), 0, 0.0, false, calc_gradient, false);
  y_ = PMMTrajectory(from.p(1), from.v(1), to.p(1), to.v(1), max_acc1(1),
                     max_acc2(1), max_vel(1), 1, 0.0, false, calc_gradient, false);
  z_ = PMMTrajectory(from.p(2), from.v(2), to.p(2), to.v(2), max_acc1(2),
                     max_acc2(2), max_vel(2), 2, 0.0, false, calc_gradient, false);

  // keep track of min trajectory duration for each axis
  min_axis_trajectory_duration_[0] = x_.time();
  min_axis_trajectory_duration_[1] = y_.time();
  min_axis_trajectory_duration_[2] = z_.time();

  // std::cout << "min_axis_trajectory_duration_ " << min_axis_trajectory_duration_.transpose() << std::endl;

  if (equalize_time && x_.exists_ && y_.exists_ && z_.exists_) {
    // compute sync trajectories according to bang-bang or bang-singular-bang approach
    Scalar tr_time = time();
    for(int i = 0; i<3; i++){
      auto axis_tr = get_axis_trajectory(i);
      if (fabs(axis_tr.time() - tr_time) > PRECISION_PMM_VALUES){
        // recomputation needed
        Scalar res = axis_tr.computeSyncTrajectory(tr_time, max_acc1(i), max_acc2(i), max_vel(i));
        if (fabs(res) > PRECISION_PMM_VALUES){
          if (res > tr_time){
            // set new trajectory time and recompute all axis
            tr_time = res;
            min_axis_trajectory_duration_[i] = tr_time;
            set_axis_trajectory(i, axis_tr);
            i = -1;
            continue;
          } else {
            axis_tr.exists_=false;
          }
        }
        set_axis_trajectory(i, axis_tr);
      }
    }
  }
}

Vector<3> PointMassTrajectory3D::acc_in_time(const Scalar time_in_tr) const {
  Vector<3> acc;
  acc[0] = x_.acc_in_time(time_in_tr);
  acc[1] = y_.acc_in_time(time_in_tr);
  acc[2] = z_.acc_in_time(time_in_tr);
  return acc;
}

Vector<3> PointMassTrajectory3D::vel_in_time(const Scalar time_in_tr) const {
  Vector<3> vel;
  vel[0] = x_.vel_in_time(time_in_tr);
  vel[1] = y_.vel_in_time(time_in_tr);
  vel[2] = z_.vel_in_time(time_in_tr);
  return vel;
}

QuadState PointMassTrajectory3D::state_in_time(const Scalar time_in_tr) const {
  QuadState ds;
  ds.setZero();
  Vector<3> x = x_.state_in_time(time_in_tr);
  Vector<3> y = y_.state_in_time(time_in_tr);
  Vector<3> z = z_.state_in_time(time_in_tr);

  ds.p(0) = x(0);
  ds.p(1) = y(0);
  ds.p(2) = z(0);

  ds.v(0) = x(1);
  ds.v(1) = y(1);
  ds.v(2) = z(1);

  ds.a(0) = x(2);
  ds.a(1) = y(2);
  ds.a(2) = z(2);
  return ds;
}

QuadState PointMassTrajectory3D::get_start_state() const {
  return state_in_time(0);
}

QuadState PointMassTrajectory3D::get_end_state() const {
  return state_in_time(time());
}

Vector<3> PointMassTrajectory3D::start_acc() const {
  return Vector<3>(x_.a_(0), y_.a_(0), z_.a_(0));
}

Vector<3> PointMassTrajectory3D::end_acc() const {
  return Vector<3>(x_.a_(1), y_.a_(1), z_.a_(1));
}

Vector<3> PointMassTrajectory3D::switch_times() const {
  Vector<3> t = Vector<3>(x_.t_(0), y_.t_(0), z_.t_(0)) - Vector<3>::Constant(TIME_PRECISION);
  for (int a=0;a<3;a++) {
    if (t(a) < TIME_PRECISION) {
      t(a) = TIME_PRECISION;
    }
  }
  return t;
}

Vector<3> PointMassTrajectory3D::zero_time() const {
  // return Vector<3>(0, 0, 0);
  return Vector<3>(x_.t_(1), y_.t_(1), z_.t_(1)) / time();
}

Vector<3> PointMassTrajectory3D::max_velocity() const {
  Scalar v_x = std::max(std::max(fabs(x_.v_(0)), fabs(x_.v_(1))), fabs(x_.v_(2)));
  Scalar v_y = std::max(std::max(fabs(y_.v_(0)), fabs(y_.v_(1))), fabs(y_.v_(2)));
  Scalar v_z = std::max(std::max(fabs(z_.v_(0)), fabs(z_.v_(1))), fabs(z_.v_(2)));
  return Vector<3>(v_x, v_y, v_z);
}

Vector<3> PointMassTrajectory3D::drag(Vector<3> acc, const Scalar time_in_tr) const {

  Vector<3> thrust_acceleration = acc;
  double heading = 0;
  double acc_norm = thrust_acceleration.norm();
  Vector<3> z_body = thrust_acceleration.normalized();
  Vector<3> x_C = Vector<3>(cos(heading), sin(heading), 0);
  Vector<3> y_body = (z_body.cross(x_C)).normalized();
  Vector<3> x_body = y_body.cross(z_body);
  Eigen::Matrix3d rotation_matrix;
  rotation_matrix << x_body, y_body, z_body;
  Eigen::Quaterniond quat = Eigen::Quaterniond(rotation_matrix);

  Vector<3> body_vel = quat.inverse() * vel_in_time(time_in_tr);
  Vector<3> drag_acc_body = D_coeffs * body_vel;
  Vector<3> drag_acc_world = quat * drag_acc_body;
  return drag_acc_world;
}

void PointMassTrajectory3D::set_axis_trajectory(const int i,
                                                const PMMTrajectory tr) {
  switch (i) {
    case 0:
      x_ = tr;
      break;
    case 1:
      y_ = tr;
      break;
    case 2:
      z_ = tr;
      break;
    default:
      std::cout << "bad axis index " << i << std::endl;
  }
}
PMMTrajectory &PointMassTrajectory3D::get_axis_trajectory(const int i) {
  switch (i) {
    case 0:
      return x_;
    case 1:
      return y_;
    case 2:
      return z_;
    default:
      std::cout << "bad axis index " << i << std::endl;
      return x_;
  }
}

Scalar PointMassTrajectory3D::get_axis_switch_time(const int i) const {
  switch (i) {
    case 0:
      return x_.t_(0);
    case 1:
      return y_.t_(0);
    case 2:
      return z_.t_(0);
    default:
      exit(1);
  }
}

PointMassTrajectory3D PointMassTrajectory3D::operator=(const PointMassTrajectory3D &tr){
    this->copy_trajectory(tr);
    return *this;
}

void PointMassTrajectory3D::copy_trajectory(const PointMassTrajectory3D &tr){

    this->x_ = tr.x_;
    this->y_ = tr.y_;
    this->z_ = tr.z_;
    this->min_axis_trajectory_duration_ = tr.min_axis_trajectory_duration_;
    
}


std::ostream &operator<<(std::ostream &o, const PointMassTrajectory3D &t) {
  o << "pmm3d: t:" << t.time() << ";exists:" << t.exists();
  o << "\n\tx: " << t.x_;
  o << "\n\ty: " << t.y_;
  o << "\n\tz: " << t.z_;
  return o;
}
} // namespace pmm