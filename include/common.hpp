/**
 * @file common.hpp
 * @author Krystof Teissing (teisskry@gmail.com)
 * @version 0.1
 * @date 2023-05-25
 * 
 */

#pragma once

#include <Eigen/Dense>
#include "yaml-cpp/yaml.h"

namespace pmm {

// TYPES -----------------------------------------------------------------/
// Define the scalar type used.
using Scalar = double;
static constexpr Scalar INF = std::numeric_limits<Scalar>::infinity();

// Define `Dynamic` matrix size.
static constexpr int Dynamic = Eigen::Dynamic;

// Using shorthand for `Matrix<rows, cols>` with scalar type.
template<int rows = Dynamic, int cols = rows>
using Matrix = Eigen::Matrix<Scalar, rows, cols>;

// Using shorthand for `Vector<rows>` with scalar type.
template<int rows = Dynamic>
using Vector = Matrix<rows, 1>;

// Using shorthand for `Array<rows, cols>` with scalar type.
template<int rows = Dynamic, int cols = rows>
using Array = Eigen::Array<Scalar, rows, cols>;

template<int rows = Dynamic>
using ArrayVector = Array<rows, 1>;

// STRUCTS ----------------------------------------------------------------/
struct QuadState {
  Vector<3> p;
  Vector<3> v;
  Vector<3> a;
  void setZero(){
    p.setZero();
    v.setZero();
  }
};

// GRAVITY ----------------------------------------------------------------/
// Gravity Value [m/s^2]
extern Scalar G;

// Gravity Vector [m/s^2]
extern Vector<3> GVEC;

// DRAG -------------------------------------------------------------------/
extern Eigen::Matrix3d D_coeffs;


// MISC -------------------------------------------------------------------/
void operator>>(const YAML::Node& node, bool& value);
void operator>>(const YAML::Node& node, int& value);
void operator>>(const YAML::Node& node, Scalar& value);
void operator>>(const YAML::Node& node, Vector<3>& value);
void operator>>(const YAML::Node& node, std::string& value);

void redistribute(Eigen::Vector3d& vec);

template<typename DataType>
bool parseArrayParam(const YAML::Node& config, std::string param,
                     std::vector<DataType>& array) {
  if (config[param.c_str()]) {
    YAML::Node vector = config[param.c_str()];
    for (YAML::const_iterator ti = vector.begin(); ti != vector.end(); ++ti) {
      const YAML::Node& datayaml = *ti;
      DataType from_yaml;
      datayaml >> from_yaml;
      array.push_back(from_yaml);
    }
  } else {
    return false;
  }
  return true;
}
} // namespace pmm