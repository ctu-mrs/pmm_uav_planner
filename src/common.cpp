/**
 * @file common.cpp
 * @author Krystof Teissing (teisskry@gmail.com)
 * @version 0.1
 * @date 2023-05-25
 * 
 */

#include "common.hpp"

namespace pmm {

// GRAVITY ----------------------------------------------------------------/
// Gravity Value [m/s^2]
Scalar G = 9.8066;

// Gravity Vector [m/s^2]
Vector<3> GVEC = {0, 0, -G};

// DRAG -------------------------------------------------------------------/
Eigen::Matrix3d D_coeffs = -(Eigen::Matrix3d()<< 0.0, 0, 0, 0, 0.0, 0, 0, 0, 0.0).finished();

void operator>>(const YAML::Node& node, bool& value) {
  value = node.as<bool>();
}

void operator>>(const YAML::Node& node, int& value) {
  value = node.as<int>();
}

void operator>>(const YAML::Node& node, Scalar& value) {
  assert(node.IsScalar());
  value = node.as<Scalar>();
}
void operator>>(const YAML::Node& node, std::string& value) {
  value = node.as<std::string>();
}

void operator>>(const YAML::Node& node, Vector<3>& value) {
    std::vector<Scalar> tmp;
    assert(node.IsSequence());
    tmp = node.as<std::vector<Scalar>>();
    value[0] = tmp[0];
    value[1] = tmp[1];
    value[2] = tmp[2];
}

} // namespace pmm