//
// Created by Hao-Xin on 2022/10/1.
//

#ifndef HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCK_DOF_H
#define HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCK_DOF_H

#include <array>
#include <cmath>
#include "../MonteCarlo_src/local_dof.h"

template<size_t q>
class ClockDOF {
  using LocalDOFT = LocalDOF<2>;
 public:
  ClockDOF(void) = default;
  ClockDOF(const size_t state) : state_(state) {}
  ClockDOF(const ClockDOF &clock_dof) : state_(clock_dof.state_) {}

  ClockDOF &operator=(const ClockDOF &rhs) {
    state_ = rhs.state_;
    return *this;
  }

  /// inner product
  double operator*(const ClockDOF &rhs) const {
//    size_t delta_s = (state_ > rhs.state_) ? (state_ - rhs.state_) : (state_ + q - rhs.state_);
//    return sx_set[delta_s];

    return sx_set[state_] * sx_set[rhs.state_] + sy_set[state_] * sy_set[rhs.state_];
  }

  LocalDOFT operator+(const ClockDOF &rhs) const {
    LocalDOFT clock({sx_set[state_], sy_set[state_]});
    LocalDOFT clock_rhs({sx_set[rhs.state_], sy_set[rhs.state_]});
    return clock + clock_rhs;
  }

  LocalDOFT operator+(const LocalDOFT &rhs) const {
    LocalDOFT clock({sx_set[state_], sy_set[state_]});
    return clock + rhs;
  }

  double Cos2Theta() const {
    return sx_set[state_] * sx_set[state_] - sy_set[state_] * sy_set[state_];
  }

  double Sin2Theta() const {
    return 2 * sx_set[state_] * sy_set[state_];
  }
/// Flip the spin_ according to vector axis.
/// note this is different from the Flip in LocalDOF
  void Flip2(const ClockDOF axis) {
    state_ = (2 * axis.state_ + q - state_) % q;
  }

  size_t GetCoor() const {
    return state_;
  }

  void Random() {
    state_ = u_int(random_engine);
  }

  void Show() {
    std::cout << state_ << std::endl;
  }

 private:
  static const std::array<double, q> angles;
  static const std::array<double, q> sx_set;
  static const std::array<double, q> sy_set;
  static std::uniform_int_distribution<int> u_int;
  size_t state_; // 0 <= state < q




};


template<> std::uniform_int_distribution<int> ClockDOF<6>::u_int = std::uniform_int_distribution<int>(0, 6-1);


template<size_t q>
std::array<double, q> InitializeAngleSet() {
  std::array<double, q> angles;
  for (size_t i = 0; i < q; i++) {
    angles[i] = 2 * M_PI * i / q;
  }
  return angles;
}

template<size_t q>
std::array<double, q> InitializeSxSet() {
  std::array<double, q> sx_set;
  for (size_t i = 0; i < q; i++) {
    sx_set[i] = cos(2 * M_PI * i / q);
  }
  return sx_set;
}

template<size_t q>
std::array<double, q> InitializeSySet() {
  std::array<double, q> sy_set;
  for (size_t i = 0; i < q; i++) {
    sy_set[i] = sin(2 * M_PI * i / q);
  }
  return sy_set;
}

template<> const std::array<double, 6> ClockDOF<6>::angles = InitializeAngleSet<6>();
template<> const std::array<double, 6> ClockDOF<6>::sx_set = InitializeSxSet<6>();
template<> const std::array<double, 6> ClockDOF<6>::sy_set = InitializeSySet<6>();

#endif //HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCK_DOF_H
