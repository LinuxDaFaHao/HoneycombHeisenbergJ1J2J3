//
// Created by Hao-Xin on 2022/5/22.
//

/**
@file local_dof.h
@brief The local degree of freedom
*/
#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_LOCAL_DOF_H
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_LOCAL_DOF_H

#include <array>
#include <assert.h>
#include <random>
#include <iostream>
#include "./const.h"

std::default_random_engine random_engine;

template <size_t DimDof>
class LocalDOF {
 public:
  LocalDOF(void) = default;
  LocalDOF(const LocalDOF &local_dof) : spin_(local_dof.spin_) {}
  LocalDOF(const std::array<double, DimDof> &spin) : spin_(spin) {}
  LocalDOF &operator=(const LocalDOF &rhs) {
    spin_ = rhs.spin_;
    return *this;
  }

  LocalDOF operator*(const double scalar) const {
    LocalDOF result;
    for(size_t i = 0; i < DimDof; i++) {
      result.spin_[i] = spin_[i] * scalar;
    }
    return result;
  }

  LocalDOF &operator*=(const double scalar) {
    for(size_t i = 0; i < DimDof; i++) {
      spin_[i] = spin_[i] * scalar;
    }
    return (*this);
  }

  ///  inner product
  double operator*(const LocalDOF &rhs) const {
    double inner_prod(0.0);
    for(size_t i = 0; i < DimDof; i++) {
      inner_prod += spin_[i] * rhs.spin_[i];
    }
    return inner_prod;
  }

  LocalDOF operator+(const LocalDOF &rhs) const {
    LocalDOF res;
    for(size_t i = 0; i < DimDof; i++) {
      res.spin_[i] = spin_[i] + rhs.spin_[i];
    }
    return res;
  }

  LocalDOF& operator+=(const LocalDOF &rhs) {
    for(size_t i = 0; i < DimDof; i++) {
      spin_[i] += rhs.spin_[i];
    }
    return *this;
  }

  LocalDOF operator-(void) const {
    LocalDOF res;
    for(size_t i = 0; i < DimDof; i++) {
      res.spin_[i] = -spin_[i];
    }
    return res;
  }

  LocalDOF &operator-=(const LocalDOF &rhs) {
    (*this) += (-rhs);
    return *this;
  }

  LocalDOF DecompPara(const LocalDOF& n) const {
    if(DimDof == 1) { //for performance
      LocalDOF res(*this);
      return res;
    } else {
      return n * ((*this) * n);
    }
  }

  LocalDOF DecompPerpen(const LocalDOF& n) const {
    if(DimDof == 1) {
      return LocalDOF<1>({0.0});
    } else {
      return *this - this->DecompPara(n);
    }
  }

  std::pair<LocalDOF, LocalDOF> DecompSpin(const LocalDOF& n) const {
    std::pair<LocalDOF, LocalDOF> res;
    res.first = DecompPara(n);
    res.second = *this + (- res.first);
    return res;
  }

  int Sign(const LocalDOF& n) const {
    if((*this) * n > 0) { return 1;}
    else { return -1;}
  }

  /// Flip the spin_ according to the plain perpendicular to vector n.
  void Flip(const LocalDOF& n) {
    assert( n.IsNormlized() && this->IsNormlized());
    if(DimDof == 1) {
      spin_[0] = -spin_[0];
    } else{
      LocalDOF s_para = DecompPara(n);
      (*this) += s_para * (-2.0);
      assert(this->IsNormlized());
      this -> Normalize();
    }
  }
  

  double GetCoor(const size_t i) const {
    assert(i < DimDof);
    return spin_[i];
  }

  double Norm2() const {
    double sum2(0.0);
    for(size_t i = 0; i < DimDof; i++) {
      sum2 += spin_[i] * spin_[i];
    }
    return std::sqrt(sum2);
  }

  double Normalize() {
    double norm2 = Norm2();
    for(size_t i = 0; i < DimDof; i++) {
      spin_[i] /= norm2;
    }
    return norm2;
  }

  bool IsNormlized() const {
    return (this->Norm2() - 1.0) < 1e-13;
  }
  /// Generate a random local_dof with length 1.
  /// TODO: more effective algorithm, e.g. Box–Muller transform, see in wikipedia
  void Random(){
    std::uniform_real_distribution<double> u(-1, 1);
    do {
      for(size_t i = 0; i < DimDof; i++) {
        spin_[i] = u(random_engine);
      }
    }while(Norm2() > 1);
    Normalize();
  }

  void Zero() {
    for(size_t i = 0; i < DimDof; i++) {
      spin_[i] = 0;
    }
  }

  void Show() {
    std::cout << "[ ";
    for(size_t i = 0; i < DimDof - 1; i++) {
      std::cout << spin_[i] << ", ";
    }
    std::cout << spin_[DimDof-1] << "]" << std::endl;
  }

 private :
  std::array<double, DimDof> spin_;

  template <size_t DimDof2>
  friend double Project2DInnerProduct(
      const LocalDOF<DimDof2> &s1,
      const LocalDOF<DimDof2> &s2
  );

  template <size_t DimDof2>
  friend double ZCompOfCrossProduct(
      const LocalDOF<DimDof2> &s1,
      const LocalDOF<DimDof2> &s2
  );
};

/// sin(theta1 - theta2)
template <size_t DimDof>
double SinDiff(
    const LocalDOF<DimDof> &theta1,
    const LocalDOF<DimDof> &theta2
    ) {
  if(DimDof == 2) {
    double res;
    res = theta2.GetCoor(0) * theta1.GetCoor(1) - theta1.GetCoor(0) * theta2.GetCoor(1);
#ifndef NDEBUG
    double angle1 = ::atan2(theta1.GetCoor(1), theta1.GetCoor(0));
    double angle2 = ::atan2(theta2.GetCoor(1), theta2.GetCoor(0));
    double res2 = ::sin(angle1 - angle2);
    assert( fabs(res - res2) < 1e-14);
#endif
    return res;
  } else {
    std::cout << "I'm not sure how to deal with sin(theta1 - theta2) with DimDof != 2" <<std::endl;
    return 0.0;
  }
}

template <size_t DimDof>
double Project2DInnerProduct(
    const LocalDOF<DimDof> &s1,
    const LocalDOF<DimDof> &s2
    ) {
  assert(DimDof >=2);
  return s1.spin_[0] * s2.spin_[0] + s1.spin_[1] * s2.spin_[1];
}

template <size_t DimDof>
double ZCompOfCrossProduct(
    const LocalDOF<DimDof> &s1,
    const LocalDOF<DimDof> &s2
    ) {
  assert(DimDof >=2 );
  return s1.spin_[0] * s2.spin_[1] - s1.spin_[1] * s2.spin_[0];
}


#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_LOCAL_DOF_H
