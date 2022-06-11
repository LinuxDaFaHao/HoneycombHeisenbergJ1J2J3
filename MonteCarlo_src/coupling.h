//
// Created by Hao-Xin on 2022/5/24.
//

#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_COUPLING_H
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_COUPLING_H


#include <array>
#include "local_dof.h"
#include "./const.h"


template <size_t DimDof>
class CouplingStructure {
  using LocalDOFT = LocalDOF<DimDof>;
 public:
  CouplingStructure(void) = default;
  CouplingStructure(const CouplingStructure& coupling_structure) {
    for(size_t i = 0; i < DimDof; i++) {
      for (size_t j = 0; j < DimDof; j++) {
        coupling_const[i][j] = coupling_structure.coupling_const[i][j];
      }
    }
  }

  CouplingStructure(const double coupling_const_rhs[DimDof][DimDof]) {
    for (size_t i = 0; i < DimDof; i++) {
      for (size_t j = 0; j < DimDof; j++) {
        coupling_const[i][j] = coupling_const_rhs[i][j];
      }
    }
  }

  CouplingStructure &operator=(const CouplingStructure &rhs) {
    for(size_t i = 0; i < DimDof; i++) {
      for (size_t j = 0; j < DimDof; j++) {
        coupling_const[i][j] = rhs.coupling_const[i][j];
      }
    }
    return *this;
  }

  bool operator==(const CouplingStructure &rhs) {
    for(size_t i = 0; i < DimDof; i++) {
      for (size_t j = 0; j < DimDof; j++) {
        if(coupling_const[i][j] != rhs.coupling_const[i][j]){
          return false;
        }
      }
    }
    return true;
  }


  CouplingStructure &operator*=(const double scalar) {
    for(size_t i = 0; i < DimDof; i++) {
      for (size_t j = 0; j < DimDof; j++) {
        coupling_const[i][j] *= scalar;
      }
    }
    return *this;
  }

  CouplingStructure operator*(const double scalar) const {
    CouplingStructure res(*this);
    res *= scalar;
    return res;
  }

  LocalDOFT operator*(const LocalDOFT& n) const {
    std::array<double, DimDof> result;
    for(size_t i = 0; i < DimDof; i++) {
      double sum(0.0);
      for(size_t j = 0; j < DimDof; j++) {
        sum += coupling_const[i][j] * n.GetCoor(j);
      }
      result[i] = sum;
    }
    return LocalDOF(result);
  }

  bool IsIsometry() const {
    double coupling00 = coupling_const[0][0];
    for(size_t i = 0; i < DimDof; i++) {
      for (size_t j = 0; j < DimDof; j++) {
        if(i == j) {
          if(coupling_const[i][j] != coupling00) {
            return false;
          }
        } else {
          if(coupling_const[i][j] != 0.0) {
            return false;
          }
        }
      }
    }
    return true;
  }

  //coupling constant, Column first (the second index traverse columns of matrix).
  double coupling_const[DimDof][DimDof];
};



#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_COUPLING_H