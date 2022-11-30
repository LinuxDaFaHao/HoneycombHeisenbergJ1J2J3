//
// Created by Hao-Xin on 2022/5/24.
//

#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICELINK_H
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICELINK_H

#include <vector>
#include "./coupling.h"

/// define the Hamiltonian
template<size_t DimDof, size_t NumOfCouplingType>
class LatticeLink {
  using CST = CouplingStructure<DimDof>;
  using OrientedBond = std::pair<size_t, size_t>;
 public:
  LatticeLink(void) = default;

  LatticeLink(const std::array<CST, NumOfCouplingType> &coupling_structures) :
      coupling_structures(coupling_structures), area_(0.0) {}

  LatticeLink(const std::array<CST, NumOfCouplingType> &coupling_structures,
              const double area) :
      coupling_structures(coupling_structures), area_(area) {}

  LatticeLink(const LatticeLink &lattice_link) :
      coupling_structures(lattice_link.coupling_structures),
      links(lattice_link.links),
      area_(lattice_link.area_) {}

  ~LatticeLink() = default;

  std::vector<std::pair<size_t, size_t>> GetBondLinkedTo(const size_t site, const size_t InterType) {
    std::vector<std::pair<size_t, size_t>> bonds;
    const auto &site2_set = links[InterType][site];
    for (auto &site2: site2_set) {
      bonds.push_back(std::make_pair(site, site2));
    }
    return bonds;
  }

  const std::vector<size_t> &GetSitesLinkedTo(const size_t site, const size_t InterType) const {
    return links[InterType][site];
  }

  size_t GetNumofInteraction(const size_t inter_type) const {
    size_t num = 0;
    for (auto &site2_set: links[inter_type]) {
      num += site2_set.size();
    }
    if (inter_type > 0) {
      num = num / 2;
    }
    return num;
  }

  double GetLatticeArea() const {
    return area_;
  }

  virtual size_t Distance(const size_t i, const size_t j) const { return 0; };

  ///0, single ion anisotropy D;
  ///1, J1; 2, J2; 3, J3;
  std::array<CST, NumOfCouplingType> coupling_structures;
  std::array<std::vector<std::vector<size_t>>, NumOfCouplingType> links;
  /// links[i][j] is a std::vector<size_t>, which is the set of under interaction Ji(D) the sites link to site j.
  std::array<std::vector<OrientedBond>, NumOfCouplingType> oriented_links;
  std::array<std::vector<double>, NumOfCouplingType>
      proj2x_link_length; // corresponding to the elements of oriented_links
  std::array<std::vector<double>, NumOfCouplingType>
      proj2y_link_length; // corresponding to the elements of oriented_links

  size_t lx;
  size_t ly;
  size_t n;
 private:
  const double area_; // set distance of NN atoms as 1
};

/// all link
template<size_t DimDof, size_t NumOfCouplingType>
class SVWLatticeLink : public LatticeLink<DimDof, NumOfCouplingType> {
  using CST = CouplingStructure<DimDof>;
 public:
  SVWLatticeLink(const size_t N, const std::array<CST, NumOfCouplingType> &coupling_structures) :
      LatticeLink<DimDof, NumOfCouplingType>(coupling_structures) {
    this->n = N;
    for (size_t i = 0; i < NumOfCouplingType; i++) {
      this->links[i].resize(N);
    }
    /// SIA
    for (size_t i = 0; i < N; i++) {
      this->links[0][i] = {i};
    }
    /// Interaction
    for (size_t i = 0; i < N; i++) {
      this->links[1][i] = {};
      this->links[1][i].reserve(N - 1);
      for (size_t j = 0; j < N; j++) {
        if (i != j)
          this->links[1][i].push_back(j);
      }
    }
  }

  ~SVWLatticeLink() = default;

  size_t Distance(const size_t i, const size_t j) const override;
};

template<size_t DimDof, size_t NumOfCouplingType>
size_t SVWLatticeLink<DimDof, NumOfCouplingType>::Distance(const size_t i, const size_t j) const {
  return (i > j) ? (i - j) : (j - i);
}

template<size_t DimDof, size_t NumOfCouplingType>
class SquareTorusLatticeLink : public LatticeLink<DimDof, NumOfCouplingType> {
  using CST = CouplingStructure<DimDof>;
 public:
  SquareTorusLatticeLink(const size_t Lx, const size_t Ly,
                         const std::array<CST, NumOfCouplingType> &coupling_structures) :
      LatticeLink<DimDof, NumOfCouplingType>(coupling_structures, Lx * Ly) {
    this->lx = Lx;
    this->ly = Ly;
    const size_t N = Lx * Ly;
    this->n = N;
    for (size_t i = 0; i < NumOfCouplingType; i++) {
      this->links[i].resize(N);
      if (i > 0) this->oriented_links[i].reserve(2 * N);
    }
    /// SIA
    for (size_t i = 0; i < N; i++) {
      this->links[0][i] = {i};
    }
    /// NN
    if (!coupling_structures[1].IsZero()) {
      for (size_t i = 0; i < N; i++) {
        size_t x = i % (Lx);
        size_t y = i / (Lx);
        size_t Tx = (x + 1) % (Lx) + y * (Lx);
        size_t Tmx = (x + Lx - 1) % (Lx) + y * (Lx);
        size_t Ty = x + ((y + 1) % Ly) * (Lx);
        size_t Tmy = x + ((y + Ly - 1) % Ly) * (Lx);
        this->links[1][i] = {Tx, Ty, Tmx, Tmy};

        this->oriented_links[1].push_back(std::make_pair(i, Tx));
        this->proj2x_link_length[1].push_back(1.0);
        this->proj2y_link_length[1].push_back(0.0);
        this->oriented_links[1].push_back(std::make_pair(i, Ty));
        this->proj2x_link_length[1].push_back(0.0);
        this->proj2y_link_length[1].push_back(1.0);
      }
    }
    /// NNN
    if (NumOfCouplingType > 2 && !coupling_structures[2].IsZero()) {
      for (size_t i = 0; i < N; i++) {
        size_t x = i % (Lx);
        size_t y = i / (Lx);
        size_t Txy = (x + 1) % (Lx) + ((y + 1) % Ly) * (Lx);
        size_t Txmy = (x + 1) % (Lx) + ((y + Ly - 1) % Ly) * (Lx);
        size_t Tmxy = (x + Lx - 1) % (Lx) + ((y + 1) % Ly) * (Lx);
        size_t Tmxmy = (x + Lx - 1) % (Lx) + ((y + Ly - 1) % Ly) * (Lx);

        this->links[2][i] = {Txy, Txmy, Tmxy, Tmxmy}; //sort?

        this->oriented_links[2].push_back(std::make_pair(i, Txy));
        this->proj2x_link_length[2].push_back(1.0);
        this->proj2y_link_length[2].push_back(1.0);

        this->oriented_links[2].push_back(std::make_pair(i, Txmy));
        this->proj2x_link_length[2].push_back(1.0);
        this->proj2y_link_length[2].push_back(-1.0);
      }
    }


    /// TNN
    if (NumOfCouplingType > 3 && !coupling_structures[3].IsZero()) {
      for (size_t i = 0; i < N; i++) {
        size_t x = i % (Lx);
        size_t y = i / (Lx);
        size_t Tx2 = (x + 2) % (Lx) + y * (Lx);
        size_t Tmx2 = (x + Lx - 2) % (Lx) + y * (Lx);
        size_t Ty2 = x + ((y + 2) % Ly) * (Lx);
        size_t Tmy2 = x + ((y + Ly - 2) % Ly) * (Lx);
        this->links[3][i] = {Tx2, Ty2, Tmx2, Tmy2};

        this->oriented_links[3].push_back(std::make_pair(i, Tx2));
        this->oriented_links[3].push_back(std::make_pair(i, Ty2));
        this->proj2x_link_length[3].push_back(2.0);
        this->proj2y_link_length[3].push_back(2.0);
      }
    }
  }

  ~SquareTorusLatticeLink() = default;
  size_t Distance(const size_t i, const size_t j) const override;
};

template<size_t DimDof, size_t NumOfCouplingType>
size_t SquareTorusLatticeLink<DimDof, NumOfCouplingType>::Distance(const size_t i, const size_t j) const {
  size_t x_i = i % (this->lx);
  size_t y_i = i / (this->lx);
  size_t x_j = j % (this->lx);
  size_t y_j = j / (this->lx);
  size_t delta_x = (x_i > x_j) ? (x_i - x_j) : (x_j - x_i);
  size_t delta_y = (y_i > y_j) ? (y_i - y_j) : (y_j - y_i);
  return delta_x + delta_y;
}

/// YC-torus
template<size_t DimDof, size_t NumOfCouplingType>
class HoneyCombTorusLatticeLink : public LatticeLink<DimDof, NumOfCouplingType> {
  using CST = CouplingStructure<DimDof>;
 public:
  HoneyCombTorusLatticeLink(const size_t Lx, const size_t Ly,
                            const std::array<CST, NumOfCouplingType> &coupling_structures) :
      LatticeLink<DimDof, NumOfCouplingType>(coupling_structures, Lx * Ly * 3 * sqrt(3) / 2),
      Lx_(Lx), Ly_(Ly), N_(2 * Lx * Ly) {
    const size_t N = Lx * Ly * 2;
    this->lx = Lx;
    this->ly = Ly;
    this->n = N;
    for (size_t i = 0; i < NumOfCouplingType; i++) {
      this->links[i].resize(N);
      if (i > 0) this->oriented_links[i].reserve(3 * N); //reserve more for J1, J3
    }
    /// for single ion anisotropy D, on-site
    for (size_t i = 0; i < N; i++) {
      this->links[0][i] = {i};
    }
    /// for nearest neighbor (NN) J1
    for (size_t i = 0; i < N; i++) {
      size_t x = i % (2 * Lx);
      size_t y = i / (2 * Lx);
      size_t Tx = (x + 1) % (2 * Lx) + y * (2 * Lx);
      size_t Tmx = (x + 2 * Lx - 1) % (2 * Lx) + y * (2 * Lx);
      size_t Ty = x + ((y + 1) % Ly) * (2 * Lx);
      size_t Tmy = x + ((y + Ly - 1) % Ly) * (2 * Lx);
      if ((x + y) % 2 == 0) {
        this->links[1][i] = {Tmx, Tx, Ty}; //sort?
        this->oriented_links[1].push_back(std::make_pair(i, Ty));
        this->proj2x_link_length[1].push_back(0.0);
        this->proj2y_link_length[1].push_back(1.0);
        this->oriented_links[1].push_back(std::make_pair(i, Tx));
        this->proj2x_link_length[1].push_back(sqrt(3.0) / 2.0);
        this->proj2y_link_length[1].push_back(-1.0 / 2.0);
      } else {
        this->links[1][i] = {Tmx, Tx, Tmy};
        this->oriented_links[1].push_back(std::make_pair(i, Tx));
        this->proj2x_link_length[1].push_back(sqrt(3.0) / 2.0);
        this->proj2y_link_length[1].push_back(1.0 / 2.0);
      }
    }

    /// for NNN J2
    if (NumOfCouplingType > 2 && !coupling_structures[2].IsZero()) {
      for (size_t i = 0; i < N; i++) {
        size_t x = i % (2 * Lx);
        size_t y = i / (2 * Lx);
        size_t Tx2 = (x + 2) % (2 * Lx) + y * (2 * Lx);
        size_t Tmx2 = (x + 2 * Lx - 2) % (2 * Lx) + y * (2 * Lx);
        size_t Txy = (x + 1) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
        size_t Txmy = (x + 1) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);
        size_t Tmxy = (x + 2 * Lx - 1) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
        size_t Tmxmy = (x + 2 * Lx - 1) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);

        this->links[2][i] = {Tx2, Tmx2, Txy, Txmy, Tmxy, Tmxmy}; //sort?

        this->oriented_links[2].push_back(std::make_pair(i, Tx2));
        this->proj2x_link_length[2].push_back(sqrt(3.0));
        this->proj2y_link_length[2].push_back(0.0);

        this->oriented_links[2].push_back(std::make_pair(i, Txy));
        this->proj2x_link_length[2].push_back(sqrt(3.0) / 2.0);
        this->proj2y_link_length[2].push_back(3.0 / 2.0);

        this->oriented_links[2].push_back(std::make_pair(i, Txmy));
        this->proj2x_link_length[2].push_back(sqrt(3.0) / 2.0);
        this->proj2y_link_length[2].push_back(-3.0 / 2.0);
      }
    }

    /// for TNN J3
    if (NumOfCouplingType > 3 && !coupling_structures[3].IsZero()) {
      for (size_t i = 0; i < N; i++) {
        size_t x = i % (2 * Lx);
        size_t y = i / (2 * Lx);
        if ((x + y) % 2 == 0) {
          size_t Tmy = x + ((y + Ly - 1) % Ly) * (2 * Lx);
          size_t Tx2y = (x + 2) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
          size_t Tmx2y = (x + 2 * Lx - 2) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
          this->links[3][i] = {Tmy, Tx2y, Tmx2y};

          this->oriented_links[3].push_back(std::make_pair(i, Tx2y));
          this->proj2x_link_length[3].push_back(sqrt(3.0));
          this->proj2y_link_length[3].push_back(1.0);
        } else {
          size_t Ty = x + ((y + 1) % Ly) * (2 * Lx);
          size_t Tx2my = (x + 2) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);
          size_t Tmx2my = (x + 2 * Lx - 2) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);
          this->links[3][i] = {Ty, Tx2my, Tmx2my};
          this->oriented_links[3].push_back(std::make_pair(i, Ty));
          this->proj2x_link_length[3].push_back(0.0);
          this->proj2y_link_length[3].push_back(2.0);

          this->oriented_links[3].push_back(std::make_pair(i, Tx2my));
          this->proj2x_link_length[3].push_back(sqrt(3.0));
          this->proj2y_link_length[3].push_back(-1.0);
        }
      }
    }
  }

  bool IsSublatticeA(const size_t site) const {
    size_t x = site % (2 * Lx_);
    size_t y = site / (2 * Lx_);
    return (x + y) % 2 == 0;
  }
  size_t Distance(const size_t, const size_t) const override;
 private:
  const size_t Lx_;
  const size_t Ly_;
  const size_t N_;

};

template<size_t DimDof, size_t NumOfCouplingType>
size_t HoneyCombTorusLatticeLink<DimDof, NumOfCouplingType>::Distance(const size_t i, const size_t j) const {
  size_t x_i = i % (2 * this->lx);
  size_t y_i = i / (2 * this->lx);
  size_t x_j = j % (2 * this->lx);
  size_t y_j = j / (2 * this->lx);
  size_t delta_x = (x_i > x_j) ? (x_i - x_j) : (x_j - x_i);
  size_t delta_y = (y_i > y_j) ? (y_i - y_j) : (y_j - y_i);
  return delta_x + delta_y;
}

/**
 * For the case the couplings is different in different direction
 * Please refer to Figure 1 in [Nano Lett. 2021, 21, 10114−10121] for the definition.
 * coupling_structures[0] is single ion anisotropy,
 * coupling_structures[1,2,3] is the J1
 * coupling_structures[4,5,6] is the J2, 4:(1,2), 5:(3,6), 6:(4,5)
 * coupling_structures[7,8,9] is the J3, 7:1, 8:2, 9:3
 * @tparam DimDof
 * @tparam NumOfCouplingType, 9 actually.
 */
template<size_t DimDof, size_t NumOfCouplingType>
class HoneyComb3TorusLatticeLink : public LatticeLink<DimDof, NumOfCouplingType> {
  using CST = CouplingStructure<DimDof>;
 public:
  HoneyComb3TorusLatticeLink(const size_t Lx, const size_t Ly,
                             const std::array<CST, NumOfCouplingType> &coupling_structures) :
      LatticeLink<DimDof, NumOfCouplingType>(coupling_structures, Lx * Ly * 3 * sqrt(3) / 2) {
    const size_t N = Lx * Ly * 2;
    this->lx = Lx;
    this->ly = Ly;
    this->n = N;
    for (size_t i = 0; i < NumOfCouplingType; i++) {
      this->links[i].resize(N);
    }
    /// for single ion anisotropy D, on-site
    for (size_t i = 0; i < N; i++) {
      this->links[0][i] = {i};
    }
    /// for nearest neighbor (NN) J1
    for (size_t i = 0; i < N; i++) {
      size_t x = i % (2 * Lx);
      size_t y = i / (2 * Lx);
      size_t Tx = (x + 1) % (2 * Lx) + y * (2 * Lx);
      size_t Tmx = (x + 2 * Lx - 1) % (2 * Lx) + y * (2 * Lx);
      size_t Ty = x + ((y + 1) % Ly) * (2 * Lx);
      size_t Tmy = x + ((y + Ly - 1) % Ly) * (2 * Lx);
      if ((x + y) % 2 == 0) {
        this->links[1][i] = {Ty};
        this->links[2][i] = {Tx};
        this->links[3][i] = {Tmx};

        this->oriented_links[1].push_back(std::make_pair(i, Ty));
        this->proj2x_link_length[1].push_back(0.0);
        this->proj2y_link_length[1].push_back(1.0);

        this->oriented_links[2].push_back(std::make_pair(i, Tx));
        this->proj2x_link_length[2].push_back(sqrt(3.0) / 2.0);
        this->proj2y_link_length[2].push_back(-1.0 / 2.0);

      } else {
        this->links[1][i] = {Tmy};
        this->links[2][i] = {Tmx};
        this->links[3][i] = {Tx};
        this->oriented_links[3].push_back(std::make_pair(i, Tx));
        this->proj2x_link_length[3].push_back(sqrt(3.0) / 2.0);
        this->proj2y_link_length[3].push_back(1.0 / 2.0);
      }
    }

    /// for NNN J2
    bool is_j2_zero;
    if (NumOfCouplingType <= 4) {
      is_j2_zero = true;
    } else {
      is_j2_zero = coupling_structures[4].IsZero()
          && coupling_structures[5].IsZero()
          && coupling_structures[6].IsZero();
    }

    if (NumOfCouplingType > 4 && !is_j2_zero) {
      for (size_t i = 0; i < N; i++) {
        size_t x = i % (2 * Lx);
        size_t y = i / (2 * Lx);
        size_t Tx2 = (x + 2) % (2 * Lx) + y * (2 * Lx);
        size_t Tmx2 = (x + 2 * Lx - 2) % (2 * Lx) + y * (2 * Lx);
        size_t Txy = (x + 1) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
        size_t Txmy = (x + 1) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);
        size_t Tmxy = (x + 2 * Lx - 1) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
        size_t Tmxmy = (x + 2 * Lx - 1) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);

        this->links[4][i] = {Tx2, Tmx2};
        this->links[5][i] = {Txy, Tmxmy};
        this->links[6][i] = {Txmy, Tmxy};
        this->oriented_links[4].push_back(std::make_pair(i, Tx2));
        this->proj2x_link_length[4].push_back(sqrt(3.0));
        this->proj2y_link_length[4].push_back(0.0);

        this->oriented_links[5].push_back(std::make_pair(i, Txy));
        this->proj2x_link_length[5].push_back(sqrt(3.0) / 2.0);
        this->proj2y_link_length[5].push_back(3.0 / 2.0);

        this->oriented_links[6].push_back(std::make_pair(i, Txmy));
        this->proj2x_link_length[6].push_back(sqrt(3.0) / 2.0);
        this->proj2y_link_length[6].push_back(-3.0 / 2.0);
      }
    }

    /// for TNN J3
    bool is_j3_zero;
    if (NumOfCouplingType <= 7) {
      is_j3_zero = true;
    } else {
      is_j3_zero = coupling_structures[7].IsZero()
          && coupling_structures[8].IsZero()
          && coupling_structures[9].IsZero();
    }
    if (NumOfCouplingType > 7 && !is_j3_zero) {
      for (size_t i = 0; i < N; i++) {
        size_t x = i % (2 * Lx);
        size_t y = i / (2 * Lx);
        if ((x + y) % 2 == 0) {
          size_t Tmy = x + ((y + Ly - 1) % Ly) * (2 * Lx);
          size_t Tx2y = (x + 2) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
          size_t Tmx2y = (x + 2 * Lx - 2) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
          this->links[7][i] = {Tmy};
          this->links[8][i] = {Tx2y};
          this->links[9][i] = {Tmx2y};
          this->oriented_links[7].push_back(std::make_pair(i, Tmy));
          this->proj2x_link_length[7].push_back(0.0);
          this->proj2y_link_length[7].push_back(-2.0);

          this->oriented_links[8].push_back(std::make_pair(i, Tx2y));
          this->proj2x_link_length[8].push_back(sqrt(3.0));
          this->proj2y_link_length[8].push_back(1.0);

        } else {
          size_t Ty = x + ((y + 1) % Ly) * (2 * Lx);
          size_t Tx2my = (x + 2) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);
          size_t Tmx2my = (x + 2 * Lx - 2) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);
          this->links[7][i] = {Ty};
          this->links[8][i] = {Tmx2my};
          this->links[9][i] = {Tx2my};
          this->oriented_links[9].push_back(std::make_pair(i, Tx2my));
          this->proj2x_link_length[9].push_back(sqrt(3.0));
          this->proj2y_link_length[9].push_back(-1.0);
        }
      }
    }
  }
  size_t Distance(size_t i, size_t j) const override;
};

template<size_t DimDof, size_t NumOfCouplingType>
size_t HoneyComb3TorusLatticeLink<DimDof, NumOfCouplingType>::Distance(size_t i, size_t j) const {
  size_t x_i = i % (2 * this->lx);
  size_t y_i = i / (2 * this->lx);
  size_t x_j = j % (2 * this->lx);
  size_t y_j = j / (2 * this->lx);
  size_t delta_x = (x_i > x_j) ? (x_i - x_j) : (x_j - x_i);
  size_t delta_y = (y_i > y_j) ? (y_i - y_j) : (y_j - y_i);
  return delta_x + delta_y;
}
#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICELINK_H
