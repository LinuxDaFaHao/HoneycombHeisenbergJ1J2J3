//
// Created by Hao-Xin on 2022/5/24.
//

#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICELINK_H
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICELINK_H

#include <vector>
#include "./coupling.h"

/// define the Hamiltonian
template <size_t DimDof, size_t NumOfCouplingType>
class LatticeLink {
  using CST = CouplingStructure<DimDof>;
  using OrientedBond = std::pair<size_t, size_t>;
 public:
  LatticeLink(void) = default;

  LatticeLink(const std::array<CST, NumOfCouplingType> &coupling_structures) :
      coupling_structures(coupling_structures) {}

  LatticeLink(const LatticeLink &lattice_link) :
      coupling_structures(lattice_link.coupling_structures),
      links(lattice_link.links) {}

  std::vector<std::pair<size_t, size_t>> GetBondLinkedTo(const size_t site, const size_t InterType){
    std::vector<std::pair<size_t, size_t>> bonds;
    const auto& site2_set = links[InterType][site];
    for(auto& site2 : site2_set) {
      bonds.push_back(std::make_pair(site, site2));
    }
    return bonds;
  }

  const std::vector<size_t>&  GetSitesLinkedTo(const size_t site, const size_t InterType) const {
    return links[InterType][site];
  }

  size_t GetNumofInteraction(const size_t inter_type) const {
    size_t num = 0;
    for(auto& site2_set : links[inter_type]){
      num += site2_set.size();
    }
    if(inter_type > 0) {
      num = num / 2;
    }
    return num;
  }

  ///0, single ion anisotropy D;
  ///1, J1; 2, J2; 3, J3;
  std::array<CST, NumOfCouplingType> coupling_structures;
  std::array<std::vector<std::vector<size_t>>, NumOfCouplingType> links;
  /// links[i][j] is a std::vector<size_t>, which is the set of under interaction Ji(D) the sites link to site j.
  std::array<std::vector<OrientedBond>, NumOfCouplingType> oriented_links;

// private:
//  const size_t Lx_;
//  const size_t Ly_;
//  const size_t N_;
};

/// all link
template <size_t DimDof, size_t NumOfCouplingType>
class SVWLatticeLink : public LatticeLink<DimDof, NumOfCouplingType> {
  using CST = CouplingStructure<DimDof>;
 public:
  SVWLatticeLink(const size_t N, const std::array<CST, NumOfCouplingType> &coupling_structures) :
      LatticeLink<DimDof, NumOfCouplingType>(coupling_structures) {
    for(size_t i = 0; i < NumOfCouplingType; i++) {
      this->links[i].resize(N);
    }
    /// SIA
    for(size_t i = 0; i < N; i++) {
      this->links[0][i] = {i};
    }
    /// Interaction
    for(size_t i = 0; i < N; i++) {
      this->links[1][i] = {};
      this->links[1][i].reserve(N - 1);
      for(size_t j = 0; j < N; j++){
        if(i!=j)
          this->links[1][i].push_back(j);
      }
    }
  }
};

template <size_t DimDof, size_t NumOfCouplingType>
class SquareTorusLatticeLink : public LatticeLink<DimDof, NumOfCouplingType> {
  using CST = CouplingStructure<DimDof>;
 public:
  SquareTorusLatticeLink(const size_t Lx, const size_t Ly,
                         const std::array<CST, NumOfCouplingType> &coupling_structures) :
      LatticeLink<DimDof, NumOfCouplingType>(coupling_structures) {
    const size_t N = Lx * Ly;
    for(size_t i = 0; i < NumOfCouplingType; i++) {
      this->links[i].resize(N);
      if(i > 0) this->oriented_links[i].reserve(2 * N);
    }
    /// SIA
    for(size_t i = 0; i < N; i++) {
      this->links[0][i] = {i};
    }
    /// NN
    for(size_t i = 0; i < N; i++) {
      size_t x = i%(Lx);
      size_t y = i/(Lx);
      size_t Tx = (x + 1) % (Lx) + y * (Lx);
      size_t Tmx = (x + Lx - 1) % (Lx) + y *(Lx);
      size_t Ty = x + ( (y + 1) % Ly ) * (Lx);
      size_t Tmy = x + ( (y + Ly - 1) % Ly ) * (Lx);
      this->links[1][i] = {Tx, Ty, Tmx, Tmy};

      this->oriented_links[1].push_back(std::make_pair(i, Tx));
      this->oriented_links[1].push_back(std::make_pair(i, Ty));
    }
    /// NNN
    if( NumOfCouplingType > 2) {
      for(size_t i = 0; i < N; i++) {
        size_t x = i%(Lx);
        size_t y = i/(Lx);
        size_t Txy = (x + 1) % (Lx) + ( (y + 1) % Ly ) * (Lx);
        size_t Txmy = (x + 1) % (Lx) + ( (y + Ly - 1) % Ly ) * (Lx);
        size_t Tmxy = (x + Lx - 1) % (Lx) + ( (y + 1) % Ly ) * (Lx);
        size_t Tmxmy = (x + Lx - 1) % (Lx) + ( (y + Ly - 1) % Ly ) * (Lx);

        this->links[2][i] = {Txy, Txmy, Tmxy, Tmxmy}; //sort?

        this->oriented_links[2].push_back(std::make_pair(i, Txy));
        this->oriented_links[2].push_back(std::make_pair(i, Txmy));
      }
    }


    /// TNN
    if( NumOfCouplingType > 3) {
      for (size_t i = 0; i < N; i++) {
        size_t x = i % (Lx);
        size_t y = i / (Lx);
        size_t Tx2 = (x + 2) % (Lx) + y * (Lx);
        size_t Tmx2 = (x + Lx - 2) % (Lx) + y * (Lx);
        size_t Ty2 = x + ((y + 2) % Ly) * (Lx);
        size_t Tmy2 = x + ((y + Ly - 2) % Ly) * (Lx);
        this->links[3][i] = {Tx2, Ty2, Tmx2,  Tmy2};

        this->oriented_links[3].push_back(std::make_pair(i, Tx2));
        this->oriented_links[3].push_back(std::make_pair(i, Ty2));
      }
    }
  }
};

template <size_t DimDof, size_t NumOfCouplingType>
class HoneyCombTorusLatticeLink : public LatticeLink<DimDof, NumOfCouplingType> {
  using CST = CouplingStructure<DimDof>;
 public:
  HoneyCombTorusLatticeLink(const size_t Lx, const size_t Ly,
                            const std::array<CST, NumOfCouplingType> &coupling_structures) :
      LatticeLink<DimDof, NumOfCouplingType>(coupling_structures),
      Lx_(Lx), Ly_(Ly), N_(2 * Lx * Ly) {
    const size_t N = Lx * Ly * 2;
    for(size_t i = 0; i < NumOfCouplingType; i++) {
      this->links[i].resize(N);
      if(i > 0) this->oriented_links[i].reserve(3 * N); //reserve more for J1, J3
    }
    /// for single ion anisotropy D, on-site
    for(size_t i = 0; i < N; i++) {
      this->links[0][i] = {i};
    }
    /// for nearest neighbor (NN) J1
    for(size_t i = 0; i < N; i++) {
      size_t x = i%(2*Lx);
      size_t y = i/(2*Lx);
      size_t Tx = (x + 1) % (2 * Lx) + y * (2*Lx);
      size_t Tmx = (x + 2 * Lx - 1) % (2 * Lx) + y *(2*Lx);
      size_t Ty = x + ( (y + 1) % Ly ) * (2*Lx);
      size_t Tmy = x + ( (y + Ly - 1) % Ly ) * (2*Lx);
      if( (x+y)%2 == 0) {
        this->links[1][i] = {Tmx, Tx, Ty}; //sort?
        this->oriented_links[1].push_back(std::make_pair(i, Tx));
        this->oriented_links[1].push_back(std::make_pair(i, Ty));
      } else {
        this->links[1][i] = {Tmx, Tx, Tmy};
        this->oriented_links[1].push_back(std::make_pair(i, Tx));
      }
    }

    /// for NNN J2
    if( NumOfCouplingType > 2) {
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
        this->oriented_links[2].push_back(std::make_pair(i, Txy));
        this->oriented_links[2].push_back(std::make_pair(i, Txmy));
      }
    }

    /// for TNN J3
    if( NumOfCouplingType > 3) {
      for (size_t i = 0; i < N; i++) {
        size_t x = i % (2 * Lx);
        size_t y = i / (2 * Lx);
        if ((x + y) % 2 == 0) {
          size_t Tmy = x + ((y + Ly - 1) % Ly) * (2 * Lx);
          size_t Tx2y = (x + 2) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
          size_t Tmx2y = (x + 2 * Lx - 2) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
          this->links[3][i] = {Tmy, Tx2y, Tmx2y};

          this->oriented_links[3].push_back(std::make_pair(i, Tx2y));
          this->oriented_links[3].push_back(std::make_pair(i, Tmx2y));
        } else {
          size_t Ty = x + ((y + 1) % Ly) * (2 * Lx);
          size_t Tx2my = (x + 2) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);
          size_t Tmx2my = (x + 2 * Lx - 2) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);
          this->links[3][i] = {Ty, Tx2my, Tmx2my};

          this->oriented_links[3].push_back(std::make_pair(i, Ty));
        }


      }
    }
  }


  bool IsSublatticeA(const size_t site) const {
    size_t x = site % (2 * Lx_);
    size_t y = site / (2 * Lx_);
    return (x+y)%2 == 0;
  }
 private:
  const size_t Lx_;
  const size_t Ly_;
  const size_t N_;
};


/** TODO class
 * For the case the couplings is different in different direction
 * coupling_structures[0] is single ion anisotropy,
 * coupling_structures[1,2,3] is the J1
 * coupling_structures[4,5,6] is the J2
 * coupling_structures[7,8,9] is the J3
 * @tparam DimDof
 * @tparam NumOfCouplingType, 9 actually.
 */
template <size_t DimDof, size_t NumOfCouplingType>
class HoneyComb3TorusLatticeLink : public LatticeLink<DimDof, NumOfCouplingType> {
  using CST = CouplingStructure<DimDof>;
 public:
  HoneyComb3TorusLatticeLink(const size_t Lx, const size_t Ly,
                            const std::array<CST, NumOfCouplingType> &coupling_structures) :
      LatticeLink<DimDof, NumOfCouplingType>(coupling_structures) {
    const size_t N = Lx * Ly * 2;
    for(size_t i = 0; i < NumOfCouplingType; i++) {
      this->links[i].resize(N);
    }
    /// for single ion anisotropy D, on-site
    for(size_t i = 0; i < N; i++) {
      this->links[0][i] = {i};
    }
    /// for nearest neighbor (NN) J1
    for(size_t i = 0; i < N; i++) {
      size_t x = i%(2*Lx);
      size_t y = i/(2*Lx);
      size_t Tx = (x + 1) % (2 * Lx) + y * (2*Lx);
      size_t Tmx = (x + 2 * Lx - 1) % (2 * Lx) + y *(2*Lx);
      size_t Ty = x + ( (y + 1) % Ly ) * (2*Lx);
      size_t Tmy = x + ( (y + Ly - 1) % Ly ) * (2*Lx);
      if( (x+y)%2 == 0) {
        this->links[1][i] = {Ty};
        this->links[2][i] = {Tx};
        this->links[3][i] = {Tmx};
      } else {
        this->links[1][i] = {Tmy};
        this->links[2][i] = {Tmx};
        this->links[3][i] = {Tx};
      }
    }

    /// for NNN J2
    if( NumOfCouplingType > 2) {
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
      }
    }

    /// for TNN J3
    if( NumOfCouplingType > 3) {
      for (size_t i = 0; i < N; i++) {
        size_t x = i % (2 * Lx);
        size_t y = i / (2 * Lx);
        if ((x + y) % 2 == 0) {
          size_t Tmy = x + ((y + Ly - 1) % Ly) * (2 * Lx);
          size_t Tx2y = (x + 2) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
          size_t Tmx2y = (x + 2 * Lx - 2) % (2 * Lx) + ((y + 1) % Ly) * (2 * Lx);
          this->links[3][i] = {Tmy, Tx2y, Tmx2y};
        } else {
          size_t Ty = x + ((y + 1) % Ly) * (2 * Lx);
          size_t Tx2my = (x + 2) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);
          size_t Tmx2my = (x + 2 * Lx - 2) % (2 * Lx) + ((y + Ly - 1) % Ly) * (2 * Lx);
          this->links[3][i] = {Ty, Tx2my, Tmx2my};
        }
      }
    }
  }
};


#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICELINK_H
