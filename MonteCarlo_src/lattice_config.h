//
// Created by Hao-Xin on 2022/5/22.
//



#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICE_CONFIG_H
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICE_CONFIG_H
#include "local_dof.h"
#include "lattice_link.h"
#include "swcluster.h"

template <size_t DimDof>
using SpinConfig = std::vector<LocalDOF<DimDof>>;

template <size_t DimDof>
class LatticeConfig {
  using LocalDOFT = LocalDOF<DimDof>;
  using SpinConfigT = SpinConfig<DimDof>;
 public:
  LatticeConfig(const size_t N) : N_(N), config_(N) {}
  LatticeConfig(const LatticeConfig& lattice_config) :
      N_(lattice_config.N_), config_(lattice_config.config_) {}

  LatticeConfig &operator=(const LatticeConfig& rhs) {
    N_ = rhs.N_;
    config_ = rhs.config_;
    return *this;
  }


  void Random() {
    for(auto& local_dof : config_) {
      local_dof.Random();
    }
  }

  void Flip(const LocalDOFT& axis) {
    for(auto& local_dof : config_) {
      local_dof.Flip(axis);
    }
  }

  void Flip(const LocalDOFT& axis, const SWCluster& cluster) {
    for(auto site : cluster) {
      config_[site].Flip(axis);
    }
  }

  void SetSiteSpin(const size_t site, const LocalDOFT &spin) {
    config_[site] = spin;
  }
  /* calculate the energy according to the lattice link (including coupling structures)
   * The Hamiltonian is defined with extra minus sign
   */
  template <size_t InteractionRange>
  double Energy(const LatticeLink<DimDof, InteractionRange> & lattice_link) const {
    auto& coupling_structures = lattice_link.coupling_structures;
    auto& links = lattice_link.links;

    ///single ion anistropy
    double e0(0.0);
    auto& D = coupling_structures[0];
    for(auto& spin : config_) {
      e0 += spin * (D * spin);
    }

    ///two-site interaction, calculate twice so half finally
    double e2(0.0);
    for(size_t inter_range = 1; inter_range < InteractionRange; inter_range++) {
      auto& J = coupling_structures[inter_range];
      auto& link_set = links[inter_range];
      assert(N_ == link_set.size());
      for(size_t site1 = 0; site1 < N_; site1++) {
        auto& spin1 = config_[site1];
        auto& site2_set = link_set[site1];
        for(auto site2 : site2_set) {
          auto& spin2 = config_[site2];
          e2 += spin1 * (J * spin2);
        }
      }
    }
    return -( e0 + e2 / 2 );
  }

  double SumOverComponent(const size_t i) const {
    double sum = 0.0;
    for(auto& spin : config_) {
      sum += spin.GetCoor(i);
    }
    return sum;
  }

  void PrintSign(const LocalDOFT& n, const size_t lx) const {
    for(size_t i = 0; i < N_; i++) {
      if(i%lx == 0) {
        std::cout <<"\n";
      }
      if(config_[i] * n > 0) {
        std::cout <<"+";
      } else{
        std::cout << "-";
      }
    }
    std::cout << std::endl;
  }

  /**
   * using when growing the cluster
   *
   * @param site1
   * @param site2
   * @param Keff = beta * J_eff, where J_eff = n*J*n
   * @param n    the axis
   * @return
   */
  bool Active(const size_t site1,
              const size_t site2,
              const double Keff,
              const LocalDOFT& n) const {
    double DeltaE2 = 2 * Keff * ( config_[site1] * n ) * ( config_[site2] * n );
    if(DeltaE2 <= 0) {
      return false;
    } else {
      double P = exp(-DeltaE2);
      std::uniform_real_distribution<double> u(0, 1);
      if( u(random_engine) > P) {
        return true;
      } else {
        return false;
      }
    }
  }

  //TODO: optimize
  template <size_t InteractionRange>
  double EnergyRelevantToSite(
      const size_t site,
      const LatticeLink<DimDof, InteractionRange> &lattice_link
      ) const {
    auto& coupling_structures = lattice_link.coupling_structures;
    auto& links = lattice_link.links;
    auto& spin = config_[site];
    double e0 = spin * ( coupling_structures[0] * spin);

    double e2 = 0.0;
    for(size_t inter_range = 1; inter_range < InteractionRange; inter_range++) {
      auto& J = coupling_structures[inter_range];
      auto& link_set = links[inter_range];
      assert(N_ == link_set.size());
      auto& site2_set = link_set[site];
      for(auto site2 : site2_set) {
        auto& spin2 = config_[site2];
        e2 += spin * (J * spin2);
      }
    }
    return -( e0 + e2 );
  }

  /// the case we assume the spin in the site is replaced by another spin
  template <size_t InteractionRange>
  double EnergyRelevantToSite(
      const size_t site,
      const LatticeLink<DimDof, InteractionRange> &lattice_link,
      const LocalDOFT &spin
  ) const {
    auto& coupling_structures = lattice_link.coupling_structures;
    auto& links = lattice_link.links;
    double e0 = spin * ( coupling_structures[0] * spin);

    double e2 = 0.0;
    for(size_t inter_range = 1; inter_range < InteractionRange; inter_range++) {
      auto& J = coupling_structures[inter_range];
      auto& link_set = links[inter_range];
      assert(N_ == link_set.size());
      auto& site2_set = link_set[site];
      for(auto site2 : site2_set) {
        auto& spin2 = config_[site2];
        e2 += spin * (J * spin2);
      }
    }
    return -( e0 + e2 );
  }

  template <size_t InteractionRange>
  double EnergyDifferenceFlipSite(
      const size_t site,
      const LocalDOFT &fliped_spin,
      const LatticeLink<DimDof, InteractionRange>& lattice_link
      ) const {
    return - EnergyRelevantToSite(site, lattice_link) + EnergyRelevantToSite(site, lattice_link, fliped_spin);
  }

  template <size_t InteractionRange>
  double EnergyDifferenceFlipCluster(
      const LocalDOFT& n,
      const SWCluster& cluster,
      const LatticeLink<DimDof, InteractionRange>& lattice_link
      ) const {
    auto [para_config, perp_config] = ComponentDecompose_(n);
    double h2(0.0), sia(0.0); //H2, SIA
    auto& D = lattice_link.coupling_structures[0];
    for(size_t i : cluster) {
      sia += (para_config[i] * (D * perp_config[i]));
    }
    sia *= 2;
    for(size_t inter_type = 1; inter_type < InteractionRange; inter_type++){
      auto& J = lattice_link.coupling_structures[inter_type];
      for(size_t i : cluster) {
        auto j_set = lattice_link.GetSitesLinkedTo(i, inter_type);

        LocalDOFT sum_j_perp; sum_j_perp.Zero();
        for(auto j : j_set) {
          sum_j_perp += perp_config[j];
        }
        h2 += sum_j_perp * (J * para_config[i]);
       }
    }

    return   2 * (h2 + sia);
  }

  /// interaction support upto J3, and Heisenberg or XY like (no anisotropy except SIA terms)
  template <size_t NumOfCouplingType>
  double StiffnessHoneycomb(
      const LatticeLink<DimDof, NumOfCouplingType>& lattice_link,
      const double beta
      ) const {
    double inner_product_sum[NumOfCouplingType - 1];
    double cross_product_sum[NumOfCouplingType - 1];
    for(size_t inter_range = 1; inter_range < NumOfCouplingType; inter_range++) {
      auto& oriented_link = lattice_link.oriented_links[inter_range];
      double sum_a(0.0), sum_b(0.0);
      for(auto [site1, site2] : oriented_link) {
        sum_a += Project2DInnerProduct(config_[site1], config_[site2]);

#ifndef NDEBUG
        double z_comp_of_cross_product = ZCompOfCrossProduct(config_[site1], config_[site2]);
        if(DimDof == 2) {
          double sin_value = ::SinDiff(config_[site1], config_[site2]);
          assert(sin_value = -z_comp_of_cross_product);
        }
        sum_b += z_comp_of_cross_product;
#else
        sum_b += ZCompOfCrossProduct(config_[site1], config_[site2]);
#endif
      }
      inner_product_sum[inter_range-1] = sum_a ;
      cross_product_sum[inter_range-1] = sum_b ;
    }
    double coefficient_by_lattice[3] = {0.5, 1.5, 2.0};
    double stiffness = 0.0;
    for(size_t i = 0; i < NumOfCouplingType - 1; i++) {
      const size_t inter_range = i+1;
      double j = lattice_link.coupling_structures[inter_range].coupling_const[0][0];
      stiffness += coefficient_by_lattice[i] * ( j * inner_product_sum[i] - beta * j * j * cross_product_sum[i] * cross_product_sum[i] ) / N_;
    }
    return stiffness;
  }

  template <size_t NumOfCouplingType>
  double StiffnessSquare(
      const LatticeLink<DimDof, NumOfCouplingType>& lattice_link,
      const double beta
  ) const {
    double inner_product_sum[NumOfCouplingType - 1];
    double cross_product_sum[NumOfCouplingType - 1];
    for(size_t inter_range = 1; inter_range < NumOfCouplingType; inter_range++) {
      auto& oriented_link = lattice_link.oriented_links[inter_range];
      double sum_a(0.0), sum_b(0.0);
      for(auto [site1, site2] : oriented_link) {
        sum_a += Project2DInnerProduct(config_[site1], config_[site2]);

#ifndef NDEBUG
        double z_comp_of_cross_product = ZCompOfCrossProduct(config_[site1], config_[site2]);
        if(DimDof == 2) {
          double sin_value = ::SinDiff(config_[site1], config_[site2]);
          assert(sin_value = -z_comp_of_cross_product);
        }
        sum_b += z_comp_of_cross_product;
#else
        sum_b += ZCompOfCrossProduct(config_[site1], config_[site2]);
#endif
      }
      inner_product_sum[inter_range-1] = sum_a ;
      cross_product_sum[inter_range-1] = sum_b ;
    }
    double coefficient_by_lattice[3] = {0.5, 1, 2.0}; // the only line different with Honeycomb
    double stiffness = 0.0;
    for(size_t i = 0; i < NumOfCouplingType - 1; i++) {
      const size_t inter_range = i+1;
      double j = lattice_link.coupling_structures[inter_range].coupling_const[0][0];
      stiffness += coefficient_by_lattice[i] * ( j * inner_product_sum[i] - beta * j * j * cross_product_sum[i] * cross_product_sum[i] ) / N_;
    }
    return stiffness;
  }

 private:

  std::pair<SpinConfigT, SpinConfigT> ComponentDecompose_(
      const LocalDOFT& n
      ) const {
    std::pair<SpinConfigT, SpinConfigT> res;
    res.first.reserve(N_);
    res.second.reserve(N_);
    for(size_t i = 0; i < N_; i++) {
//      [ res.first[i],  res.send[i] ] = config_[i]
      auto [para, perp] =  config_[i].DecompSpin(n);
      res.first[i] = para;
      res.second[i] = perp;
    }
    return res;
  }

  size_t N_; //lattice size
  SpinConfigT config_;
};




#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICE_CONFIG_H
