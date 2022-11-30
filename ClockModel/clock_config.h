//
// Created by Hao-Xin on 2022/10/2.
//

#ifndef HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCK_CONFIG_H_
#define HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCK_CONFIG_H_

#include "clock_dof.h"
#include "../MonteCarlo_src/lattice_link.h"
#include "../MonteCarlo_src/swcluster.h"

inline std::uniform_real_distribution<double> u_double(0, 1);

template<size_t q>
using ClockVector = std::vector<ClockDOF<q>>;
using RandomField = std::vector<double>;

inline const size_t clock_dim_dof = 1;
inline const size_t clock_interaction_range = 2;

template<size_t q>
class ClockConfig {
  using ClockVectorT = ClockVector<q>;
 public:
  ClockConfig(const size_t N) : N_(N), config_(N) {}
  ClockConfig(const ClockConfig &clock_config) :
      N_(clock_config.N_),
      config_(clock_config.N_) {}
  ClockConfig &operator=(const ClockConfig &rhs) {
    N_ = rhs.N_;
    config_ = rhs.config_;
    return *this;
  }

  void Random() {
    for (auto &local_dof: config_) {
      local_dof.Random();
    }
  }

  void SetSiteSpin(const size_t site, const ClockDOF<q> &spin) {
    config_[site] = spin;
  }

  void Flip2(const ClockDOF<q> &axis) {
    for (auto &local_dof: config_) {
      local_dof.Flip2(axis);
    }
  }

  void Flip2(const ClockDOF<q> &axis, const SWCluster &cluster) {
    for (auto site: cluster) {
      config_[site].Flip2(axis);
    }
  }

  double Energy(const LatticeLink<clock_dim_dof, clock_interaction_range> &lattice_link,
                const RandomField &hp_x = {},
                const RandomField &hp_y = {}) const {
    auto &coupling_structures = lattice_link.coupling_structures;
    auto &links = lattice_link.links;
    double e2 = 0.0;
    for (size_t inter_range = 1; inter_range < clock_interaction_range; inter_range++) {
      double j = coupling_structures[inter_range].coupling_const[0][0];
      auto &link_set = links[inter_range];
      assert(N_ == link_set.size());
      for (size_t site1 = 0; site1 < N_; site1++) {
        auto &spin1 = config_[site1];
        auto &site2_set = link_set[site1];
        for (auto site2: site2_set) {
          auto &spin2 = config_[site2];
          e2 += j * (spin1 * spin2);
        }
      }
    }

    double e_h = 0.0;
    if (hp_x.size() > 0) {
      for (size_t site = 0; site < N_; site++) {
        e_h += config_[site].Cos2Theta() * hp_x[site] + config_[site].Sin2Theta() * hp_y[site];
      }
    }
    return -e2 / 2 - e_h;
  }

  double EnergyRelevantToSite(
      const size_t site,
      const LatticeLink<clock_dim_dof, clock_interaction_range> &lattice_link,
      const RandomField &hp_x = {},
      const RandomField &hp_y = {}
  ) const {
    auto &coupling_structures = lattice_link.coupling_structures;
    auto &links = lattice_link.links;
    auto &spin = config_[site];

    double e2 = 0.0;
    for (size_t inter_range = 1; inter_range < clock_interaction_range; inter_range++) {
      double j = coupling_structures[inter_range].coupling_const[0][0];
      auto &link_set = links[inter_range];
      assert(N_ == link_set.size());
      auto &site2_set = link_set[site];
      for (auto site2: site2_set) {
        auto &spin2 = config_[site2];
        e2 += j * (spin * spin2);
      }
    }

    double e_h = 0.0;
    if (hp_x.size() > 0) {
      e_h += spin.Cos2Theta() * hp_x[site] + spin.Sin2Theta() * hp_y[site];
    }

    return -e2 - e_h;
  }

  /// the case we assume the spin in the site is replaced by another spin
  double ReplacedEnergyRelevantToSite(
      const size_t site,
      const LatticeLink<clock_dim_dof, clock_interaction_range> &lattice_link,
      const ClockDOF<q> &spin,
      const RandomField &hp_x = {},
      const RandomField &hp_y = {}
  ) const {
    auto &coupling_structures = lattice_link.coupling_structures;
    auto &links = lattice_link.links;

    double e2 = 0.0;
    for (size_t inter_range = 1; inter_range < clock_interaction_range; inter_range++) {
      double j = coupling_structures[inter_range].coupling_const[0][0];
      auto &link_set = links[inter_range];
      assert(N_ == link_set.size());
      auto &site2_set = link_set[site];
      for (auto site2: site2_set) {
        auto &spin2 = config_[site2];
        e2 += j * (spin * spin2);
      }
    }

    double e_h = 0.0;
    if (hp_x.size() > 0) {
      e_h += spin.Cos2Theta() * hp_x[site] + spin.Sin2Theta() * hp_y[site];
    }

    return -e2 - e_h;
  }

  double EnergyDifferenceFlipSite(
      const size_t site,
      const ClockDOF<q> &flipped_spin,
      const LatticeLink<clock_dim_dof, clock_interaction_range> &lattice_link,
      const RandomField &hp_x = {},
      const RandomField &hp_y = {}
  ) const {
    return -EnergyRelevantToSite(site, lattice_link, hp_x, hp_y)
        + ReplacedEnergyRelevantToSite(site, lattice_link, flipped_spin, hp_x, hp_y);
  }

  double NotActiveProbability(const size_t site1,
                         const size_t site2,
                         const double k_eff,
                         const ClockDOF<q> &axis) const {
    double
        DeltaE2 = 2 * k_eff * ((config_[site1] * config_[site2])
        - (config_[site1] * axis) * (config_[site2] * axis));
//    if(DeltaE2 <=0) return 1.0;
//    else return exp(-DeltaE2);
    return exp(-DeltaE2);
  }

  bool Active(const size_t site1,
              const size_t site2,
              const double k_eff,
              const ClockDOF<q> &axis) const {
    double
        DeltaE2 = 2 * k_eff * ((config_[site1] * config_[site2])
            - (config_[site1] * axis) * (config_[site2] * axis));
    if (DeltaE2 <= 0) {
      return false;
    } else {
      double P = exp(-DeltaE2);
      if (u_double(random_engine) > P) {
        return true;
      } else {
        return false;
      }
    }
  }

  double EnergyDifferenceOfRandomFieldFlipCluster(
      const ClockDOF<q> &axis,
      const SWCluster &cluster,
      const RandomField &h_x,
      const RandomField &h_y
  ) {
    double delta_e(0.0);
    for (size_t site: cluster) {
      auto clock = config_[site];
      double e_original = -(clock.Cos2Theta() * h_x[site] + clock.Sin2Theta() * h_y[site]);
      clock.Flip2(axis);
      double e_flipped = -(clock.Cos2Theta() * h_x[site] + clock.Sin2Theta() * h_y[site]);
      delta_e += (e_flipped - e_original);
    }
    return delta_e;
  }

  double FerromagneticMSquare() {
    LocalDOF<2> sum_spin({0.0, 0.0});
    for (size_t site = 0; site < N_; site++) {
      sum_spin = config_[site] + sum_spin;
    }
    return sum_spin.Norm();
  }

 private:
  size_t N_; //site number
  ClockVectorT config_;

};

#endif //HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCK_CONFIG_H_
