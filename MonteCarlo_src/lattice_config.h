//
// Created by Hao-Xin on 2022/5/22.
//



#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICE_CONFIG_H
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICE_CONFIG_H
#include "local_dof.h"
#include "lattice_link.h"
#include "swcluster.h"
#include <complex>
#include <cmath>

template<size_t DimDof>
using SpinConfig = std::vector<LocalDOF<DimDof>>;
inline std::uniform_real_distribution<double> u_double(0, 1);

template<size_t DimDof>
class LatticeConfig {
  using LocalDOFT = LocalDOF<DimDof>;
  using SpinConfigT = SpinConfig<DimDof>;
 public:
  LatticeConfig(const size_t N) : N_(N), config_(N), parallel_components_(N), perpendicular_components_(N) {}
  LatticeConfig(const LatticeConfig &lattice_config) :
      N_(lattice_config.N_),
      config_(lattice_config.config_),
      parallel_components_(lattice_config.N_),
      perpendicular_components_(lattice_config.N_) {}

  LatticeConfig &operator=(const LatticeConfig &rhs) {
    N_ = rhs.N_;
    config_ = rhs.config_;
    parallel_components_ = SpinConfigT(rhs.N_);
    perpendicular_components_ = SpinConfigT(rhs.N_);
    return *this;
  }

  void Random(std::default_random_engine& generator) {
    for (auto &local_dof: config_) {
      local_dof.Random(generator);
    }
  }

  void Random() {
    Random(random_engine);
  }

  void Flip(const LocalDOFT &axis) {
    for (auto &local_dof: config_) {
      local_dof.Flip(axis);
    }
  }

  void Flip(const LocalDOFT &axis, const SWCluster &cluster) {
    for (auto site: cluster) {
      config_[site].Flip(axis);
    }
  }

  void SetSiteSpin(const size_t site, const LocalDOFT &spin) {
    config_[site] = spin;
  }
  /* calculate the energy according to the lattice link (including coupling structures)
   * The Hamiltonian is defined with extra minus sign
   */
  template<size_t InteractionRange>
  double Energy(const LatticeLink<DimDof, InteractionRange> &lattice_link) const {
    auto &coupling_structures = lattice_link.coupling_structures;
    auto &links = lattice_link.links;

    ///single ion anistropy
    double e0(0.0);
    auto &D = coupling_structures[0];
    for (auto &spin: config_) {
      e0 += spin * (D * spin);
    }

    ///two-site interaction, calculate twice so half finally
    double e2(0.0);
    for (size_t inter_range = 1; inter_range < InteractionRange; inter_range++) {
      auto &J = coupling_structures[inter_range];
      auto &link_set = links[inter_range];
      assert(N_ == link_set.size());
      for (size_t site1 = 0; site1 < N_; site1++) {
        auto &spin1 = config_[site1];
        auto &site2_set = link_set[site1];
        for (auto site2: site2_set) {
          auto &spin2 = config_[site2];
          e2 += spin1 * (J * spin2);
        }
      }
    }
    return -(e0 + e2 / 2);
  }

  template<size_t InteractionRange>
  double Energy(const LatticeLink<DimDof, InteractionRange> &lattice_link,
                const double hp) const {
    double e0 = Energy(lattice_link);
    double e3 = 0.0;
    for (auto &spin: config_) {
      e3 += spin.Cos6Theta();
    }
    return e0 + hp * e3;
  }

  double SumOverComponent(const size_t i) const {
    double sum = 0.0;
    for (auto &spin: config_) {
      sum += spin.GetCoor(i);
    }
    return sum;
  }

  /// For honeycomb lattice, anisotropy, D3d symmetry, complex order parameter
  std::complex<double> LocalComplexOrderParameter(const size_t lx,
                                                  const size_t ly,
                                                  const size_t center_site,
                                                  LocalDOFT &gsoa,
                                                  LocalDOFT &gsob,
                                                  LocalDOFT &gsoc) {
    size_t y = center_site / (2 * lx);
    size_t x = center_site % (2 * lx);
    size_t Tmx = (x + 2 * lx - 1) % (2 * lx) + y * (2 * lx);
    size_t Tx = (x + 1) % (2 * lx) + y * (2 * lx);
    size_t Tmy = x + ((y + ly - 1) % ly) * (2 * lx);

    size_t site1 = Tmx, site2 = Tmy, site3 = Tx, site4 = center_site;

    double term1 = (config_[site1] + (-config_[site2]) + config_[site3] + config_[site4]) * gsoa;
    double term2 = (-config_[site1] + config_[site2] + config_[site3] + config_[site4]) * gsob;
    double term3 = (config_[site1] + config_[site2] + (-config_[site3]) + config_[site4]) * gsoc;
    using complexT = std::complex<double>;
    complexT fA(1.0), fB(exp(complexT(0.0, 2 * M_PI / 3))), fC(exp(complexT(0.0, 4 * M_PI / 3)));
    return term1 * fA + term2 * fB + term3 * fC;
  }

  double LocalComplexOrderParameterCorrelation(
      const size_t lx, const size_t ly,
      const size_t center_site1, const size_t center_site2,
      LocalDOFT &gsoa, LocalDOFT &gsob, LocalDOFT &gsoc
  ) {
    std::complex<double> m1 = LocalComplexOrderParameter(lx, ly, center_site1, gsoa, gsob, gsoc);
    std::complex<double> m2 = LocalComplexOrderParameter(lx, ly, center_site2, gsoa, gsob, gsoc);
    return m1.real() * m2.real() + m1.imag() * m2.imag();
  }

  double LocalComplexOrderParameterCorrelationLover2(
      const size_t lx, const size_t ly,
      LocalDOFT &gsoa, LocalDOFT &gsob, LocalDOFT &gsoc
  ) {
    double correlation(0.0);
    for (size_t site = 0; site < N_ / 2; site++) { //assum mod(ly,8)==0
      size_t y = site / (2 * lx);
      size_t x = site % (2 * lx);
      if ((y % 4 == 1 && x % 4 == 2) || (y % 4 == 3 && x % 4 == 0)) {
        correlation += LocalComplexOrderParameterCorrelation(lx, ly, site, site + N_ / 2, gsoa, gsob, gsoc);
      }
    }
    correlation = correlation / (N_ / 8);
    return correlation;
  }

  double LocalComplexOrderParameterCorrelationLover4(
      const size_t lx, const size_t ly,
      LocalDOFT &gsoa, LocalDOFT &gsob, LocalDOFT &gsoc
  ) {
    double correlation(0.0);
    for (size_t site = 0; site < N_ / 2; site++) { //assum mod(ly,16)==0
      size_t y = site / (2 * lx);
      size_t x = site % (2 * lx);
      if ((y % 4 == 1 && x % 4 == 2) || (y % 4 == 3 && x % 4 == 0)) {
        correlation += LocalComplexOrderParameterCorrelation(lx, ly, site, site + N_ / 4, gsoa, gsob, gsoc);
      }
    }
    correlation = correlation / (N_ / 8);
    return correlation;
  }

  std::array<double, DimDof> SpinCorrelationHoneycombLover2(const size_t lx) const {
    std::array<double, DimDof> correlations;
    for (size_t alpha = 0; alpha < DimDof; alpha++) {
      correlations[alpha] = 0.0;
    }
    for (size_t site = 0; site < N_ / 2; site++) {
      size_t y = site / (2 * lx);
      size_t x = site % (2 * lx);
      if ((x + y) % 2 == 0) { // A sublattice
        for (size_t alpha = 0; alpha < DimDof; alpha++) {
          correlations[alpha] += config_[site].GetCoor(alpha) * config_[site + N_ / 2].GetCoor(alpha);
        }
      }
    }

    for (size_t alpha = 0; alpha < DimDof; alpha++) {
      correlations[alpha] /= (N_ / 4);
    }
    return correlations;
  }

  std::array<double, DimDof> SpinCorrelationHoneycombLover4(const size_t lx) const {
    std::array<double, DimDof> correlations;
    for (size_t alpha = 0; alpha < DimDof; alpha++) {
      correlations[alpha] = 0.0;
    }
    for (size_t site = 0; site < N_; site++) {
      size_t y = site / (2 * lx);
      size_t x = site % (2 * lx);
      if ((x + y) % 2 == 0) { // A sublattice
        for (size_t alpha = 0; alpha < DimDof; alpha++) {
          size_t site2 = (site + N_ / 4) % N_;
          correlations[alpha] += config_[site].GetCoor(alpha) * config_[site2].GetCoor(alpha);
        }
      }
    }
    for (size_t alpha = 0; alpha < DimDof; alpha++) {
      correlations[alpha] /= (N_ / 2);
    }
    return correlations;
  }

  double HoneycombZigzagComplexOrderParameter(const size_t lx) const {
    using complexT = std::complex<double>;
    complexT m(0.0), fA(1.0 / 2.0), fB(exp(complexT(0.0, 2 * M_PI / 3)) / 2.0),
        fC(exp(complexT(0.0, 4 * M_PI / 3)) / 2.0);
    for (size_t site = 0; site < N_; site++) {
      size_t y = site / (2 * lx);
      size_t x = site % (2 * lx);
      if (y % 2 == 0) {//A, C, (-B)
        if ((y / 2) % 2 == 0) {
          if (x % 4 == 1) {//C
            m += config_[site].GetCoor(0) * fC;
          } else if (x % 4 == 3) {//A
            m += config_[site].GetCoor(0) * fA;
          } else if (x % 4 == 2) {// -B
            m = m - config_[site].GetCoor(0) * fB;
          }
        } else {
          if (x % 4 == 1) {//A
            m += config_[site].GetCoor(0) * fA;
          } else if (x % 4 == 3) {//C
            m += config_[site].GetCoor(0) * fC;
          } else if (x % 4 == 0) {// -B
            m = m - config_[site].GetCoor(0) * fB;
          }
        }
      } else { //B, (-A), (-C)
        if (((y - 1) / 2) % 2 == 0) {
          if (x % 4 == 0) {//B
            m += config_[site].GetCoor(0) * fB;
          } else if (x % 4 == 1) {//-C
            m -= config_[site].GetCoor(0) * fC;
          } else if (x % 4 == 3) {//-A
            m -= config_[site].GetCoor(0) * fA;
          }
        } else {
          if (x % 4 == 2) {
            m += config_[site].GetCoor(0) * fB;
          } else if (x % 4 == 1) {//-A
            m -= config_[site].GetCoor(0) * fA;
          } else if (x % 4 == 3) {//-C
            m -= config_[site].GetCoor(0) * fC;
          }
        }
      }
    }
    return abs(m * m) / N_ * 4 / N_ * 4;
  }

  std::array<double, DimDof> HoneycombZigzagAntiferromagneticMagnetization(const size_t lx) const {
    std::array<double, DimDof> m_q;
    for (size_t alpha = 0; alpha < DimDof; alpha++) {
      m_q[alpha] = 0.0;
    }
    for (size_t site = 0; site < N_; site++) {
      size_t y = site / (2 * lx);

      if (y % 2 == 0) {
        for (size_t alpha = 0; alpha < DimDof; alpha++) {
          m_q[alpha] += config_[site].GetCoor(alpha);
        }
      } else {
        for (size_t alpha = 0; alpha < DimDof; alpha++) {
          m_q[alpha] = m_q[alpha] - config_[site].GetCoor(alpha);
        }
      }
    }
//    for (size_t alpha = 0; alpha < DimDof; alpha++) {
//      m_q[alpha] *= m_q[alpha];
//    }
    return m_q;
  }

  /// another direction of zig-zag order
  std::array<double, DimDof> HoneycombZigzagAntiferromagneticMagnetization2(const size_t lx) const {
    std::array<double, DimDof> m_q;
    for (size_t alpha = 0; alpha < DimDof; alpha++) {
      m_q[alpha] = 0.0;
    }
    for (size_t site = 0; site < N_; site++) {
      size_t y = site / (2 * lx);
      size_t x = site % (2 * lx);
      size_t x_mod = x % 4;

      if ((y + 1) % 4 == x_mod || (y + 2) % 4 == x_mod) {
        for (size_t alpha = 0; alpha < DimDof; alpha++) {
          m_q[alpha] -= config_[site].GetCoor(alpha);
        }
      } else {
        for (size_t alpha = 0; alpha < DimDof; alpha++) {
          m_q[alpha] += config_[site].GetCoor(alpha);
        }
      }
    }
    return m_q;
  }

  /// another direction of zig-zag order
  std::array<double, DimDof> HoneycombZigzagAntiferromagneticMagnetization3(const size_t lx) const {
    std::array<double, DimDof> m_q;
    for (size_t alpha = 0; alpha < DimDof; alpha++) {
      m_q[alpha] = 0.0;
    }
    for (size_t site = 0; site < N_; site++) {
      size_t y = site / (2 * lx);
      size_t x = site % (2 * lx);
      size_t y_mod = y % 4;
      size_t x_mod = x % 4;

      if ((4 - y_mod) % 4 == x_mod || (5 - y_mod) % 4 == x_mod) {
        for (size_t alpha = 0; alpha < DimDof; alpha++) {
          m_q[alpha] += config_[site].GetCoor(alpha);
        }
      } else {
        for (size_t alpha = 0; alpha < DimDof; alpha++) {
          m_q[alpha] -= config_[site].GetCoor(alpha);
        }
      }
    }
    return m_q;
  }

  void PrintSign(const LocalDOFT &n, const size_t lx) const {
    for (size_t i = 0; i < N_; i++) {
      if (i % lx == 0) {
        std::cout << "\n";
      }
      if (config_[i] * n > 0) {
        std::cout << "+";
      } else {
        std::cout << "-";
      }
    }
    std::cout << std::endl;
  }

  void PrintProjection(const LocalDOFT &n, const size_t lx) const {
    for (size_t i = 0; i < N_; i++) {
      if (i % lx == 0) {
        std::cout << "\n";
      }
      double projection = config_[i] * n;
      if (projection > 0) {
        std::cout << "+" << std::fixed << std::setprecision(1) << projection << " ";
      } else {
        std::cout << std::fixed << std::setprecision(1) << projection << " ";
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
              const LocalDOFT &n) const {
    double DeltaE2 = 2 * Keff * (config_[site1] * n) * (config_[site2] * n);
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

  /// One more argument. Generate the random number outsite
  bool Active(const size_t site1,
              const size_t site2,
              const double Keff,
              const LocalDOFT &n,
              const double random_number) const {
    double DeltaE2 = 2 * Keff * (config_[site1] * n) * (config_[site2] * n);
    if (DeltaE2 <= 0) {
      return false;
    } else {
      double P = exp(-DeltaE2);
      if (random_number > P) {
        return true;
      } else {
        return false;
      }
    }
  }

  //TODO: optimize
  template<size_t InteractionRange>
  double EnergyRelevantToSite(
      const size_t site,
      const LatticeLink<DimDof, InteractionRange> &lattice_link
  ) const {
    auto &coupling_structures = lattice_link.coupling_structures;
    auto &links = lattice_link.links;
    auto &spin = config_[site];
    double e0 = spin * (coupling_structures[0] * spin);

    double e2 = 0.0;
    for (size_t inter_range = 1; inter_range < InteractionRange; inter_range++) {
      auto &J = coupling_structures[inter_range];
      auto &link_set = links[inter_range];
      assert(N_ == link_set.size());
      auto &site2_set = link_set[site];
      for (auto site2: site2_set) {
        auto &spin2 = config_[site2];
        e2 += spin * (J * spin2);
      }
    }
    return -(e0 + e2);
  }
  /// same with above function but consider the hp break Zp symmetry field, usually in XY model
  template<size_t InteractionRange>
  double EnergyRelevantToSite(
      const size_t site,
      const LatticeLink<DimDof, InteractionRange> &lattice_link,
      const double hp
  ) const {
    double e0 = EnergyRelevantToSite(site, lattice_link);
    auto &spin = config_[site];
    return e0 + hp * spin.Cos6Theta();
  }

  /// the case we assume the spin in the site is replaced by another spin
  template<size_t InteractionRange>
  double EnergyRelevantToSite(
      const size_t site,
      const LatticeLink<DimDof, InteractionRange> &lattice_link,
      const LocalDOFT &spin
  ) const {
    auto &coupling_structures = lattice_link.coupling_structures;
    auto &links = lattice_link.links;
    double e0 = spin * (coupling_structures[0] * spin);

    double e2 = 0.0;
    for (size_t inter_range = 1; inter_range < InteractionRange; inter_range++) {
      auto &J = coupling_structures[inter_range];
      auto &link_set = links[inter_range];
      assert(N_ == link_set.size());
      auto &site2_set = link_set[site];
      for (auto site2: site2_set) {
        auto &spin2 = config_[site2];
        e2 += spin * (J * spin2);
      }
    }
    return -(e0 + e2);
  }
  /// The same with above function but consider the hp field
  template<size_t InteractionRange>
  double EnergyRelevantToSite(
      const size_t site,
      const LatticeLink<DimDof, InteractionRange> &lattice_link,
      const LocalDOFT &spin,
      const double hp
  ) const {
    double e0 = EnergyRelevantToSite(site, lattice_link, spin);
    return e0 + hp * spin.Cos6Theta();
  }

  template<size_t InteractionRange>
  double EnergyDifferenceFlipSite(
      const size_t site,
      const LocalDOFT &fliped_spin,
      const LatticeLink<DimDof, InteractionRange> &lattice_link
  ) const {
    return -EnergyRelevantToSite(site, lattice_link) + EnergyRelevantToSite(site, lattice_link, fliped_spin);
  }

  template<size_t InteractionRange>
  double EnergyDifferenceFlipSite(
      const size_t site,
      const LocalDOFT &fliped_spin,
      const LatticeLink<DimDof, InteractionRange> &lattice_link,
      const double hp
  ) const {
    return -EnergyRelevantToSite(site, lattice_link, hp) + EnergyRelevantToSite(site, lattice_link, fliped_spin, hp);
  }

  double EnergyDifferenceOfSIAFlipCluster(
      const LocalDOFT &n,
      const SWCluster &cluster,
      const CouplingStructure<DimDof> &D
  ) {
    double sia(0.0);
    for (size_t i: cluster) {
      LocalDOF<DimDof> &spin = config_[i];
      auto [para, perp] = spin.DecompSpin(n);
      sia += para * (D * perp);
    }
    return 4 * sia;
  }

  double EnergyDifferenceOfHpFlipCluster(
      const LocalDOFT &n,
      const SWCluster &cluster,
      const double hp
  ) {
    double ehp(0.0);//energy of hp field
    for (size_t i: cluster) {
      LocalDOF<DimDof> flipped_spin = config_[i];
      flipped_spin.Flip(n);
      ehp += flipped_spin.Cos6Theta() - config_[i].Cos6Theta();
    }
    return ehp;
  }

  template<size_t InteractionRange>
  double EnergyDifferenceFlipCluster(
      const LocalDOFT &n,
      const SWCluster &cluster,
      const LatticeLink<DimDof, InteractionRange> &lattice_link,
      const double hp = 0.0
  ) {
    ComponentDecompose_(n);
    auto &para_config = parallel_components_;
    auto &perp_config = perpendicular_components_;
    double h2(0.0), sia(0.0); //H2, SIA
    auto &D = lattice_link.coupling_structures[0];
    for (size_t i: cluster) {
      sia += (para_config[i] * (D * perp_config[i]));
    }
    sia *= 2;
    for (size_t inter_type = 1; inter_type < InteractionRange; inter_type++) {
      auto &J = lattice_link.coupling_structures[inter_type];
      if (!J.IsIsometry()) {
        for (size_t i: cluster) {
          const std::vector<size_t> &j_set = lattice_link.GetSitesLinkedTo(i, inter_type);

          LocalDOFT sum_j_perp;
          sum_j_perp.Zero();
          for (auto j: j_set) {
            sum_j_perp += perp_config[j];
          }
          h2 += sum_j_perp * (J * para_config[i]);
        }
      }
    }

    double ehp(0.0);//energy of hp field
    if (hp != 0.0) {
      for (size_t i: cluster) {
        LocalDOF<DimDof> flipped_spin = config_[i];
        flipped_spin.Flip(n);
        ehp += flipped_spin.Cos6Theta() - config_[i].Cos6Theta();
      }
    }
    return 2 * (h2 + sia) + ehp;
  }

  /// interaction support upto J3, and Heisenberg or XY like (no anisotropy except SIA terms)
  template<size_t NumOfCouplingType>
  double StiffnessIsotropy(
      const LatticeLink<DimDof, NumOfCouplingType> &lattice_link,
      const double beta
  ) const {
    double a_term(0.0), b_term(0.0);
    for (size_t inter_range = 1; inter_range < NumOfCouplingType; inter_range++) {
      double j = lattice_link.coupling_structures[inter_range].coupling_const[0][0];
      auto &oriented_link = lattice_link.oriented_links[inter_range];
      auto &proj2y_distance = lattice_link.proj2y_link_length[inter_range];
      for (size_t i = 0; i < oriented_link.size(); i++) {
        auto [site1, site2] = oriented_link[i];
        double distance = proj2y_distance[i];
        a_term += j * Project2DInnerProduct(config_[site1], config_[site2]) * distance * distance;
        b_term += j * ZCompOfCrossProduct(config_[site1], config_[site2]) * distance;
      }
    }
    double
        stiffness = (a_term - beta * b_term * b_term) / lattice_link.GetLatticeArea();
    return stiffness;
  }

  template<size_t NumOfCouplingType>
  double StiffnessAnisotropy(
      const LatticeLink<DimDof, NumOfCouplingType> &lattice_link,
      const double beta,
      double &b_term // output
  ) const {
    double a_term(0.0);
    b_term = 0.0;
    for (size_t inter_range = 1; inter_range < NumOfCouplingType; inter_range++) {
      auto &J = lattice_link.coupling_structures[inter_range];
      auto &oriented_link = lattice_link.oriented_links[inter_range];
      auto &proj2y_distance = lattice_link.proj2y_link_length[inter_range];
      for (size_t i = 0; i < oriented_link.size(); i++) {
        auto [site1, site2] = oriented_link[i];
        double distance = proj2y_distance[i];
        LocalDOF<DimDof> projected_spin2 = config_[site2].Project2XYPlane();
        a_term += config_[site1] * (J * projected_spin2) * distance * distance;
        LocalDOF<DimDof> orthognoal_of_projected_spin2 = config_[site2].OrthogonalProject2XYPlane();
        b_term += config_[site1] * (J * orthognoal_of_projected_spin2) * distance;
      }
    }
    double
        stiffness = (a_term - beta * b_term * b_term) / lattice_link.GetLatticeArea();
    return stiffness;
  }

  template<size_t NumOfCouplingType>
  double StiffnessSquare(
      const LatticeLink<DimDof, NumOfCouplingType> &lattice_link,
      const double beta
  ) const {
    double inner_product_sum[NumOfCouplingType - 1];
    double cross_product_sum[NumOfCouplingType - 1];
    for (size_t inter_range = 1; inter_range < NumOfCouplingType; inter_range++) {
      auto &oriented_link = lattice_link.oriented_links[inter_range];
      double sum_a(0.0), sum_b(0.0);
      for (auto [site1, site2]: oriented_link) {
        sum_a += Project2DInnerProduct(config_[site1], config_[site2]);

#ifndef NDEBUG
        double z_comp_of_cross_product = ZCompOfCrossProduct(config_[site1], config_[site2]);
        if (DimDof == 2) {
          double sin_value = ::SinDiff(config_[site1], config_[site2]);
          assert(sin_value = -z_comp_of_cross_product);
        }
        sum_b += z_comp_of_cross_product;
#else
        sum_b += ZCompOfCrossProduct(config_[site1], config_[site2]);
#endif
      }
      inner_product_sum[inter_range - 1] = sum_a;
      cross_product_sum[inter_range - 1] = sum_b;
    }
    double coefficient_by_lattice[3] = {0.5, 1, 2.0}; // the only line different with Honeycomb
    double stiffness = 0.0;
    for (size_t i = 0; i < NumOfCouplingType - 1; i++) {
      const size_t inter_range = i + 1;
      double j = lattice_link.coupling_structures[inter_range].coupling_const[0][0];
      stiffness += coefficient_by_lattice[i]
          * (j * inner_product_sum[i] - beta * j * j * cross_product_sum[i] * cross_product_sum[i])
          / lattice_link.GetLatticeArea();
    }
    return stiffness;
  }

 private:

//  std::pair<SpinConfigT, SpinConfigT> ComponentDecompose_(
//      const LocalDOFT &n
//  ) const {
//    std::pair<SpinConfigT, SpinConfigT> res;
//    res.first.reserve(N_);
//    res.second.reserve(N_);
//    for (size_t i = 0; i < N_; i++) {
//      auto [para, perp] = config_[i].DecompSpin(n);
//      res.first[i] = para;
//      res.second[i] = perp;
//    }
//    return res;
//  }

  void ComponentDecompose_(
      const LocalDOFT &n
  ) {
    for (size_t i = 0; i < N_; i++) {
      auto [para, perp] = config_[i].DecompSpin(n);
      parallel_components_[i] = para;
      perpendicular_components_[i] = perp;
    }
  }

  size_t N_; //lattice size
  SpinConfigT config_;

  SpinConfigT parallel_components_;
  SpinConfigT perpendicular_components_;
};
#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_LATTICE_CONFIG_H
