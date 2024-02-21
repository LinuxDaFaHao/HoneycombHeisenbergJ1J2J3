//
// Created by Hao-Xin on 2022/5/25.
//

#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_XYMCEXECUTOR_H
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_XYMCEXECUTOR_H

#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
#include "qlten/framework/bases/executor.h"
#include "qlten/utility/timer.h"
#include "swcluster.h"
#include "lattice_link.h"
#include "lattice_config.h"
#include "mpi.h"
#include "common.h"

const size_t dim_dof_xy = 2;

struct MCParams {
  size_t sweeps;
  size_t sample_interval; // for how many times of wolf update between two samples. 0 means sweep lattice
  size_t print_interval;
  std::string filename_postfix;
};

template<size_t dim_dof_xy, size_t NumOfCouplingType>
struct XYPhysParams {
  XYPhysParams(const std::string &geometry,
               size_t lx,
               size_t ly,
               size_t n,
               double beta,
               const std::array<CouplingStructure<dim_dof_xy>, NumOfCouplingType> &coupling_structures,
               const double hp)
      : geometry(geometry), Lx(lx), Ly(ly), N(n), beta(beta), coupling_structures(coupling_structures),
        hp(hp) {}

  XYPhysParams(const std::string &geometry,
               size_t lx,
               size_t ly,
               double beta,
               const std::array<CouplingStructure<dim_dof_xy>, NumOfCouplingType> &coupling_structures,
               const double hp)
      : geometry(geometry), Lx(lx), Ly(ly), beta(beta), coupling_structures(coupling_structures),
        hp(hp) {
    if (geometry == "Honeycomb" || geometry == "Honeycomb3") {
      N = lx * ly * 2;
    } else if (geometry == "Square") {
      N = lx * ly;
    } else if (geometry == "SVW") {
      N = lx * ly;
    }
  }

  std::string geometry; //Honeycomb
  size_t Lx;
  size_t Ly;
  size_t N;
  double beta; // 1/T
  std::array<CouplingStructure<dim_dof_xy>, NumOfCouplingType> coupling_structures;
  double hp; // Zp broken field, we assume p = 6 now.
};

template<size_t NumOfCouplingType>
class XYMCExecutor : public qlten::Executor {
  using LocalDOFT = LocalDOF<dim_dof_xy>;
 public:
  XYMCExecutor(const MCParams &,
               const XYPhysParams<dim_dof_xy, NumOfCouplingType> &);
  void Execute() override;

  const std::vector<double> &GetEnergyData() const {
    return energy_;
  }

  const std::array<std::vector<double>, dim_dof_xy> &GetMagneticSusceptibility() const {
    return sum_spin_;
  }

  ~XYMCExecutor() {
    delete plattice_link_;
  }

  MCParams mc_params;
  XYPhysParams<dim_dof_xy, NumOfCouplingType> phys_params;

 private:

  void WolfSweep_();
  void WolfRandomSweep_(size_t sample_interval);
  void WolfGrowCluster_(const LocalDOFT &, const size_t);
  bool WolfFlipCluster_(const LocalDOFT &);

  void MetropolisSweep_();
  bool MetropolisSingleStep_(const size_t site);

  void Measure_(size_t sweep);
  void StatisticAndDumpData_();

  size_t geometry_id_; // 0: square, 1: honeycomb, 2: SVW
  bool interaction_isotropy_; // not consider SIA
  bool isotropy_;
  LatticeConfig<dim_dof_xy> config_;
  SWCluster cluster_;
  LatticeLink<dim_dof_xy, NumOfCouplingType> *plattice_link_;

  /// result data
  std::vector<double> energy_;
  std::array<std::vector<double>, dim_dof_xy> sum_spin_;
  std::vector<double> correlation_lover2_;
  std::vector<double> correlation_lover4_;
  std::vector<double> stiffness_;         // if U(1) symmetric model

  std::vector<double> res_; //average value, using the latter half part data
  // persite energy, specific heat, magnetic susceptibility, stiffness,
};

template<size_t NumOfCouplingType>
XYMCExecutor<NumOfCouplingType>::XYMCExecutor(const MCParams &mc_params,
                                              const XYPhysParams<dim_dof_xy, NumOfCouplingType> &phys_params) :
    qlten::Executor(),
    mc_params(mc_params),
    phys_params(phys_params),
    config_(phys_params.N),
    cluster_() {
  if (phys_params.geometry == "Honeycomb") {
    plattice_link_ = new HoneyCombTorusLatticeLink(phys_params.Lx, phys_params.Ly, phys_params.coupling_structures);
    geometry_id_ = 1;
  } else if (phys_params.geometry == "Honeycomb3") {
    plattice_link_ = new HoneyComb3TorusLatticeLink(phys_params.Lx, phys_params.Ly, phys_params.coupling_structures);
    geometry_id_ = 3;
  } else if (phys_params.geometry == "Square") {
    plattice_link_ = new SquareTorusLatticeLink(phys_params.Lx, phys_params.Ly, phys_params.coupling_structures);
    geometry_id_ = 0;
  } else if (phys_params.geometry == "SVW") {
    plattice_link_ = new SVWLatticeLink(phys_params.N, phys_params.coupling_structures);
    geometry_id_ = 2;
  } else {
    std::cout << "do not support now. " << std::endl;
    exit(0);
  }

  interaction_isotropy_ = true;
  for (size_t i = 1; i < NumOfCouplingType; i++) {
    interaction_isotropy_ = interaction_isotropy_ && phys_params.coupling_structures[i].IsIsometry();
  }

  isotropy_ = interaction_isotropy_ && phys_params.coupling_structures[0].IsIsometry() && (phys_params.hp == 0.0);

  config_.Random();
  size_t sweeps = mc_params.sweeps;
  energy_.reserve(sweeps);

  if (dim_dof_xy >= 2 && interaction_isotropy_) {
    stiffness_.reserve(sweeps);
  }

  for (size_t i = 0; i < dim_dof_xy; i++) {
    sum_spin_[i].reserve(sweeps);
  }

  correlation_lover2_.reserve(sweeps);
  correlation_lover4_.reserve(sweeps);

  SetStatus(qlten::INITED);
}

template<size_t NumOfCouplingType>
void XYMCExecutor<NumOfCouplingType>::Execute() {
  SetStatus(qlten::EXEING);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  qlten::Timer wolf_mc_execute_timer("wolf_mc_execute");
  for (size_t sweep = 0; sweep < mc_params.sweeps; sweep++) {
    MetropolisSweep_();

    if (mc_params.sample_interval != phys_params.N) {
      WolfRandomSweep_(mc_params.sample_interval);
    } else {
      WolfSweep_();
    }

    Measure_(sweep);

    if (world_rank == 0 && sweep % mc_params.print_interval == 0) {
      double execute_time = wolf_mc_execute_timer.Elapsed();
      std::cout << "[ sweep = " << sweep << " ]"
                << "time = " << execute_time
                << std::endl;
    }
  }

  StatisticAndDumpData_();

  SetStatus(qlten::FINISH);
}

template<size_t NumOfCouplingType>
void XYMCExecutor<NumOfCouplingType>::WolfSweep_() {
  LocalDOF<dim_dof_xy> axis;
  for (size_t starting_site = 0; starting_site < phys_params.N; starting_site++) {
    axis.Random();
    WolfGrowCluster_(axis, starting_site);
    WolfFlipCluster_(axis);
  }
}

template<size_t NumOfCouplingType>
void XYMCExecutor<NumOfCouplingType>::WolfRandomSweep_(size_t sample_interval) {
  LocalDOF<dim_dof_xy> axis;
  std::uniform_int_distribution<int> u_int(0, phys_params.N - 1);
  for (size_t i = 0; i < sample_interval; i++) {
    axis.Random();
    WolfGrowCluster_(axis, u_int(random_engine));
    WolfFlipCluster_(axis);
  }
}

template<size_t NumOfCouplingType>
void XYMCExecutor<NumOfCouplingType>::WolfGrowCluster_(const LocalDOF<dim_dof_xy> &axis,
                                                       const size_t starting_site) {
  cluster_ = {starting_site};//reserve?
  cluster_.reserve(phys_params.N);
  std::vector<bool> in_cluster(phys_params.N, false);
  in_cluster[starting_site] = true;
// config_.PrintSign(axis, phys_params.Lx);
// Note here the effective single ion anisotropy is not important.
  double K_eff[NumOfCouplingType - 1];
  for (size_t inter_type = 1; inter_type < NumOfCouplingType; inter_type++) {
    double eff_coupling = axis * (phys_params.coupling_structures[inter_type] * axis);
    K_eff[inter_type - 1] = eff_coupling * phys_params.beta;
  }

  for (size_t j = 0; j < cluster_.size(); j++) {
    for (size_t inter_type = 1; inter_type < NumOfCouplingType; inter_type++) {
      const std::vector<size_t> &site2_set = plattice_link_->GetSitesLinkedTo(cluster_[j], inter_type);
      for (auto site2: site2_set) {
        if (!in_cluster[site2] && config_.Active(cluster_[j], site2, K_eff[inter_type - 1], axis)) {
          in_cluster[site2] = true;
          cluster_.push_back(site2);
        }
      }
    }
  }
}

/**
 *
 * @param axis
 * @return flipped or not
 */
template<size_t NumOfCouplingType>
bool XYMCExecutor<NumOfCouplingType>::WolfFlipCluster_(const LocalDOF<dim_dof_xy> &axis) {

  if (isotropy_) {
    config_.Flip(axis, cluster_);
    return true;
  }
  double delta_e;
  if (interaction_isotropy_) {
    delta_e = config_.EnergyDifferenceOfSIAFlipCluster(axis, cluster_, plattice_link_->coupling_structures[0])
        + config_.EnergyDifferenceOfHpFlipCluster(axis, cluster_, phys_params.hp);
  } else {
    delta_e = config_.EnergyDifferenceFlipCluster(axis, cluster_, *plattice_link_, phys_params.hp);
  }
  if (delta_e <= 0.0) {
    config_.Flip(axis, cluster_);
    return true;
  } else {
    std::uniform_real_distribution<double> u(0, 1);
    if (u(random_engine) <= std::exp(-phys_params.beta * delta_e)) {
      config_.Flip(axis, cluster_);
      return true;
    } else {
      return false;
    }
  }
}

template<size_t NumOfCouplingType>
void XYMCExecutor<NumOfCouplingType>::MetropolisSweep_() {
  for (size_t i = 0; i < phys_params.N; i++) {
    MetropolisSingleStep_(i);
  }
}

template<size_t NumOfCouplingType>
bool XYMCExecutor<NumOfCouplingType>::MetropolisSingleStep_(const size_t site) {
  LocalDOF<dim_dof_xy> flipto;
  flipto.Random();
  double delta_e =
      config_.template EnergyDifferenceFlipSite<NumOfCouplingType>(site, flipto, *plattice_link_, phys_params.hp);
  if (delta_e <= 0.0) {
    config_.SetSiteSpin(site, flipto);
    return true;
  } else {
    std::uniform_real_distribution<double> u(0, 1);
    if (u(random_engine) <= std::exp(-phys_params.beta * delta_e)) {
      config_.SetSiteSpin(site, flipto);
      return true;
    } else {
      return false;
    }
  }
}


template<size_t NumOfCouplingType>
void XYMCExecutor<NumOfCouplingType>::Measure_(size_t sweep) {
//  config_.PrintProjection(LocalDOF<dim_dof_xy>({1.0,0,0}), phys_params.Lx * 2);
//  std::cout << "\n" << std::endl;

  energy_.push_back(config_.Energy(*plattice_link_, phys_params.hp));
  for (size_t i = 0; i < dim_dof_xy; i++) {
    sum_spin_[i].push_back(config_.SumOverComponent(i));
  }

  if (dim_dof_xy >= 2 && interaction_isotropy_ && (phys_params.hp != 0.0)) {
    stiffness_.push_back(config_.template StiffnessIsotropy<NumOfCouplingType>(*plattice_link_,
                                                                               phys_params.beta));
  }
}


template<size_t NumOfCouplingType>
void XYMCExecutor<NumOfCouplingType>::StatisticAndDumpData_() {
  qlten::Timer dump_data_timer("dump_data");
  DumpData("energy" + mc_params.filename_postfix, energy_);
  for (size_t i = 0; i < dim_dof_xy; i++) {
    DumpData("sum_spin"
                 + std::to_string(i)
                 + mc_params.filename_postfix,
             sum_spin_[i]);
  }

  if (dim_dof_xy >= 2 && interaction_isotropy_) {
    DumpData("stiffness" + mc_params.filename_postfix, stiffness_);
  }
  dump_data_timer.PrintElapsed();

  qlten::Timer statistic_timer("statistic");
  // energy
  auto half_data_of_energy = std::vector(energy_.begin() + mc_params.sweeps / 2, energy_.end());
  res_.push_back(Mean(half_data_of_energy) / phys_params.N);
  // specific heat
  res_.push_back(
      Variance(half_data_of_energy, res_[0] * phys_params.N) * phys_params.beta * phys_params.beta / phys_params.N);
  // magnetic susceptibility
  for (size_t i = 0; i < dim_dof_xy; i++) {
    auto half_data_of_sum_spin = std::vector(sum_spin_[i].begin() + mc_params.sweeps / 2, sum_spin_[i].end());
    res_.push_back(Variance(half_data_of_sum_spin) / phys_params.N);
  }

  if (dim_dof_xy > 1 && interaction_isotropy_ && (phys_params.hp == 0.0)) {
    // stiffness
    res_.push_back(Mean(std::vector(stiffness_.begin() + mc_params.sweeps / 2, stiffness_.end())));
  }

  //FM magnetization
  size_t sweeps = mc_params.sweeps;
  std::vector<double> magnetization_square(mc_params.sweeps / 2, 0.0), magnetization_quartic(mc_params.sweeps / 2, 0.0);
  for (size_t j = 0; j < sweeps / 2; j++) {
    magnetization_square[j] = (sum_spin_[0][j + sweeps / 2] * sum_spin_[0][j + sweeps / 2]
        + sum_spin_[1][j + sweeps / 2] * sum_spin_[1][j + sweeps / 2]) / phys_params.N / phys_params.N;
    magnetization_quartic[j] = magnetization_square[j] * magnetization_square[j];
  }
  double M2 = Mean(magnetization_square);
  res_.push_back(M2);
  //binder ratio
  double M4 = Mean(magnetization_quartic);
  res_.push_back(M4 / M2 / M2);

  DumpData("summary" + mc_params.filename_postfix, res_);
  statistic_timer.
      PrintElapsed();
}

#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_XYMCEXECUTOR_H
