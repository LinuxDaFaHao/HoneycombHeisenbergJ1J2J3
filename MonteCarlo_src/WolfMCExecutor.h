//
// Created by Hao-Xin on 2022/5/25.
//

#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_WOLFMCEXECUTOR_H
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_WOLFMCEXECUTOR_H

#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
#include "qlten/framework/bases/executor.h"
#include "qlten/utility/timer.h"
#include "swcluster.h"
#include "lattice_link.h"
#include "lattice_config.h"
//#include "mpi.h"
#include "common.h"

struct MCParams {
  size_t sweeps;
  size_t sample_interval; // for how many times of wolf update between two samples. 0 means sweep lattice
  size_t print_interval;
  std::string filename_postfix;
};

template<size_t DimDof, size_t NumOfCouplingType>
struct PhysParams {
  PhysParams(const std::string &geometry,
             size_t lx,
             size_t ly,
             size_t n,
             double beta,
             const std::array<CouplingStructure<DimDof>, NumOfCouplingType> &coupling_structures)
      : geometry(geometry), Lx(lx), Ly(ly), N(n), beta(beta), coupling_structures(coupling_structures) {
    zigzag_gs_oritation_1.Zero();
    zigzag_gs_oritation_2.Zero();
    zigzag_gs_oritation_3.Zero();
  }

  PhysParams(const std::string &geometry,
             size_t lx,
             size_t ly,
             size_t n,
             double beta,
             const std::array<CouplingStructure<DimDof>, NumOfCouplingType> &coupling_structures,
             const LocalDOF<DimDof> &zigzag_gs_oritation_1,
             const LocalDOF<DimDof> &zigzag_gs_oritation_2,
             const LocalDOF<DimDof> &zigzag_gs_oritation_3)
      : geometry(geometry),
        Lx(lx),
        Ly(ly),
        N(n),
        beta(beta),
        coupling_structures(coupling_structures),
        zigzag_gs_oritation_1(zigzag_gs_oritation_1),
        zigzag_gs_oritation_2(zigzag_gs_oritation_2),
        zigzag_gs_oritation_3(zigzag_gs_oritation_3) {}

  PhysParams(const std::string &geometry,
             size_t lx,
             size_t ly,
             double beta,
             const std::array<CouplingStructure<DimDof>, NumOfCouplingType> &coupling_structures)
      : geometry(geometry), Lx(lx), Ly(ly), beta(beta), coupling_structures(coupling_structures) {
    if (geometry == "Honeycomb" || geometry == "Honeycomb3") {
      N = lx * ly * 2;
    } else if (geometry == "Square") {
      N = lx * ly;
    } else if (geometry == "SVW") {
      N = lx * ly;
    }
    zigzag_gs_oritation_1.Zero();
    zigzag_gs_oritation_2.Zero();
    zigzag_gs_oritation_3.Zero();
  }

  PhysParams(const std::string &geometry,
             size_t lx,
             size_t ly,
             double beta,
             const std::array<CouplingStructure<DimDof>, NumOfCouplingType> &coupling_structures,
             const LocalDOF<DimDof> &zigzag_gs_oritation_1,
             const LocalDOF<DimDof> &zigzag_gs_oritation_2,
             const LocalDOF<DimDof> &zigzag_gs_oritation_3)
      : geometry(geometry),
        Lx(lx),
        Ly(ly),
        beta(beta),
        coupling_structures(coupling_structures),
        zigzag_gs_oritation_1(zigzag_gs_oritation_1),
        zigzag_gs_oritation_2(zigzag_gs_oritation_2),
        zigzag_gs_oritation_3(zigzag_gs_oritation_3) {
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
  std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures;
  LocalDOF<DimDof> zigzag_gs_oritation_1;
  LocalDOF<DimDof> zigzag_gs_oritation_2;
  LocalDOF<DimDof> zigzag_gs_oritation_3;
};

template<size_t DimDof, size_t NumOfCouplingType>
class WolfMCExecutor : public qlten::Executor {
  using LocalDOFT = LocalDOF<DimDof>;
 public:
  WolfMCExecutor(const MCParams &,
                 const PhysParams<DimDof, NumOfCouplingType> &,
                 const size_t);
  void Execute() override;

  const std::vector<double> &GetEnergyData() const {
    return energy_;
  }

  const std::array<std::vector<double>, DimDof> &GetMagneticSusceptibility() const {
    return sum_spin_;
  }

  ~WolfMCExecutor() {
    delete plattice_link_;
  }

  MCParams mc_params;
  PhysParams<DimDof, NumOfCouplingType> phys_params;

 private:

  void WolfSweep_();
  void WolfRandomSweep_(size_t sample_interval);
  void WolfGrowCluster_(const LocalDOFT &, const size_t);
  bool WolfFlipCluster_(const LocalDOFT &);

  void MetropolisSweep_();
  bool MetropolisSingleStep_(const size_t site);

  void Measure_(size_t sweep);
  void StatisticAndDumpData_();

  int mpi_world_rank_;

  size_t geometry_id_; // 0: square, 1: honeycomb, 2: SVW
  bool interaction_isotropy_; // not consider SIA
  bool isotropy_;
  LatticeConfig<DimDof> config_;
  SWCluster cluster_;
  LatticeLink<DimDof, NumOfCouplingType> *plattice_link_;

  LocalDOF<DimDof> ground_state_oritation1_; //corresponding to the zig-zag order 1
  LocalDOF<DimDof> ground_state_oritation2_;
  LocalDOF<DimDof> ground_state_oritation3_;
  bool has_gs_data;

  size_t flip_times_;

  /// result data
  std::vector<double> energy_;
  std::array<std::vector<double>, DimDof> sum_spin_;
  std::array<std::vector<double>, DimDof> zigzag_af_magnetization_;
  std::array<std::vector<double>, DimDof> zigzag_af_magnetization2_;
  std::array<std::vector<double>, DimDof> zigzag_af_magnetization3_;
  std::vector<double> correlation_lover2_;
  std::vector<double> correlation_lover4_;
  std::vector<double> stiffness_;         // if U(1) symmetric model
//  std::vector<double> stiffness_correction_; // come from the first derivative term of free energy
  std::vector<double> complex_xy_order_parameter_real_;
  std::vector<double> complex_xy_order_parameter_imag_;

  std::vector<double> res_; //average value, using the latter half part data
  // persite energy, specific heat, magnetic susceptibility, stiffness,
};

template<size_t DimDof, size_t NumOfCouplingType>
WolfMCExecutor<DimDof, NumOfCouplingType>::WolfMCExecutor(const MCParams &mc_params,
                                                          const PhysParams<DimDof, NumOfCouplingType> &phys_params,
                                                          const size_t world_rank) :
    qlten::Executor(),
    mc_params(mc_params),
    phys_params(phys_params),
    mpi_world_rank_(world_rank),
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

  isotropy_ = interaction_isotropy_ && phys_params.coupling_structures[0].IsIsometry();

  config_.Random();
  size_t sweeps = mc_params.sweeps;
  energy_.reserve(sweeps);

  if (DimDof >= 2 && interaction_isotropy_) {
    stiffness_.reserve(sweeps);
//    stiffness_correction_ = std::vector<double>(sweeps, 0.0);
  }
  ground_state_oritation1_ = phys_params.zigzag_gs_oritation_1;
  ground_state_oritation2_ = phys_params.zigzag_gs_oritation_2;
  ground_state_oritation3_ = phys_params.zigzag_gs_oritation_3;
  if (ground_state_oritation1_.IsZero() && ground_state_oritation2_.IsZero() && ground_state_oritation3_.IsZero()) {
    has_gs_data = false;
  } else {
    has_gs_data = true;
  }

  flip_times_ = 0;

  for (size_t i = 0; i < DimDof; i++) {
    sum_spin_[i].reserve(sweeps);
    if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {
      zigzag_af_magnetization_[i].reserve(sweeps);
      zigzag_af_magnetization2_[i].reserve(sweeps);
      zigzag_af_magnetization3_[i].reserve(sweeps);
//      correlation_lover2_[i].reserve(sweeps);
//      correlation_lover4_[i].reserve(sweeps);
      if (has_gs_data) {
        complex_xy_order_parameter_real_.reserve(sweeps);
        complex_xy_order_parameter_imag_.reserve(sweeps);

        correlation_lover2_.reserve(sweeps);
        correlation_lover4_.reserve(sweeps);
      }
    }
  }

  SetStatus(qlten::INITED);
}

template<size_t DimDof, size_t NumOfCouplingType>
void WolfMCExecutor<DimDof, NumOfCouplingType>::Execute() {
  SetStatus(qlten::EXEING);
  qlten::Timer wolf_mc_execute_timer("wolf_mc_execute");
  for (size_t sweep = 0; sweep < mc_params.sweeps; sweep++) {
    MetropolisSweep_();

    if (mc_params.sample_interval != phys_params.N) {
      WolfRandomSweep_(mc_params.sample_interval);
    } else {
      WolfSweep_();
    }

    Measure_(sweep);

    if (mpi_world_rank_ == 0 && sweep % mc_params.print_interval == 0) {
      double execute_time = wolf_mc_execute_timer.Elapsed();
      std::cout << "[ sweep = " << sweep << " ]"
                << "time = " << execute_time << "\t"
                << "cluster size = " << cluster_.size()
                << std::endl;
    }
  }
  if (mpi_world_rank_ == 0) {
    std::cout << "Flip probability: " << double(flip_times_) / double(mc_params.sweeps * mc_params.sample_interval)
              << std::endl;
  }

  StatisticAndDumpData_();

  SetStatus(qlten::FINISH);
}

template<size_t DimDof, size_t NumOfCouplingType>
void WolfMCExecutor<DimDof, NumOfCouplingType>::WolfSweep_() {
  LocalDOF<DimDof> axis;
  for (size_t starting_site = 0; starting_site < phys_params.N; starting_site++) {
    axis.Random();
    WolfGrowCluster_(axis, starting_site);
    if (WolfFlipCluster_(axis)) flip_times_++;
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
void WolfMCExecutor<DimDof, NumOfCouplingType>::WolfRandomSweep_(size_t sample_interval) {
  LocalDOF<DimDof> axis;
  std::uniform_int_distribution<int> u_int(0, phys_params.N - 1);
  for (size_t i = 0; i < sample_interval; i++) {
    axis.Random();
    WolfGrowCluster_(axis, u_int(random_engine));
    if (WolfFlipCluster_(axis)) flip_times_++;
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
void WolfMCExecutor<DimDof, NumOfCouplingType>::WolfGrowCluster_(const LocalDOF<DimDof> &axis,
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
template<size_t DimDof, size_t NumOfCouplingType>
bool WolfMCExecutor<DimDof, NumOfCouplingType>::WolfFlipCluster_(const LocalDOF<DimDof> &axis) {

  if (isotropy_) {
    config_.Flip(axis, cluster_);
    return true;
  }
  double delta_e;
  if (interaction_isotropy_) {
    delta_e = config_.EnergyDifferenceOfSIAFlipCluster(axis, cluster_, plattice_link_->coupling_structures[0]);
  } else {
    delta_e = config_.EnergyDifferenceFlipCluster(axis, cluster_, *plattice_link_);
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

template<size_t DimDof, size_t NumOfCouplingType>
void WolfMCExecutor<DimDof, NumOfCouplingType>::MetropolisSweep_() {
  for (size_t i = 0; i < phys_params.N; i++) {
    MetropolisSingleStep_(i);
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
bool WolfMCExecutor<DimDof, NumOfCouplingType>::MetropolisSingleStep_(const size_t site) {
  LocalDOF<DimDof> flipto;
  flipto.Random();
  double delta_e = config_.template EnergyDifferenceFlipSite<NumOfCouplingType>(site, flipto, *plattice_link_);
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

template<size_t DimDof, size_t NumOfCouplingType>
void WolfMCExecutor<DimDof, NumOfCouplingType>::Measure_(size_t sweep) {
  // below two lines are used to debug
//  config_.PrintProjection(LocalDOF<DimDof>({1.0,0,0}), phys_params.Lx * 2);
//  std::cout << "\n" << std::endl;

  energy_.push_back(config_.Energy(*plattice_link_));
  /*
  for (size_t i = 0; i < DimDof; i++) {
    sum_spin_[i].push_back(config_.SumOverComponent(i));
  }
   */

  if (DimDof >= 2) {
    if (interaction_isotropy_) {
      stiffness_.push_back(config_.template StiffnessIsotropy<NumOfCouplingType>(*plattice_link_,
                                                                                 phys_params.beta));
    } else {
//      stiffness_.push_back(config_.template StiffnessAnisotropy<NumOfCouplingType>(*plattice_link_,
//                                                                                   phys_params.beta,
//                                                                                   stiffness_correction_[sweep]));
    }
  }

  if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {
    auto m1 = config_.HoneycombZigzagAntiferromagneticMagnetization(phys_params.Lx);
    auto m2 = config_.HoneycombZigzagAntiferromagneticMagnetization2(phys_params.Lx);
    auto m3 = config_.HoneycombZigzagAntiferromagneticMagnetization3(phys_params.Lx);

    for (size_t i = 0; i < DimDof; i++) {
      zigzag_af_magnetization_[i].push_back(m1[i]);
      zigzag_af_magnetization2_[i].push_back(m2[i]);
      zigzag_af_magnetization3_[i].push_back(m3[i]);
    }
    if (has_gs_data) {
      using complexT = std::complex<double>;
      complexT complex_order_parameter(0);
      complex_order_parameter += LocalDOF<DimDof>(m1) * ground_state_oritation1_;
      complex_order_parameter += LocalDOF<DimDof>(m2) * ground_state_oritation2_ * exp(complexT(0.0, 2 * M_PI / 3));
      complex_order_parameter += LocalDOF<DimDof>(m3) * ground_state_oritation3_ * exp(complexT(0.0, 4 * M_PI / 3));

      complex_xy_order_parameter_real_.push_back(complex_order_parameter.real() / phys_params.N);
      complex_xy_order_parameter_imag_.push_back(complex_order_parameter.imag() / phys_params.N);
    }
  }
//  correlation_lover2_.push_back(config_.LocalComplexOrderParameterCorrelationLover2(phys_params.Lx, phys_params.Ly,
//                                                                                    ground_state_oritation1_,
//                                                                                    ground_state_oritation2_,
//                                                                                    ground_state_oritation3_));
//  correlation_lover4_.push_back(config_.LocalComplexOrderParameterCorrelationLover4(phys_params.Lx, phys_params.Ly,
//                                                                                    ground_state_oritation1_,
//                                                                                    ground_state_oritation2_,
//                                                                                    ground_state_oritation3_));
}

template<size_t DimDof, size_t NumOfCouplingType>
void WolfMCExecutor<DimDof, NumOfCouplingType>::StatisticAndDumpData_() {
  qlten::Timer dump_data_timer("dump_data");
  DumpData("energy" + mc_params.filename_postfix, energy_);
  /*
   * for (size_t i = 0; i < DimDof; i++) {
    DumpData("sum_spin"
                 + std::to_string(i)
                 + mc_params.filename_postfix,
             sum_spin_[i]);
  }
   */
  if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {
    for (size_t i = 0; i < DimDof; i++) {
      DumpData("magnetization" + std::to_string(i) + mc_params.filename_postfix, zigzag_af_magnetization_[i]);
    }
    for (size_t i = 0; i < DimDof; i++) {
      DumpData("magnetization2" + std::to_string(i) + mc_params.filename_postfix, zigzag_af_magnetization2_[i]);
    }
    for (size_t i = 0; i < DimDof; i++) {
      DumpData("magnetization3" + std::to_string(i) + mc_params.filename_postfix, zigzag_af_magnetization3_[i]);
    }
    if (has_gs_data) {
      DumpData("xy_order_parameter_real" + mc_params.filename_postfix, complex_xy_order_parameter_real_);
      DumpData("xy_order_parameter_imag" + mc_params.filename_postfix, complex_xy_order_parameter_imag_);
    }
  }
  if (DimDof >= 2 && interaction_isotropy_) {
    DumpData("stiffness" + mc_params.filename_postfix, stiffness_);
  }

//  DumpData("xy_order_correlation2" + mc_params.filename_postfix, correlation_lover2_);
//  DumpData("xy_order_correlation2" + mc_params.filename_postfix, correlation_lover4_);
  dump_data_timer.PrintElapsed();

  qlten::Timer statistic_timer("statistic");
  size_t sweeps = mc_params.sweeps;
  // energy
  auto half_data_of_energy = std::vector(energy_.begin() + mc_params.sweeps / 2, energy_.end());
  res_.push_back(Mean(half_data_of_energy) / phys_params.N);
  // specific heat
  res_.push_back(
      Variance(half_data_of_energy, res_[0] * phys_params.N) * phys_params.beta * phys_params.beta / phys_params.N);
  // magnetic susceptibility
//  for (size_t i = 0; i < DimDof; i++) {
//    auto half_data_of_sum_spin = std::vector(sum_spin_[i].begin() + mc_params.sweeps / 2, sum_spin_[i].end());
//    res_.push_back(Variance(half_data_of_sum_spin) / phys_params.N);
//  }
  // anti-ferromagnetic susceptibility
  if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {
    for (size_t i = 0; i < DimDof; i++) {
      std::vector<double> valid_data_of_magnetization_abs1(mc_params.sweeps / 2, 0.0),
          valid_data_of_magnetization_abs2(mc_params.sweeps / 2, 0.0),
          valid_data_of_magnetization_abs3(mc_params.sweeps / 2, 0.0);
      for (size_t j = 0; j < valid_data_of_magnetization_abs1.size(); j++) {
        valid_data_of_magnetization_abs1[j] = abs(zigzag_af_magnetization_[i][j + sweeps / 2]);
        valid_data_of_magnetization_abs2[j] = abs(zigzag_af_magnetization2_[i][j + sweeps / 2]);
        valid_data_of_magnetization_abs3[j] = abs(zigzag_af_magnetization3_[i][j + sweeps / 2]);
      }
      res_.push_back(Variance((valid_data_of_magnetization_abs1)) / phys_params.N * phys_params.beta);
      res_.push_back(Variance((valid_data_of_magnetization_abs2)) / phys_params.N * phys_params.beta);
      res_.push_back(Variance((valid_data_of_magnetization_abs3)) / phys_params.N * phys_params.beta);
    }
  }
  if (DimDof > 1 && interaction_isotropy_) {
    // stiffness
    res_.push_back(Mean(std::vector(stiffness_.begin() + mc_params.sweeps / 2, stiffness_.end())));
  }

  if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {
    // zig-zag anti-ferromagnetic magnetization
    std::array<std::vector<double>, DimDof> valid_data_of_magnetization_squareQ1;
    std::array<std::vector<double>, DimDof> valid_data_of_magnetization_squareQ2;
    std::array<std::vector<double>, DimDof> valid_data_of_magnetization_squareQ3;
    // Q1
    double M2Q1 = 0.0;
    for (size_t dof = 0; dof < DimDof; dof++) {
      valid_data_of_magnetization_squareQ1[dof] = std::vector<double>(sweeps / 2, 0.0);
      for (size_t i = 0; i < valid_data_of_magnetization_squareQ1[dof].size(); i++) {
        valid_data_of_magnetization_squareQ1[dof][i] =
            zigzag_af_magnetization_[dof][i + sweeps / 2] * zigzag_af_magnetization_[dof][i + sweeps / 2];
      }
      M2Q1 += Mean(valid_data_of_magnetization_squareQ1[dof]);
      res_.push_back(Mean(valid_data_of_magnetization_squareQ1[dof]) / phys_params.N / phys_params.N);
    }
    // Q2
    double M2Q2 = 0.0;
    for (size_t dof = 0; dof < DimDof; dof++) {
      valid_data_of_magnetization_squareQ2[dof] = std::vector<double>(sweeps / 2, 0.0);
      for (size_t i = 0; i < valid_data_of_magnetization_squareQ2[dof].size(); i++) {
        valid_data_of_magnetization_squareQ2[dof][i] =
            zigzag_af_magnetization2_[dof][i + sweeps / 2] * zigzag_af_magnetization2_[dof][i + sweeps / 2];
      }
      M2Q2 += Mean(valid_data_of_magnetization_squareQ2[dof]);
      res_.push_back(Mean(valid_data_of_magnetization_squareQ2[dof]) / phys_params.N / phys_params.N);
    }
    // Q3
    double M2Q3 = 0.0;
    for (size_t dof = 0; dof < DimDof; dof++) {
      valid_data_of_magnetization_squareQ3[dof] = std::vector<double>(sweeps / 2, 0.0);
      for (size_t i = 0; i < valid_data_of_magnetization_squareQ3[dof].size(); i++) {
        valid_data_of_magnetization_squareQ3[dof][i] =
            zigzag_af_magnetization3_[dof][i + sweeps / 2] * zigzag_af_magnetization3_[dof][i + sweeps / 2];
      }
      M2Q3 += Mean(valid_data_of_magnetization_squareQ3[dof]);
      res_.push_back(Mean(valid_data_of_magnetization_squareQ3[dof]) / phys_params.N / phys_params.N);
    }


    //Binder ratio
    std::vector<double> valid_data_of_M2_Q1(sweeps / 2, 0.0), valid_data_of_M2_Q2(sweeps / 2, 0.0),
        valid_data_of_M2_Q3(sweeps / 2, 0.0);
    std::vector<double> valid_data_of_magnetization_quarticQ1(sweeps / 2, 0.0),
        valid_data_of_magnetization_quarticQ2(sweeps / 2, 0.0),
        valid_data_of_magnetization_quarticQ3(sweeps / 2, 0.0);

    double binder_ratio1, binder_ratio2, binder_ratio3;
    //Q1
    for (size_t i = 0; i < valid_data_of_magnetization_quarticQ1.size(); i++) {
      double sum(0.0);
      for (size_t dof = 0; dof < DimDof; dof++) {
        sum += valid_data_of_magnetization_squareQ1[dof][i];
      }
      valid_data_of_M2_Q1[i] = sum;
      valid_data_of_magnetization_quarticQ1[i] = sum * sum;
    }
    double M4Q1 = Mean(valid_data_of_magnetization_quarticQ1);
    binder_ratio1 = M4Q1 / M2Q1 / M2Q1;
    //Q2
    for (size_t i = 0; i < valid_data_of_magnetization_quarticQ2.size(); i++) {
      double sum(0.0);
      for (size_t dof = 0; dof < DimDof; dof++) {
        sum += valid_data_of_magnetization_squareQ2[dof][i];
      }
      valid_data_of_M2_Q2[i] = sum;
      valid_data_of_magnetization_quarticQ2[i] = sum * sum;
    }
    double M4Q2 = Mean(valid_data_of_magnetization_quarticQ2);
    binder_ratio2 = M4Q2 / M2Q2 / M2Q2;
    //Q3
    for (size_t i = 0; i < valid_data_of_magnetization_quarticQ3.size(); i++) {
      double sum(0.0);
      for (size_t dof = 0; dof < DimDof; dof++) {
        sum += valid_data_of_magnetization_squareQ3[dof][i];
      }
      valid_data_of_M2_Q3[i] = sum;
      valid_data_of_magnetization_quarticQ3[i] = sum * sum;
    }
    double M4Q3 = Mean(valid_data_of_magnetization_quarticQ3);
    binder_ratio3 = M4Q3 / M2Q3 / M2Q3;

    std::vector<double> valid_data_of_M2(sweeps / 2, 0.0); //m^2 = m_1^2 + m_2^2 + m_3^2
    std::vector<double> valid_data_of_magnetization_quartic(sweeps / 2, 0.0); // m^4
    for (size_t i = 0; i < valid_data_of_M2.size(); i++) {
      valid_data_of_M2[i] = valid_data_of_M2_Q1[i] + valid_data_of_M2_Q2[i] + valid_data_of_M2_Q3[i];
      valid_data_of_magnetization_quartic[i] = valid_data_of_M2[i] * valid_data_of_M2[i];
    }
    double M4 = Mean(valid_data_of_magnetization_quartic);
    double M2 = Mean(valid_data_of_M2);

    res_.push_back(binder_ratio1);
    res_.push_back(binder_ratio2);
    res_.push_back(binder_ratio3);
    res_.push_back(M4 / M2 / M2);




    //C3 symmetry breaking
    std::vector<double> data_of_c3_order_parameter(sweeps / 2, 0.0);
    std::vector<double> data_of_c3_order_parameter_square(sweeps / 2, 0.0);
    std::vector<double> data_of_c3_order_parameter_abs(sweeps / 2, 0.0);
    for (size_t i = 0; i < data_of_c3_order_parameter.size(); i++) {
      data_of_c3_order_parameter[i] = valid_data_of_M2_Q1[i] + valid_data_of_M2_Q2[i] - 2 * valid_data_of_M2_Q3[i];
      data_of_c3_order_parameter_square[i] = data_of_c3_order_parameter[i] * data_of_c3_order_parameter[i];
      data_of_c3_order_parameter_abs[i] = abs(data_of_c3_order_parameter[i]);
    }

    size_t N2 = phys_params.N * phys_params.N;
    size_t N4 = N2 * N2;
    res_.push_back(Mean(data_of_c3_order_parameter_square) / N4); // <b^2>, assume b = 0.0;
    res_.push_back(
        Variance(data_of_c3_order_parameter_abs) / N2 * phys_params.beta); // susceptibility of C3 order parameter b

  }

  //spin correlation (L/2), (L/4)
//  if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {
//    for (size_t i = 0; i < DimDof; i++) {
//      auto half_data_of_correlation2 = std::vector(correlation_lover2_[i].begin() + mc_params.sweeps/2, correlation_lover2_[i].end());
//      res_.push_back(Mean(half_data_of_correlation2));
//      auto half_data_of_correlation4 = std::vector(correlation_lover4_[i].begin() + mc_params.sweeps/2, correlation_lover4_[i].end());
//      res_.push_back(Mean(half_data_of_correlation4));
//    }
//  }


  if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {
    size_t sweeps = mc_params.sweeps;
    double msquare, mquartic;
    if (has_gs_data) {
      std::vector<double> m2_data(sweeps / 2);
      std::vector<double> m4_data(sweeps / 2);
      for (size_t i = 0; i < sweeps / 2; i++) {
        auto &real = complex_xy_order_parameter_real_[i + sweeps / 2];
        auto &imag = complex_xy_order_parameter_imag_[i + sweeps / 2];
        m2_data[i] = real * real + imag * imag;
        m4_data[i] = m2_data[i] * m2_data[i];
      }
      msquare = Mean(m2_data);
      mquartic = Mean(m4_data);
      res_.push_back(msquare);

      double binder_ratio_m = mquartic / msquare / msquare;
      res_.push_back(binder_ratio_m);
    }
  }
  DumpData("summary" + mc_params.filename_postfix, res_);

//  auto half_data_of_correlation2 =
//      std::vector(correlation_lover2_.begin() + mc_params.sweeps / 2, correlation_lover2_.end());
//  res_.push_back(Mean(half_data_of_correlation2));
//  auto half_data_of_correlation4 =
//      std::vector(correlation_lover4_.begin() + mc_params.sweeps / 2, correlation_lover4_.end());
//  res_.push_back(Mean(half_data_of_correlation4));
//  res_.push_back(res_[0] / res_[1]);
//  DumpData("correlation_summary" + mc_params.filename_postfix, res_);
  statistic_timer.PrintElapsed();
}

#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_WOLFMCEXECUTOR_H
