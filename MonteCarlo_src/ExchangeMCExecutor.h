/**
 * Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
 * Created Data: 2022/10/13.
 *
 * Description: Exchange Monte Carlo with cluster update.
 */


#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_ExchangeMCEXECUTOR_H
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_ExchangeMCEXECUTOR_H

#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <omp.h>
#include "gqten/framework/bases/executor.h"
#include "gqten/utility/timer.h"
#include "swcluster.h"
#include "lattice_link.h"
#include "lattice_config.h"
#include "common.h"

struct ExchangeMCParams {
  size_t thread_num;      // parallel thread number, it's also the number of temperatures calculated at the same time
  size_t adjust_temperature_times;
  size_t adjust_temperature_samples;
  size_t sweeps;          // how many samples/sweeps
  size_t sample_interval; // how many times of wolf update between two samples. 0 means metropolis algorithm.
  size_t exchange_interval; // how many times of measure (MC sweep) between two replica exchanges.
  size_t warmup_sample_num; // how many samples are used to thermalization
  size_t print_interval;    // how many sample interval between which the program print some informations
  std::string filename_postfix; // used to dump data
};

template<size_t DimDof, size_t NumOfCouplingType>
struct PhysParams {
  PhysParams(const std::string &geometry,
             size_t lx,
             size_t ly,
             size_t n,
             double beta_min,
             double beta_max,
             const std::array<CouplingStructure<DimDof>, NumOfCouplingType> &coupling_structures)
      : geometry(geometry), Lx(lx), Ly(ly), N(n), beta_min(beta_min),
        beta_max(beta_max), coupling_structures(coupling_structures) {
    zigzag_gs_oritation_1.Zero();
    zigzag_gs_oritation_2.Zero();
    zigzag_gs_oritation_3.Zero();
  }

  PhysParams(const std::string &geometry,
             size_t lx,
             size_t ly,
             size_t n,
             double beta_min,
             double beta_max,
             const std::array<CouplingStructure<DimDof>, NumOfCouplingType> &coupling_structures,
             const LocalDOF<DimDof> &zigzag_gs_oritation_1,
             const LocalDOF<DimDof> &zigzag_gs_oritation_2,
             const LocalDOF<DimDof> &zigzag_gs_oritation_3)
      : geometry(geometry),
        Lx(lx),
        Ly(ly),
        N(n),
        beta_min(beta_min),
        beta_max(beta_max),
        coupling_structures(coupling_structures),
        zigzag_gs_oritation_1(zigzag_gs_oritation_1),
        zigzag_gs_oritation_2(zigzag_gs_oritation_2),
        zigzag_gs_oritation_3(zigzag_gs_oritation_3) {}

  PhysParams(const std::string &geometry,
             size_t lx,
             size_t ly,
             double beta_min,
             double beta_max,
             const std::array<CouplingStructure<DimDof>, NumOfCouplingType> &coupling_structures)
      : geometry(geometry), Lx(lx), Ly(ly), beta_min(beta_min),
        beta_max(beta_max), coupling_structures(coupling_structures) {
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
             double beta_min,
             double beta_max,
             const std::array<CouplingStructure<DimDof>, NumOfCouplingType> &coupling_structures,
             const LocalDOF<DimDof> &zigzag_gs_oritation_1,
             const LocalDOF<DimDof> &zigzag_gs_oritation_2,
             const LocalDOF<DimDof> &zigzag_gs_oritation_3)
      : geometry(geometry),
        Lx(lx),
        Ly(ly),
        beta_min(beta_min),
        beta_max(beta_max),
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
  double beta_min; // 1/T_max
  double beta_max; // 1/T_min
  std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures;
  LocalDOF<DimDof> zigzag_gs_oritation_1;
  LocalDOF<DimDof> zigzag_gs_oritation_2;
  LocalDOF<DimDof> zigzag_gs_oritation_3;
};

template<size_t DimDof>
struct StatisticResults {
  StatisticResults(size_t N) : beta(N),
                               energy(N),
                               specific_heat(N),
                               stiffness(N),
                               binder_ratio_zigzag_direction1(N),
                               binder_ratio_zigzag_direction2(N),
                               binder_ratio_zigzag_direction3(N),
                               binder_ratio_zigzag_af(N),
                               zigzag_magnetization_susceptibility_b(N),
                               c3_order_parameter_square(N),
                               c3_order_parameter_susceptibility(N),
                               binder_ratio_c3(N),
                               complex_order_parameter(N),
                               binder_ratio_complex_order_parameter(N) {
    for (size_t dof = 0; dof < DimDof; dof++) {
      zigzag_af1_magnetization_square[dof] = std::vector<double>(N);
      zigzag_af2_magnetization_square[dof] = std::vector<double>(N);
      zigzag_af3_magnetization_square[dof] = std::vector<double>(N);
      zigzag_magnetization_susceptibility[dof] = std::vector<double>(N);
    }
  }

  void DumpData(const std::string &filename) {
    std::ofstream ofs(filename, std::ofstream::binary);
    ofs.write((const char *) beta.data(), beta.size() * sizeof(double));
    ofs << std::endl;
    ofs.write((const char *) energy.data(), energy.size() * sizeof(double));
    ofs << std::endl;

    ofs.write((const char *) specific_heat.data(), specific_heat.size() * sizeof(double));
    ofs << std::endl;

    ofs.write((const char *) stiffness.data(), stiffness.size() * sizeof(double));
    ofs << std::endl;

    for (size_t dof = 0; dof < DimDof; dof++) {
      ofs.write((const char *) zigzag_af1_magnetization_square[dof].data(),
                zigzag_af1_magnetization_square[dof].size() * sizeof(double));
      ofs << std::endl;
    }
    for (size_t dof = 0; dof < DimDof; dof++) {
      ofs.write((const char *) zigzag_af2_magnetization_square[dof].data(),
                zigzag_af2_magnetization_square[dof].size() * sizeof(double));
      ofs << std::endl;
    }
    for (size_t dof = 0; dof < DimDof; dof++) {
      ofs.write((const char *) zigzag_af3_magnetization_square[dof].data(),
                zigzag_af3_magnetization_square[dof].size() * sizeof(double));
      ofs << std::endl;
    }
    for (size_t dof = 0; dof < DimDof; dof++) {
      ofs.write((const char *) zigzag_magnetization_susceptibility[dof].data(),
                zigzag_magnetization_susceptibility[dof].size() * sizeof(double));
      ofs << std::endl;
    }
    ofs.write((const char *) binder_ratio_zigzag_direction1.data(),
              binder_ratio_zigzag_direction1.size() * sizeof(double));
    ofs << std::endl;
    ofs.write((const char *) binder_ratio_zigzag_direction2.data(),
              binder_ratio_zigzag_direction2.size() * sizeof(double));
    ofs << std::endl;
    ofs.write((const char *) binder_ratio_zigzag_direction3.data(),
              binder_ratio_zigzag_direction3.size() * sizeof(double));
    ofs << std::endl;
    ofs.write((const char *) binder_ratio_zigzag_af.data(), binder_ratio_zigzag_af.size() * sizeof(double));
    ofs << std::endl;
    ofs.write((const char *) zigzag_magnetization_susceptibility_b.data(),
              zigzag_magnetization_susceptibility_b.size() * sizeof(double));
    ofs << std::endl;

    ofs.write((const char *) c3_order_parameter_square.data(), c3_order_parameter_square.size() * sizeof(double));
    ofs << std::endl;
    ofs.write((const char *) c3_order_parameter_susceptibility.data(),
              c3_order_parameter_susceptibility.size() * sizeof(double));
    ofs << std::endl;
    ofs.write((const char *) complex_order_parameter.data(), complex_order_parameter.size() * sizeof(double));
    ofs << std::endl;
    ofs.write((const char *) binder_ratio_complex_order_parameter.data(),
              binder_ratio_complex_order_parameter.size() * sizeof(double));
    ofs << std::endl;

    ofs.write((const char *) binder_ratio_c3.data(),
              binder_ratio_c3.size() * sizeof(double));
    ofs << std::endl;
    ofs.close();
  }

  std::vector<double> beta;
  // the index of vectors corresponds to the beta
  std::vector<double> energy;         // energy per site
  std::vector<double> specific_heat;
  std::vector<double> stiffness;
  std::array<std::vector<double>, DimDof> zigzag_af1_magnetization_square;
  // <M(Q_1)^2>, M(Q_1) is per-site zig-zag antiferromagnetic magnetization in the direction corresponding to momentum Q_1.
  std::array<std::vector<double>, DimDof> zigzag_af2_magnetization_square;
  std::array<std::vector<double>, DimDof> zigzag_af3_magnetization_square;
  std::array<std::vector<double>, DimDof> zigzag_magnetization_susceptibility;
  std::vector<double> binder_ratio_zigzag_direction1;   // <M(Q_1)^4>/<M(Q_1)^2>^2
  std::vector<double> binder_ratio_zigzag_direction2;
  std::vector<double> binder_ratio_zigzag_direction3;
  std::vector<double> binder_ratio_zigzag_af;     // <(M(Q_1)^2+M(Q_2)^2+M(Q_3)^2)^2>/<M(Q_1)^2+M(Q_2)^2+M(Q_3)^2>^2
  std::vector<double> zigzag_magnetization_susceptibility_b;

  std::vector<double> c3_order_parameter_square;  // <C3> = <(M(Q_1)^2+M(Q_2)^2-2*M(Q_3)^2)^2>
  std::vector<double> c3_order_parameter_susceptibility;
  std::vector<double> binder_ratio_c3; // <C3^2>/<|C3|>^2

  std::vector<double> complex_order_parameter;
  std::vector<double> binder_ratio_complex_order_parameter;
};

template<size_t DimDof, size_t NumOfCouplingType>
class ExchangeMCExecutor : public gqten::Executor {
  using LocalDOFT = LocalDOF<DimDof>;
 public:
  ExchangeMCExecutor(const ExchangeMCParams &,
                     const PhysParams<DimDof, NumOfCouplingType> &);
  void Execute() override;

  ~ExchangeMCExecutor() {
    for (size_t i = 0; i < mc_params.thread_num; i++) {
      delete configs_[i];
    }
    delete plattice_link_;
  }

  ExchangeMCParams mc_params;
  PhysParams<DimDof, NumOfCouplingType> phys_params;

 private:

  void WolfSweep_(const size_t replica_id);
  void WolfRandomSweep_(const size_t replica_id, size_t sample_interval);
  void WolfGrowCluster_(const size_t replica_id, const LocalDOFT &, const size_t);
  bool WolfFlipCluster_(const size_t replica_id, const LocalDOFT &);

  void MetropolisSweep_(const size_t replica_id);
  bool MetropolisSingleStep_(const size_t replica_id, const size_t site);

  void ExchangeSweep_();
  bool ExchangeReplica_(size_t, size_t);
  bool AdjustTemperatureSpace();

  void MeasureEnergy_(size_t replica_id);
  void Measure_(size_t replica_id);

  void DumpSamplesData_();
  void StatisticAndDumpResults_();

  size_t geometry_id_; // 0: square, 1: honeycomb, 2: SVW
  bool interaction_isotropy_; // not consider SIA
  bool isotropy_;
  LatticeLink<DimDof, NumOfCouplingType> *plattice_link_;
  LocalDOF<DimDof> ground_state_oritation1_; //corresponding to the zig-zag order 1
  LocalDOF<DimDof> ground_state_oritation2_;
  LocalDOF<DimDof> ground_state_oritation3_;
  bool has_gs_data;

  std::vector<double> beta_set_; //from small beta to large beta
  std::vector<LatticeConfig<DimDof> *> configs_;
  std::vector<SWCluster> clusters_;
  std::vector<double> energys_now_;
  std::vector<size_t> cluster_flip_times_;
  std::vector<size_t> replica_exchage_times_;
  std::vector<double> exchange_probability_;

  /// result data
  std::vector<std::vector<double>> energy_set_; // index of outer vector corresponds beta
  std::vector<std::array<std::vector<double>, DimDof>> sum_spin_set_;
  std::vector<std::array<std::vector<double>, DimDof>> zigzag_af_magnetization1_set_;
  std::vector<std::array<std::vector<double>, DimDof>> zigzag_af_magnetization2_set_;
  std::vector<std::array<std::vector<double>, DimDof>> zigzag_af_magnetization3_set_;
  std::vector<double> correlation_lover2_;
  std::vector<double> correlation_lover4_;
  std::vector<std::vector<double>> stiffness_set_;         // if U(1) symmetric model
  std::vector<std::vector<double>> complex_xy_order_parameter_real_set_;
  std::vector<std::vector<double>> complex_xy_order_parameter_imag_set_;

  StatisticResults<DimDof> results_;

  // persite energy, specific heat, magnetic susceptibility, stiffness,

  std::vector<std::default_random_engine> random_number_generators;
};

template<size_t DimDof, size_t NumOfCouplingType>
ExchangeMCExecutor<DimDof, NumOfCouplingType>::ExchangeMCExecutor(const ExchangeMCParams &mc_params,
                                                                  const PhysParams<DimDof,
                                                                                   NumOfCouplingType> &phys_params) :
    gqten::Executor(),
    mc_params(mc_params),
    phys_params(phys_params),
    beta_set_(mc_params.thread_num),
    configs_(mc_params.thread_num, nullptr),
    clusters_(mc_params.thread_num),
    energys_now_(mc_params.thread_num),
    cluster_flip_times_(mc_params.thread_num, 0),
    replica_exchage_times_(mc_params.thread_num, 0),
    exchange_probability_(mc_params.thread_num, 0.0),
    energy_set_(mc_params.thread_num),
    sum_spin_set_(mc_params.thread_num),
    zigzag_af_magnetization1_set_(mc_params.thread_num),
    zigzag_af_magnetization2_set_(mc_params.thread_num),
    zigzag_af_magnetization3_set_(mc_params.thread_num),
    stiffness_set_(mc_params.thread_num),
    complex_xy_order_parameter_real_set_(mc_params.thread_num),
    complex_xy_order_parameter_imag_set_(mc_params.thread_num),
    results_(mc_params.thread_num),
    random_number_generators(mc_params.thread_num) {
  if (phys_params.geometry == "Square") {
    plattice_link_ = new SquareTorusLatticeLink(phys_params.Lx, phys_params.Ly, phys_params.coupling_structures);
    geometry_id_ = 0;
  } else if (phys_params.geometry == "Honeycomb") {
    plattice_link_ = new HoneyCombTorusLatticeLink(phys_params.Lx, phys_params.Ly, phys_params.coupling_structures);
    geometry_id_ = 1;
  } else if (phys_params.geometry == "SVW") {
    plattice_link_ = new SVWLatticeLink(phys_params.N, phys_params.coupling_structures);
    geometry_id_ = 2;
  } else if (phys_params.geometry == "Honeycomb3") {
    plattice_link_ = new HoneyComb3TorusLatticeLink(phys_params.Lx, phys_params.Ly, phys_params.coupling_structures);
    geometry_id_ = 3;
  } else {
    std::cout << "do not support now. " << std::endl;
    exit(0);
  }

  interaction_isotropy_ = true;
  for (size_t i = 1; i < NumOfCouplingType; i++) {
    interaction_isotropy_ = interaction_isotropy_ && phys_params.coupling_structures[i].IsIsometry();
  }

  isotropy_ = interaction_isotropy_ && phys_params.coupling_structures[0].IsIsometry();

  ground_state_oritation1_ = phys_params.zigzag_gs_oritation_1;
  ground_state_oritation2_ = phys_params.zigzag_gs_oritation_2;
  ground_state_oritation3_ = phys_params.zigzag_gs_oritation_3;
  if (ground_state_oritation1_.IsZero() && ground_state_oritation2_.IsZero() && ground_state_oritation3_.IsZero()) {
    has_gs_data = false;
  } else {
    has_gs_data = true;
  }

  size_t sweeps = mc_params.sweeps;
  double d_beta = (phys_params.beta_max - phys_params.beta_min) / (mc_params.thread_num - 1);
  for (size_t replica_id = 0; replica_id < mc_params.thread_num; replica_id++) {
    beta_set_[replica_id] = phys_params.beta_min + d_beta * replica_id;

    configs_[replica_id] = new LatticeConfig<DimDof>(phys_params.N);
    configs_[replica_id]->Random(random_number_generators[replica_id]);

    energy_set_[replica_id].reserve(sweeps);
    if (DimDof >= 2 && interaction_isotropy_) {
      stiffness_set_[replica_id].reserve(sweeps);
    }

    for (size_t i = 0; i < DimDof; i++) {
      sum_spin_set_[replica_id][i].reserve(sweeps);
      if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {
        zigzag_af_magnetization1_set_[replica_id][i].reserve(sweeps);
        zigzag_af_magnetization2_set_[replica_id][i].reserve(sweeps);
        zigzag_af_magnetization3_set_[replica_id][i].reserve(sweeps);
        if (has_gs_data) {
          complex_xy_order_parameter_real_set_[replica_id].reserve(sweeps);
          complex_xy_order_parameter_imag_set_[replica_id].reserve(sweeps);
        }
      }
    }
  }

  for (size_t replica_id = 0; replica_id < mc_params.thread_num; replica_id++) {
    random_number_generators[replica_id].seed((int) clock() + replica_id * 97);
  }

  SetStatus(gqten::INITED);
}

template<size_t DimDof, size_t NumOfCouplingType>
void ExchangeMCExecutor<DimDof, NumOfCouplingType>::Execute() {
  SetStatus(gqten::EXEING);
  gqten::Timer wolf_mc_execute_timer("wolf_mc_execute");
  for (size_t adjust_time = 0; adjust_time < mc_params.adjust_temperature_times; adjust_time++) {
    std::cout << "adjust " << adjust_time << std::endl;
    for (size_t sweep = 0; sweep < mc_params.adjust_temperature_samples; sweep++) {
#pragma omp parallel for default(none) \
                shared(mc_params, phys_params)\
                num_threads(mc_params.thread_num)\
                schedule(static)
      for (size_t replica_id = 0; replica_id < mc_params.thread_num; replica_id++) {
        MetropolisSweep_(replica_id);

        if (mc_params.sample_interval != phys_params.N) {
          WolfRandomSweep_(replica_id, mc_params.sample_interval);
        } else {
          WolfSweep_(replica_id);
        }

      }
      if (sweep % mc_params.exchange_interval == 0) {
#pragma omp parallel for default(none) \
                shared(mc_params, phys_params)\
                num_threads(mc_params.thread_num)\
                schedule(static)
        for (size_t replica_id = 0; replica_id < mc_params.thread_num; replica_id++) {
          MeasureEnergy_(replica_id);
        }
        ExchangeSweep_();
      }
      if (sweep % mc_params.print_interval == 0) {
        double execute_time = wolf_mc_execute_timer.Elapsed();
        std::cout << "[ sweep = " << sweep << " ]"
                  << "time = " << execute_time << "\t"
                  << "cluster size = " << clusters_.back().size() // for the lowest temperature
                  << std::endl;
      }
    }
    if (AdjustTemperatureSpace()) break;
  }

  for (size_t sweep = 0; sweep < mc_params.sweeps; sweep++) {
#pragma omp parallel for default(none) \
                shared(mc_params, phys_params)\
                num_threads(mc_params.thread_num)\
                schedule(static)
    for (size_t replica_id = 0; replica_id < mc_params.thread_num; replica_id++) {
      MetropolisSweep_(replica_id);

      if (mc_params.sample_interval != phys_params.N) {
        WolfRandomSweep_(replica_id, mc_params.sample_interval);
      } else {
        WolfSweep_(replica_id);
      }

      Measure_(replica_id);
    }
    if (sweep % mc_params.exchange_interval == 0) {
      ExchangeSweep_();
    }
    if (sweep % mc_params.print_interval == 0) {
      double execute_time = wolf_mc_execute_timer.Elapsed();
      std::cout << "[ sweep = " << sweep << " ]"
                << "time = " << execute_time << "\t"
                << "cluster size = " << clusters_.back().size() // for the lowest temperature
                << std::endl;
    }
  }
  DumpSamplesData_();
  StatisticAndDumpResults_();
  for (size_t replica_id = 0; replica_id < mc_params.thread_num; replica_id++) {
    std::cout << "beta = " << beta_set_[replica_id]
              << ", Energy = " << results_.energy[replica_id]
              << ", Flip probability: "
              << double(cluster_flip_times_[replica_id]) / double(mc_params.sweeps * mc_params.sample_interval)
              << ", Exchange probability: "
              << double(replica_exchage_times_[replica_id]) / double(mc_params.sweeps)
                  * double(mc_params.exchange_interval)
              << std::endl;
  }
  SetStatus(gqten::FINISH);
}

template<size_t DimDof, size_t NumOfCouplingType>
void ExchangeMCExecutor<DimDof, NumOfCouplingType>::WolfSweep_(const size_t replica_id) {
  LocalDOF<DimDof> axis;
  for (size_t starting_site = 0; starting_site < phys_params.N; starting_site++) {
    axis.Random(random_number_generators[replica_id]);
    WolfGrowCluster_(replica_id, axis, starting_site);
    if (WolfFlipCluster_(replica_id, axis)) cluster_flip_times_[replica_id]++;
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
void ExchangeMCExecutor<DimDof, NumOfCouplingType>::WolfRandomSweep_(const size_t replica_id,
                                                                     size_t sample_interval) {
  LocalDOF<DimDof> axis;
  std::uniform_int_distribution<int> u_int(0, phys_params.N - 1);
  for (size_t i = 0; i < sample_interval; i++) {
    axis.Random(random_number_generators[replica_id]);
    WolfGrowCluster_(replica_id, axis, u_int(random_number_generators[replica_id]));
    if (WolfFlipCluster_(replica_id, axis)) cluster_flip_times_[replica_id]++;
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
void ExchangeMCExecutor<DimDof, NumOfCouplingType>::WolfGrowCluster_(const size_t replica_id,
                                                                     const LocalDOF<DimDof> &axis,
                                                                     const size_t starting_site) {
  SWCluster cluster_local = {starting_site};
  cluster_local.reserve(phys_params.N);
  std::vector<bool> in_cluster(phys_params.N, false);
  in_cluster[starting_site] = true;
// config_.PrintSign(axis, phys_params.Lx);
// Note here the effective single ion anisotropy is not important.
  double K_eff[NumOfCouplingType - 1];
  for (size_t inter_type = 1; inter_type < NumOfCouplingType; inter_type++) {
    double eff_coupling = axis * (phys_params.coupling_structures[inter_type] * axis);
    K_eff[inter_type - 1] = eff_coupling * beta_set_[replica_id];
  }

  for (size_t j = 0; j < cluster_local.size(); j++) {
    for (size_t inter_type = 1; inter_type < NumOfCouplingType; inter_type++) {
      const std::vector<size_t> &site2_set = plattice_link_->GetSitesLinkedTo(cluster_local[j], inter_type);
      for (auto site2: site2_set) {
        if (!in_cluster[site2]
            && configs_[replica_id]->Active(cluster_local[j],
                                            site2,
                                            K_eff[inter_type - 1],
                                            axis,
                                            u_double(random_number_generators[replica_id]))) {
          in_cluster[site2] = true;
          cluster_local.push_back(site2);
        }
      }
    }
  }
  clusters_[replica_id] = std::move(cluster_local);
}

/**
 *
 * @param axis
 * @return flipped or not
 */
template<size_t DimDof, size_t NumOfCouplingType>
bool ExchangeMCExecutor<DimDof, NumOfCouplingType>::WolfFlipCluster_(const size_t replica_id,
                                                                     const LocalDOF<DimDof> &axis) {

  if (isotropy_) {
    configs_[replica_id]->Flip(axis, clusters_[replica_id]);
    return true;
  }
  double delta_e;
  if (interaction_isotropy_) {
    delta_e =
        configs_[replica_id]->EnergyDifferenceOfSIAFlipCluster(axis,
                                                               clusters_[replica_id],
                                                               plattice_link_->coupling_structures[0]);
  } else {
    delta_e = configs_[replica_id]->EnergyDifferenceFlipCluster(axis, clusters_[replica_id], *plattice_link_);
  }
  if (delta_e <= 0.0) {
    configs_[replica_id]->Flip(axis, clusters_[replica_id]);
    return true;
  } else {
    std::uniform_real_distribution<double> u(0, 1);
    if (u(random_number_generators[replica_id]) <= std::exp(-beta_set_[replica_id] * delta_e)) {
      configs_[replica_id]->Flip(axis, clusters_[replica_id]);
      return true;
    } else {
      return false;
    }
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
void ExchangeMCExecutor<DimDof, NumOfCouplingType>::MetropolisSweep_(const size_t replica_id) {
  for (size_t i = 0; i < phys_params.N; i++) {
    MetropolisSingleStep_(replica_id, i);
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
bool ExchangeMCExecutor<DimDof, NumOfCouplingType>::MetropolisSingleStep_(const size_t replica_id,
                                                                          const size_t site) {
  LocalDOF<DimDof> flipto;
  flipto.Random(random_number_generators[replica_id]);
  double delta_e =
      configs_[replica_id]->template EnergyDifferenceFlipSite<NumOfCouplingType>(site, flipto, *plattice_link_);
  if (delta_e <= 0.0) {
    configs_[replica_id]->SetSiteSpin(site, flipto);
    return true;
  } else {
    std::uniform_real_distribution<double> u(0, 1);
    if (u(random_number_generators[replica_id]) <= std::exp(-beta_set_[replica_id] * delta_e)) {
      configs_[replica_id]->SetSiteSpin(site, flipto);
      return true;
    } else {
      return false;
    }
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
void ExchangeMCExecutor<DimDof, NumOfCouplingType>::ExchangeSweep_() {
#pragma omp parallel for default(none) \
                shared(replica_exchage_times_)\
                num_threads(mc_params.thread_num/2)\
                schedule(static)
  for (size_t replica_id = 0; replica_id < mc_params.thread_num - 1; replica_id = replica_id + 2) {
    if (ExchangeReplica_(replica_id, replica_id + 1)) replica_exchage_times_[replica_id]++;
  }
#pragma omp parallel for default(none) \
                shared(replica_exchage_times_)\
                num_threads(mc_params.thread_num/2)\
                schedule(static)
  for (size_t replica_id = 1; replica_id < mc_params.thread_num - 1; replica_id = replica_id + 2) {
    if (ExchangeReplica_(replica_id, replica_id + 1)) replica_exchage_times_[replica_id]++;
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
bool ExchangeMCExecutor<DimDof, NumOfCouplingType>::ExchangeReplica_(size_t replica_ida, size_t replica_idb) {
  double delta =
      (beta_set_[replica_idb] - beta_set_[replica_ida]) * (energys_now_[replica_ida] - energys_now_[replica_idb]);
  if (delta < 0) {
    std::swap(configs_[replica_ida], configs_[replica_idb]);
    return true;
  } else {
    std::uniform_real_distribution<double> u(0, 1);
    if (u(random_number_generators[replica_ida]) <= std::exp(-delta)) {
      std::swap(configs_[replica_ida], configs_[replica_idb]);
      return true;
    } else {
      return false;
    }
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
bool ExchangeMCExecutor<DimDof, NumOfCouplingType>::AdjustTemperatureSpace() {
  double alpha = 0.1; //tolerance
  double sum_p = 0.0;
  double p_max(0.0), p_min(1.0);
  std::cout << "Exchange probability: [";
  for (size_t i = 0; i < mc_params.thread_num - 1; i++) {
    exchange_probability_[i] = double(replica_exchage_times_[i]) / double(mc_params.adjust_temperature_samples)
        * double(mc_params.exchange_interval);
    std::cout << exchange_probability_[i] << ", ";
    sum_p += exchange_probability_[i];
    p_max = std::max(p_max, exchange_probability_[i]);
    p_min = std::min(p_min, exchange_probability_[i]);
    replica_exchage_times_[i] = 0;
  }
  std::cout << "0.0].\n";
  if (((p_max - p_min) / p_max) < alpha) {
    std::cout << "Probability almostly even." << std::endl;
    return true;
  }

  double c = sum_p / (mc_params.thread_num - 1);
  std::cout << "Update beta as: [" << beta_set_[0];
  for (size_t i = 1; i < mc_params.thread_num; i++) {
    beta_set_[i] = beta_set_[i - 1] + (beta_set_[i] - beta_set_[i - 1]) * exchange_probability_[i - 1] / c;
    std::cout << ", " << beta_set_[i];
  }
  std::cout << "]." << std::endl;
  return false;
}

template<size_t DimDof, size_t NumOfCouplingType>
void ExchangeMCExecutor<DimDof, NumOfCouplingType>::MeasureEnergy_(size_t replica_id) {
  energys_now_[replica_id] = configs_[replica_id]->Energy(*plattice_link_);
}

template<size_t DimDof, size_t NumOfCouplingType>
void ExchangeMCExecutor<DimDof, NumOfCouplingType>::Measure_(size_t replica_id) {
  double total_energy = configs_[replica_id]->Energy(*plattice_link_);
  energys_now_[replica_id] = total_energy;
  energy_set_[replica_id].push_back(total_energy);

  if (DimDof >= 2) {
    if (interaction_isotropy_) {
      stiffness_set_[replica_id].push_back(configs_[replica_id]->template StiffnessIsotropy<NumOfCouplingType>(*plattice_link_,
                                                                                                               beta_set_[replica_id]));
    }
  }

  if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {
    auto m1 = configs_[replica_id]->HoneycombZigzagAntiferromagneticMagnetization(phys_params.Lx);
    auto m2 = configs_[replica_id]->HoneycombZigzagAntiferromagneticMagnetization2(phys_params.Lx);
    auto m3 = configs_[replica_id]->HoneycombZigzagAntiferromagneticMagnetization3(phys_params.Lx);

    for (size_t i = 0; i < DimDof; i++) {
      zigzag_af_magnetization1_set_[replica_id][i].push_back(m1[i]);
      zigzag_af_magnetization2_set_[replica_id][i].push_back(m2[i]);
      zigzag_af_magnetization3_set_[replica_id][i].push_back(m3[i]);
    }
    if (has_gs_data) {
      using complexT = std::complex<double>;
      complexT complex_order_parameter(0);
      complex_order_parameter += LocalDOF<DimDof>(m1) * ground_state_oritation1_;
      complex_order_parameter += LocalDOF<DimDof>(m2) * ground_state_oritation2_ * exp(complexT(0.0, 2 * M_PI / 3));
      complex_order_parameter += LocalDOF<DimDof>(m3) * ground_state_oritation3_ * exp(complexT(0.0, 4 * M_PI / 3));

      complex_xy_order_parameter_real_set_[replica_id].push_back(complex_order_parameter.real() / phys_params.N);
      complex_xy_order_parameter_imag_set_[replica_id].push_back(complex_order_parameter.imag() / phys_params.N);
    }
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
void ExchangeMCExecutor<DimDof, NumOfCouplingType>::DumpSamplesData_() {
  gqten::Timer dump_data_timer("dump_data");
#pragma omp parallel for default(none) \
                shared(mc_params, phys_params)\
                num_threads(mc_params.thread_num)\
                schedule(static)
  for (size_t replica_id = 0; replica_id < mc_params.thread_num; replica_id++) {
    DumpData("energy_beta" + std::to_string(beta_set_[replica_id]) + mc_params.filename_postfix,
             energy_set_[replica_id]);

    if (DimDof >= 2 && interaction_isotropy_) {
      DumpData("stiffness" + std::to_string(beta_set_[replica_id]) + mc_params.filename_postfix,
               stiffness_set_[replica_id]);
    }

    if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {
      for (size_t dof = 0; dof < DimDof; dof++) {
        DumpData(
            "magnetization1" + std::to_string(beta_set_[replica_id]) + std::to_string(dof)
                + mc_params.filename_postfix,
            zigzag_af_magnetization1_set_[replica_id][dof]);
        DumpData(
            "magnetization2" + std::to_string(beta_set_[replica_id]) + std::to_string(dof)
                + mc_params.filename_postfix,
            zigzag_af_magnetization2_set_[replica_id][dof]);
        DumpData(
            "magnetization3" + std::to_string(beta_set_[replica_id]) + std::to_string(dof)
                + mc_params.filename_postfix,
            zigzag_af_magnetization3_set_[replica_id][dof]);
        if (has_gs_data) {
          DumpData("xy_order_parameter_real" + std::to_string(beta_set_[replica_id]) + mc_params.filename_postfix,
                   complex_xy_order_parameter_real_set_[replica_id]);
          DumpData("xy_order_parameter_imag" + std::to_string(beta_set_[replica_id]) + mc_params.filename_postfix,
                   complex_xy_order_parameter_imag_set_[replica_id]);
        }
      }
    }
  }
  dump_data_timer.PrintElapsed();
}

template<size_t DimDof, size_t NumOfCouplingType>
void ExchangeMCExecutor<DimDof, NumOfCouplingType>::StatisticAndDumpResults_() {
  gqten::Timer statistic_timer("statistic");
  results_.beta = beta_set_;
  size_t sweeps = mc_params.sweeps;
  size_t warmup_sample_num = mc_params.warmup_sample_num;
  size_t N = phys_params.N;
  size_t N_square = N * N;
  size_t N_quartic = N_square * N_square;
#pragma omp parallel for default(none) \
                shared(mc_params, phys_params, N, N_square, N_quartic, sweeps, warmup_sample_num, results_)\
                num_threads(mc_params.thread_num)\
                schedule(static)
  for (size_t replica_id = 0; replica_id < mc_params.thread_num; replica_id++) {
    double beta = beta_set_[replica_id];

    auto thermalized_energy_data =
        std::vector(energy_set_[replica_id].cbegin() + warmup_sample_num, energy_set_[replica_id].cend());
    double total_energy = Mean(thermalized_energy_data);
    results_.energy[replica_id] = total_energy / phys_params.N;
    results_.specific_heat[replica_id] = Variance(thermalized_energy_data, total_energy) * beta * beta / N;

    if (DimDof > 1 && interaction_isotropy_) {
      auto thermalized_stiffness_data =
          std::vector(stiffness_set_[replica_id].cbegin() + warmup_sample_num, stiffness_set_[replica_id].cend());
      results_.stiffness[replica_id] = Mean(thermalized_stiffness_data);
    }

    if (phys_params.geometry == "Honeycomb" || phys_params.geometry == "Honeycomb3") {

      //anti-ferromagnetic magnetization
      std::array<std::vector<double>, DimDof> thermalized_magnetization_component_square_Q1_data,
          thermalized_magnetization_component_square_Q2_data, thermalized_magnetization_component_square_Q3_data;

      double zigzag_af_m2_q1(0.0), zigzag_af_m2_q2(0.0), zigzag_af_m2_q3(0.0);
      for (size_t dof = 0; dof < DimDof; dof++) {
        thermalized_magnetization_component_square_Q1_data[dof] =
            std::vector<double>(sweeps - warmup_sample_num, 0.0);
        for (size_t i = 0; i < thermalized_magnetization_component_square_Q1_data[dof].size(); i++) {
          thermalized_magnetization_component_square_Q1_data[dof][i] =
              zigzag_af_magnetization1_set_[replica_id][dof][i + warmup_sample_num]
                  * zigzag_af_magnetization1_set_[replica_id][dof][i + warmup_sample_num];
        }
        double m2_q1_component = Mean(thermalized_magnetization_component_square_Q1_data[dof]);
        results_.zigzag_af1_magnetization_square[dof][replica_id] = m2_q1_component / N_square;
        zigzag_af_m2_q1 += m2_q1_component;
      }
      for (size_t dof = 0; dof < DimDof; dof++) {
        thermalized_magnetization_component_square_Q2_data[dof] =
            std::vector<double>(sweeps - warmup_sample_num, 0.0);
        for (size_t i = 0; i < thermalized_magnetization_component_square_Q2_data[dof].size(); i++) {
          thermalized_magnetization_component_square_Q2_data[dof][i] =
              zigzag_af_magnetization2_set_[replica_id][dof][i + warmup_sample_num]
                  * zigzag_af_magnetization2_set_[replica_id][dof][i + warmup_sample_num];
        }
        double m2_q2_component = Mean(thermalized_magnetization_component_square_Q2_data[dof]);
        results_.zigzag_af2_magnetization_square[dof][replica_id] = m2_q2_component / N_square;
        zigzag_af_m2_q2 += m2_q2_component;
      }
      for (size_t dof = 0; dof < DimDof; dof++) {
        thermalized_magnetization_component_square_Q3_data[dof] =
            std::vector<double>(sweeps - warmup_sample_num, 0.0);
        for (size_t i = 0; i < thermalized_magnetization_component_square_Q3_data[dof].size(); i++) {
          thermalized_magnetization_component_square_Q3_data[dof][i] =
              zigzag_af_magnetization3_set_[replica_id][dof][i + warmup_sample_num]
                  * zigzag_af_magnetization3_set_[replica_id][dof][i + warmup_sample_num];
        }
        double m2_q3_component = Mean(thermalized_magnetization_component_square_Q3_data[dof]);
        results_.zigzag_af3_magnetization_square[dof][replica_id] = m2_q3_component / N_square;
        zigzag_af_m2_q3 += m2_q3_component;
      }
      //anti-ferromagnetic susceptibility
      for (size_t dof = 0; dof < DimDof; dof++) {
        std::vector<double> thermalized_sum_magnetization_abs_data(sweeps - warmup_sample_num, 0.0);
        for (size_t j = 0; j < thermalized_sum_magnetization_abs_data.size(); j++) {
          thermalized_sum_magnetization_abs_data[j] =
              abs(zigzag_af_magnetization1_set_[replica_id][dof][j + warmup_sample_num])
                  + abs(zigzag_af_magnetization2_set_[replica_id][dof][j + warmup_sample_num])
                  + abs(zigzag_af_magnetization3_set_[replica_id][dof][j + warmup_sample_num]);
        }
        results_.zigzag_magnetization_susceptibility[dof][replica_id] =
            Variance(thermalized_sum_magnetization_abs_data) / N * beta;
      }
      //binder ratio
      std::vector<double> thermalized_magnetization_square_q1_data(sweeps - warmup_sample_num, 0.0),
          thermalized_magnetization_square_q2_data(sweeps - warmup_sample_num, 0.0),
          thermalized_magnetization_square_q3_data(sweeps - warmup_sample_num, 0.0);
      std::vector<double> thermalized_magnetization_quartic_q1_data(sweeps - warmup_sample_num, 0.0),
          thermalized_magnetization_quartic_q2_data(sweeps - warmup_sample_num, 0.0),
          thermalized_magnetization_quartic_q3_data(sweeps - warmup_sample_num, 0.0);
      for (size_t i = 0; i < sweeps - warmup_sample_num; i++) {
        double sum(0.0);
        for (size_t dof = 0; dof < DimDof; dof++) {
          sum += thermalized_magnetization_component_square_Q1_data[dof][i];
        }
        thermalized_magnetization_square_q1_data[i] = sum;
        thermalized_magnetization_quartic_q1_data[i] = sum * sum;
      }
      double zigzag_af_m4_q1 = Mean(thermalized_magnetization_quartic_q1_data);
      results_.binder_ratio_zigzag_direction1[replica_id] = zigzag_af_m4_q1 / zigzag_af_m2_q1 / zigzag_af_m2_q1;

      for (size_t i = 0; i < sweeps - warmup_sample_num; i++) {
        double sum(0.0);
        for (size_t dof = 0; dof < DimDof; dof++) {
          sum += thermalized_magnetization_component_square_Q2_data[dof][i];
        }
        thermalized_magnetization_square_q2_data[i] = sum;
        thermalized_magnetization_quartic_q2_data[i] = sum * sum;
      }
      double zigzag_af_m4_q2 = Mean(thermalized_magnetization_quartic_q2_data);
      results_.binder_ratio_zigzag_direction2[replica_id] = zigzag_af_m4_q2 / zigzag_af_m2_q2 / zigzag_af_m2_q2;

      for (size_t i = 0; i < sweeps - warmup_sample_num; i++) {
        double sum(0.0);
        for (size_t dof = 0; dof < DimDof; dof++) {
          sum += thermalized_magnetization_component_square_Q3_data[dof][i];
        }
        thermalized_magnetization_square_q3_data[i] = sum;
        thermalized_magnetization_quartic_q3_data[i] = sum * sum;
      }
      double zigzag_af_m4_q3 = Mean(thermalized_magnetization_quartic_q3_data);
      results_.binder_ratio_zigzag_direction3[replica_id] = zigzag_af_m4_q3 / zigzag_af_m2_q3 / zigzag_af_m2_q3;

      std::vector<double>
          thermalized_magnetization_square_data(sweeps - warmup_sample_num, 0.0), //m^2 = m_q1^2 + m_q2^2 + m_q3^2
      thermalized_magnetization_quartic_data(sweeps - warmup_sample_num, 0.0);
      for (size_t i = 0; i < sweeps - warmup_sample_num; i++) {
        thermalized_magnetization_square_data[i] =
            thermalized_magnetization_square_q1_data[i] + thermalized_magnetization_square_q2_data[i]
                + thermalized_magnetization_square_q3_data[i];
        thermalized_magnetization_quartic_data[i] =
            thermalized_magnetization_square_data[i] * thermalized_magnetization_square_data[i];
      }
      double zigzag_af_m4 = Mean(thermalized_magnetization_quartic_data);
      double zigzag_af_m2 = Mean(thermalized_magnetization_square_data);
      results_.binder_ratio_zigzag_af[replica_id] = zigzag_af_m4 / zigzag_af_m2 / zigzag_af_m2;
      results_.zigzag_magnetization_susceptibility_b[replica_id] =
          Variance(thermalized_magnetization_square_data) / N_square * beta;

      //C3 symmetry order parameter and its susceptibility
      using complexT = std::complex<double>;
      std::vector<complexT> thermalized_data_of_c3_order_parameter(sweeps - warmup_sample_num, 0.0);
      std::vector<double> thermalized_data_of_c3_order_parameter_abs(sweeps - warmup_sample_num, 0.0);
      std::vector<double> thermalized_data_of_c3_order_parameter_square(sweeps - warmup_sample_num, 0.0);
      for (size_t i = 0; i < sweeps - warmup_sample_num; i++) {
        thermalized_data_of_c3_order_parameter[i] = thermalized_magnetization_square_q1_data[i]
            + std::exp(complexT(0, 2 * M_PI / 3)) * thermalized_magnetization_square_q2_data[i]
            + std::exp(complexT(0, 4 * M_PI / 3)) * thermalized_magnetization_square_q3_data[i];
        thermalized_data_of_c3_order_parameter_abs[i] = abs(thermalized_data_of_c3_order_parameter[i]);
        thermalized_data_of_c3_order_parameter_square[i] =
            thermalized_data_of_c3_order_parameter_abs[i] * thermalized_data_of_c3_order_parameter_abs[i];
      }
      double B_c3_square = Mean(thermalized_data_of_c3_order_parameter_square);
      double B_c3_abs = Mean(thermalized_data_of_c3_order_parameter_abs);
      results_.c3_order_parameter_square[replica_id] = B_c3_square / N_quartic;
      results_.c3_order_parameter_susceptibility[replica_id] =
          (B_c3_square - B_c3_abs * B_c3_abs) / N_square * beta;
      results_.binder_ratio_c3[replica_id] = B_c3_square / B_c3_abs / B_c3_abs;

      if (has_gs_data) {
        std::vector<double> thermalized_data_of_complex_order_parameter_square(sweeps - warmup_sample_num, 0.0);
        std::vector<double> thermalized_data_of_complex_order_parameter_quartic(sweeps - warmup_sample_num, 0.0);
        for (size_t i = 0; i < sweeps - warmup_sample_num; i++) {
          auto &real = complex_xy_order_parameter_real_set_[replica_id][i + warmup_sample_num];
          auto &imag = complex_xy_order_parameter_imag_set_[replica_id][i + warmup_sample_num];
          thermalized_data_of_complex_order_parameter_square[i] = real * real + imag * imag;
          thermalized_data_of_complex_order_parameter_quartic[i] =
              thermalized_data_of_complex_order_parameter_square[i]
                  * thermalized_data_of_complex_order_parameter_square[i];
        }
        double mxy_square = Mean(thermalized_data_of_complex_order_parameter_square);
        double mxy_quartic = Mean(thermalized_data_of_complex_order_parameter_quartic);
        results_.complex_order_parameter[replica_id] = mxy_square;
        results_.binder_ratio_complex_order_parameter[replica_id] = mxy_quartic / mxy_square / mxy_square;
      }
    }
  }
  statistic_timer.PrintElapsed();
  results_.DumpData("results" + mc_params.filename_postfix);
}

#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_ExchangeMCEXECUTOR_H
