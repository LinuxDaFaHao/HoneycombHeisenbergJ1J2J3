//
// Created by Hao-Xin on 2022/10/2.
//

#ifndef HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCKMCEXECUTOR_H
#define HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCKMCEXECUTOR_H

#include <string>
#include "clock_config.h"
#include "gqten/framework/bases/executor.h"       //Executor
#include "gqten/utility/timer.h"
#include "../MonteCarlo_src/common.h"

struct MCParams {
  size_t sweeps;
  size_t cluster_radius;
  size_t sample_interval; // for how many times of wolf update between two samples. 0 means sweep lattice
  size_t print_interval;
  std::string filename_postfix;
};

struct ClockPhysParams {
  ClockPhysParams(const std::string &geometry,
                  size_t lx,
                  size_t ly,
                  double beta,
                  const std::array<CouplingStructure<clock_dim_dof>, clock_interaction_range> &coupling_structures,
                  double rand_field_magnitude)
      : geometry(geometry),
        Lx(lx),
        Ly(ly),
        beta(beta),
        coupling_structures(coupling_structures),
        rand_field_magnitude(rand_field_magnitude) {
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
  std::array<CouplingStructure<clock_dim_dof>, clock_interaction_range> coupling_structures;
  double rand_field_magnitude;
};

template<size_t q>
class ClockMCExecutor : public gqten::Executor {
  using LocalDOFT = ClockDOF<q>;
 public:
  ClockMCExecutor(const MCParams &,
                  const ClockPhysParams &,
                  const size_t world_rank
  );
  void Execute() override;

  MCParams mc_params;
  ClockPhysParams phys_params;
 private:
  void WolfSweep_();
  void WolfRandomSweep_(size_t sample_interval);
  void WolfGrowCluster_(const ClockDOF<q> &, const size_t);
  bool WolfFlipCluster_(const ClockDOF<q> &);

  void MetropolisSweep_();
  bool MetropolisSingleStep_(const size_t site);

  void Measure_(size_t sweep);
  void StatisticAndDumpData_();

  int mpi_world_rank_;
  size_t geometry_id_;
  ClockConfig<q> config_;
  SWCluster cluster_;
  LatticeLink<clock_dim_dof, clock_interaction_range> *plattice_link_;
  RandomField random_field_x_;
  RandomField random_field_y_;

  double b_; //The probability used to correct the deviation from the restriction of cluster size
  size_t flip_times_;

  /// result data
  std::vector<double> energy_;
  std::vector<double> fm_m_square_;

  std::vector<double> res_;
};

template<size_t q>
ClockMCExecutor<q>::ClockMCExecutor(const MCParams &mc_params,
                                    const ClockPhysParams &phys_params,
                                    size_t world_rank):
    gqten::Executor(),
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

  config_.Random();


  if (phys_params.rand_field_magnitude > 0.0) {
    double h = phys_params.rand_field_magnitude;
    random_field_x_.resize(phys_params.N);
    random_field_y_.resize(phys_params.N);
    std::uniform_real_distribution<double> u_random_field(-h, h);
    for (size_t site = 0; site < phys_params.N; site++) {
      double h_x(0.0), h_y(0.0);
      // A uniform distribution within a circular with radius h
      do {
        h_x = u_random_field(random_engine);
        h_y = u_random_field(random_engine);
      } while (h_x * h_x + h_y * h_y > h * h);

      random_field_x_[site] = h_x;
      random_field_y_[site] = h_y;
    }
  }

  flip_times_ = 0;

  size_t sweeps = mc_params.sweeps;
  energy_.reserve(sweeps);
  fm_m_square_.reserve(sweeps);

  SetStatus(gqten::INITED);
}

template<size_t q>
void ClockMCExecutor<q>::Execute() {
  SetStatus(gqten::EXEING);
  gqten::Timer wolf_mc_execute_timer("wolf_mc_execute");
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
                << "time = " << execute_time
                << std::endl;
    }
  }
  if (mpi_world_rank_ == 0) {
    std::cout << "Flip probability: " << double(flip_times_) / double(mc_params.sweeps * mc_params.sample_interval)
              << std::endl;
  }

  StatisticAndDumpData_();
  SetStatus(gqten::FINISH);
}

template<size_t q>
void ClockMCExecutor<q>::WolfSweep_() {
  ClockDOF<q> axis;
  for (size_t starting_site = 0; starting_site < phys_params.N; starting_site++) {
    axis.Random();
    WolfGrowCluster_(axis, starting_site);
    if (WolfFlipCluster_(axis)) flip_times_++;
  }
}

template<size_t q>
void ClockMCExecutor<q>::WolfRandomSweep_(size_t sample_interval) {
  ClockDOF<q> axis;
  std::uniform_int_distribution<int> u_int(0, phys_params.N - 1);
  for (size_t i = 0; i < sample_interval; i++) {
    axis.Random();
    WolfGrowCluster_(axis, u_int(random_engine));
    if (WolfFlipCluster_(axis)) flip_times_++;
  }
}

template<size_t q>
void ClockMCExecutor<q>::WolfGrowCluster_(const ClockDOF<q> &axis,
                                          const size_t starting_site) {
  cluster_ = {starting_site};
  std::vector<bool> in_cluster(phys_params.N, false);
  in_cluster[starting_site] = true;
  double k_eff = phys_params.coupling_structures[1].coupling_const[0][0] * phys_params.beta;

  b_ = 1.0;
  if (mc_params.cluster_radius == 0) { // no cluster size restriction
    for (size_t j = 0; j < cluster_.size(); j++) {
      const std::vector<size_t> &site2_set = plattice_link_->GetSitesLinkedTo(cluster_[j], 1);
      for (auto site2: site2_set) {
        if (!in_cluster[site2] && config_.Active(cluster_[j], site2, k_eff, axis)) {
          in_cluster[site2] = true;
          cluster_.push_back(site2);
        }
      }
    }
  } else {
    for (size_t j = 0; j < cluster_.size(); j++) {
      const std::vector<size_t> &site2_set = plattice_link_->GetSitesLinkedTo(cluster_[j], 1);
      for (auto site2: site2_set) {
        if (plattice_link_->Distance(j, site2) <= mc_params.cluster_radius) {
          if (!in_cluster[site2] && config_.Active(cluster_[j], site2, k_eff, axis)) {
            in_cluster[site2] = true;
            cluster_.push_back(site2);
          }
        } else {
          b_ *= config_.NotActiveProbability(cluster_[j], site2, k_eff, axis);
        }
      }
    }
  }

}

template<size_t q>
bool ClockMCExecutor<q>::WolfFlipCluster_(const ClockDOF<q> &axis) {

  std::uniform_real_distribution<double> u(0, 1);
  if (phys_params.rand_field_magnitude <= 0.0 && b_ >= 1.0) {
    config_.Flip2(axis, cluster_);
    return true;
  } else if (phys_params.rand_field_magnitude <= 0.0 && b_ < 1.0) {
    if (u(random_engine) > b_) return false;
    else {
      config_.Flip2(axis, cluster_);
      return true;
    }
  }
//  if (b_ < 1.0 && u(random_engine) > b_) { //correction for the limitation of the cluster size
//    return false;
//  }

  double delta_e = config_.EnergyDifferenceOfRandomFieldFlipCluster(axis, cluster_, random_field_x_, random_field_y_);
//  if (delta_e <= 0.0) {
//    config_.Flip2(axis, cluster_);
//    return true;
//  } else {
  if (u(random_engine) <= std::exp(-phys_params.beta * delta_e) * b_) {
    config_.Flip2(axis, cluster_);
    return true;
  } else {
    return false;
  }
//  }
}

template<size_t q>
void ClockMCExecutor<q>::MetropolisSweep_() {
  for (size_t i = 0; i < phys_params.N; i++) {
    MetropolisSingleStep_(i);
  }
}

template<size_t q>
bool ClockMCExecutor<q>::MetropolisSingleStep_(const size_t site) {
  ClockDOF<q> flipto;
  flipto.Random();
  double delta_e = config_.EnergyDifferenceFlipSite(site, flipto, *plattice_link_, random_field_x_, random_field_y_);
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

template<size_t q>
void ClockMCExecutor<q>::Measure_(size_t sweep) {
  energy_.push_back(config_.Energy(*plattice_link_, random_field_x_, random_field_y_));
  fm_m_square_.push_back(config_.FerromagneticMSquare());
}

template<size_t q>
void ClockMCExecutor<q>::StatisticAndDumpData_() {
  gqten::Timer dump_data_timer("dump_data");
  DumpData("energy" + mc_params.filename_postfix, energy_);
  DumpData("msquare" + mc_params.filename_postfix, fm_m_square_);
  dump_data_timer.PrintElapsed();

  gqten::Timer statistic_timer("statistic");
  size_t sweeps = mc_params.sweeps;
  // energy
  auto half_data_of_energy = std::vector(energy_.begin() + mc_params.sweeps / 2, energy_.end());
  res_.push_back(Mean(half_data_of_energy) / phys_params.N);
  // specific heat
  res_.push_back(
      Variance(half_data_of_energy, res_[0] * phys_params.N) * phys_params.beta * phys_params.beta / phys_params.N);
  //ferromagnetic magnetization
  auto half_data_of_msquare = std::vector(fm_m_square_.begin() + sweeps / 2, fm_m_square_.end());
  double M2 = Mean(half_data_of_msquare);
  res_.push_back(M2 / phys_params.N / phys_params.N);
  std::vector<double> half_data_of_mquartic(sweeps / 2);
  for (size_t i = 0; i < sweeps / 2; i++) {
    half_data_of_mquartic[i] = half_data_of_msquare[i] * half_data_of_msquare[i];
  }
  res_.push_back(Mean(half_data_of_mquartic) / M2 / M2);
  DumpData("summary" + mc_params.filename_postfix, res_);
  statistic_timer.PrintElapsed();
}

#endif //HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCKMCEXECUTOR_H
