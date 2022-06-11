//
// Created by Hao-Xin on 2022/5/25.
//

#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_WOLFMCEXECUTOR_H
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_WOLFMCEXECUTOR_H

#include <string>
#include <fstream>
#include "gqten/framework/bases/executor.h"
#include "gqten/utility/timer.h"
#include "swcluster.h"
#include "lattice_link.h"
#include "lattice_config.h"
#include "mpi.h"


template <typename DataType>
void DumpData(
    const std::string &filename,
    const std::vector<DataType> &data
);

struct MCParams {
  size_t sweeps;
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
      : geometry(geometry), Lx(lx), Ly(ly), N(n), beta(beta), coupling_structures(coupling_structures) {}

  PhysParams(const std::string &geometry,
             size_t lx,
             size_t ly,
             double beta,
             const std::array<CouplingStructure<DimDof>, NumOfCouplingType> &coupling_structures)
      : geometry(geometry), Lx(lx), Ly(ly), beta(beta), coupling_structures(coupling_structures) {
    if(geometry == "Honeycomb") {
      N = lx * ly * 2;
    } else if (geometry == "Square") {
      N = lx * ly;
    } else if( geometry == "SVW") {
      N = lx * ly;
    }
  }

  std::string geometry; //Honeycomb
  size_t Lx;
  size_t Ly;
  size_t N;
  double beta; // 1/T
  std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures;
};

template<size_t DimDof, size_t NumOfCouplingType>
class WolfMCExecutor : public gqten::Executor {
  using LocalDOFT = LocalDOF<DimDof>;
 public:
  WolfMCExecutor(const MCParams &,
                 const PhysParams<DimDof, NumOfCouplingType> &);
  void Execute() override;

  const std::vector<double>& GetEnergyData() const {
    return energy_;
  }

//  const std::vector<size_t>& GetClusterSizeData() const {
//    return cluster_size_;
//  }

  const std::array<std::vector<double>, DimDof>& GetMagneticSusceptibility() const {
    return sum_spin_;
  }

  ~WolfMCExecutor() {
    delete plattice_link_;
  }

  MCParams mc_params;
  PhysParams<DimDof, NumOfCouplingType> phys_params;


 private:

  void WolfSweep_();
  void WolfGrowCluster(const LocalDOFT&, const size_t);
  bool WolfFlipCluster(const LocalDOFT&);


  size_t geometry_id_; // 0: square, 1: honeycomb, 2: SVW
  bool interaction_isotropy_; // not consider SIA
  bool isotropy_;
  LatticeConfig<DimDof> config_;
  SWCluster cluster_;
  LatticeLink<DimDof, NumOfCouplingType>* plattice_link_;

  /// result data
  std::vector<double> energy_;
  std::array<std::vector<double>, DimDof> sum_spin_;
  std::vector<double> stiffness_;

  /*
  std::vector<double> ix_;  // for the spin stiffness
  std::vector<double> iy_;
  std::vector<size_t> cluster_size_;
   */
};



template<size_t DimDof, size_t NumOfCouplingType>
WolfMCExecutor<DimDof, NumOfCouplingType>::WolfMCExecutor(const MCParams &mc_params,
               const PhysParams<DimDof, NumOfCouplingType> &phys_params) :
    gqten::Executor(),
    mc_params(mc_params),
    phys_params(phys_params),
    config_(phys_params.N),
    cluster_() {
  if(phys_params.geometry == "Honeycomb") {
    plattice_link_ = new HoneyCombTorusLatticeLink(phys_params.Lx, phys_params.Ly, phys_params.coupling_structures);
    geometry_id_ = 1;
  } else if( phys_params.geometry == "Square") {
    plattice_link_ = new SquareTorusLatticeLink(phys_params.Lx, phys_params.Ly, phys_params.coupling_structures);
    geometry_id_ = 0;
  } else if( phys_params.geometry == "SVW") {
    plattice_link_ = new SVWLatticeLink(phys_params.N, phys_params.coupling_structures);
    geometry_id_ = 2;
  } else {
    std::cout << "do not support now. " << std::endl;
    exit(0);
  }

  interaction_isotropy_ = true;
  for(size_t i = 1; i < NumOfCouplingType; i++) {
    interaction_isotropy_ = interaction_isotropy_ && phys_params.coupling_structures[i].IsIsometry();
  }

  isotropy_ = interaction_isotropy_ && phys_params.coupling_structures[0].IsIsometry();

  config_.Random();
  size_t sweeps = mc_params.sweeps;
  energy_.reserve(sweeps);


  if(DimDof >=2 && interaction_isotropy_) {
    stiffness_.reserve(sweeps);
  }

  for(size_t i = 0; i < DimDof; i++) {
    sum_spin_[i].reserve(sweeps);
  }
  SetStatus(gqten::INITED);
}

template<size_t DimDof, size_t NumOfCouplingType>
void WolfMCExecutor<DimDof, NumOfCouplingType>::Execute() {
  SetStatus(gqten::EXEING);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  gqten::Timer wolf_mc_execute_timer("wolf_mc_execute");
  for(size_t sweep = 0; sweep < mc_params.sweeps; sweep++) {
    WolfSweep_();

    energy_.push_back(config_.Energy(*plattice_link_));
    for(size_t i = 0; i < DimDof; i++) {
      sum_spin_[i].push_back(config_.SumOverComponent(i));
    }

    if(DimDof >= 2 && interaction_isotropy_) {
      /*
      size_t Lx = phys_params.Lx;
      size_t Ly = phys_params.Ly;
      double ix(0.0), iy(0.0);
      for(size_t i = 0; i < Lx * Ly; i++){
        size_t x = i%(Lx);
        size_t y = i/(Lx);
        size_t Tx = (x + 1) % (Lx) + y * (Lx);
        size_t Ty = x + ( (y + 1) % Ly ) * (Lx);
        ix += config_.SinDiff(i, Tx);
        iy += config_.SinDiff(i, Ty);
      }
      ix_.push_back(ix);
      iy_.push_back(iy);
       */
      switch (geometry_id_) {
        case 1:
          stiffness_.push_back(config_.template StiffnessHoneycomb<NumOfCouplingType>(*plattice_link_, phys_params.beta));
          break;
        case 0:
          stiffness_.push_back(config_.template StiffnessSquare<NumOfCouplingType>(*plattice_link_, phys_params.beta));
          break;
      }
    }

    if( world_rank == 0 && sweep%mc_params.print_interval == 0 ) {
      double execute_time = wolf_mc_execute_timer.Elapsed();
      std::cout << "[ sweep = " << sweep << " ]"
                << "time = " << execute_time
                << std::endl;
    }
  }

  gqten::Timer dump_data_timer("dump_data");
  DumpData("energy" + mc_params.filename_postfix, energy_);
  for(size_t i = 0; i < DimDof; i++) {
    DumpData("sum_spin"
                    + std::to_string(i)
                    + mc_params.filename_postfix,
                    sum_spin_[i]);
  }
  if(DimDof >= 2 && interaction_isotropy_) {
    DumpData("stiffness" + mc_params.filename_postfix, stiffness_);
  }
  dump_data_timer.PrintElapsed();
  SetStatus(gqten::FINISH);
}
template<size_t DimDof, size_t NumOfCouplingType>
void WolfMCExecutor<DimDof, NumOfCouplingType>::WolfSweep_() {
  LocalDOF<DimDof> axis;
  for(size_t starting_site = 0; starting_site < phys_params.N; starting_site++) {
    axis.Random();
//    axis.Show();
//    config_.PrintSign(axis, phys_params.Lx);
    WolfGrowCluster(axis, starting_site);
    WolfFlipCluster(axis);
//    config_.PrintSign(axis, phys_params.Lx);
  }
}

template<size_t DimDof, size_t NumOfCouplingType>
void WolfMCExecutor<DimDof, NumOfCouplingType>::WolfGrowCluster(const LocalDOF<DimDof>& axis, const size_t starting_site) {
  cluster_ = {starting_site};//reserve?
  std::vector<bool> in_cluster(phys_params.N, false);
  in_cluster[starting_site] = true;
//  config_.PrintSign(axis, phys_params.Lx);
  // Note here the effective single ion anisotropy is not important.
  for(size_t inter_type = 1; inter_type < NumOfCouplingType; inter_type++) {
    double eff_coupling = axis * (phys_params.coupling_structures[inter_type] * axis);
    double Keff = eff_coupling * phys_params.beta;
    for(size_t j = 0; j < cluster_.size(); j++) {
      const std::vector<size_t>& site2_set = plattice_link_->GetSitesLinkedTo(cluster_[j], inter_type);
      for(auto site2 : site2_set) {
        if(!in_cluster[site2] && config_.Active(cluster_[j], site2, Keff, axis)) {
          in_cluster[site2] = true;
          cluster_.push_back(site2);
        }
      }
    }
  }

//#ifndef NDEBUG
//  auto cluster_sort = cluster_;
//  std::sort(cluster_sort.begin(), cluster_sort.end());
//  for(size_t i = 0; i < cluster_sort.size() - 1; i++) {
//    assert(cluster_sort[i] != cluster_sort[i+1]);
//  }
//#endif
}

/**
 *
 * @param axis
 * @return flipped or not
 */
template<size_t DimDof, size_t NumOfCouplingType>
bool WolfMCExecutor<DimDof, NumOfCouplingType>::WolfFlipCluster(const LocalDOF<DimDof> &axis) {

  if(isotropy_) {
    config_.Flip(axis, cluster_);
    return true;
  }
  double delta_e = config_.EnergyDifferenceFlipCluster(axis, cluster_, *plattice_link_);
#ifndef NDEBUG
  double e_i = config_.template Energy<NumOfCouplingType>(*plattice_link_);
  auto config_copy = config_;
  config_copy.Flip(axis, cluster_);
  double e_f = config_copy.template Energy<NumOfCouplingType>(*plattice_link_);
  double delta_e2 = e_f - e_i;
#endif
  if(delta_e <= 0.0) {
    config_.Flip(axis, cluster_);
    return true;
  } else {
    std::uniform_real_distribution<double> u(0, 1);
    if( u(random_engine) <= std::exp(- phys_params.beta * delta_e) ) {
      config_.Flip(axis, cluster_);
      return true;
    } else {
      return false;
    }
  }
}

template <typename DataType>
void DumpData(
    const std::string &filename,
    const std::vector<DataType> &data
    ) {
  std::ofstream ofs(filename, std::ofstream::binary);
  for(auto datum : data) {
    ofs << datum << '\n';
  }
  ofs << std::endl;
  ofs.close();
}

#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_WOLFMCEXECUTOR_H
