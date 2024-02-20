#include <iostream>
#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <time.h>
#include <vector>
#include <stdlib.h>     // system
#include "params_case.h"
#include "DefSpinOne.h"
#include "myutil.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env(mpi::threading::multiple);
  mpi::communicator world;

  CaseParams params(argv[1]);
  const size_t Lx = params.Lx;
//  const size_t Ly = params.Ly;
  const size_t N = 2 * Lx * params.Ly;
  cout << "The total number of sites: " << N << endl;
  double J1 = params.J1;
  double J2 = params.J2;
  double J3 = params.J3;
  cout << "Model parameter: J1 = " << J1 << ",\t J2 = " << J2 << ",\t J3 = " << J3 << endl;
  clock_t startTime, endTime;
  startTime = clock();

  std::vector<size_t> input_D_set;
  bool has_bond_dimension_parameter = ParserBondDimension(
      argc, argv,
      input_D_set);

  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  qlmps::FiniteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter)
  );

  bool noise_valid(false);
  for (size_t i = 0; i < params.noise.size(); i++) {
    if (params.noise[i] != 0) {
      noise_valid = true;
      break;
    }
  }

  double e0(0.0); //energy

  //const size_t Ly = params.Ly;
  using namespace spin_one_model;
  const SiteVec<TenElemT, U1QN> sites(N, pb_out);
  qlmps::MPOGenerator<TenElemT, U1QN> mpo_gen(sites, qn0);

  qlmps::MPO<Tensor> mpo(N);
  if (IsPathExist(qlmps::kMpoPath)) {
    for (size_t i = 0; i < mpo.size(); i++) {
      std::string filename = qlmps::kMpoPath + "/" +
          qlmps::kMpoTenBaseName + std::to_string(i) + "." + qlmps::kQLTenFileSuffix;
      mpo.LoadTen(i, filename);
    }
    cout << "FiniteMPO loaded." << endl;
  } else {
    cout << "No mpo directory. exiting" << std::endl;
    exit(0);
  }

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1QN>;
  FiniteMPST mps(sites);

  std::vector<long unsigned int> stat_labs(N);

  for (size_t i = 0; i < N - N % 3; ++i) {
    stat_labs[i] = i % 3;
  }
  for (size_t i = N - N % 3; i < N; ++i) {
    stat_labs[i] = 1;
  }

  if (world.rank() == 0) {
    if (IsPathExist(kMpsPath)) {
      if (N == GetNumofMps()) {
        cout << "The number of mps files is consistent with mps size." << endl;
        cout << "Directly use mps from files." << endl;
      } else {
        qlmps::DirectStateInitMps(mps, stat_labs);
        cout << "Initial mps as direct product state." << endl;
        mps.Dump(sweep_params.mps_path, true);
      }
    } else {
      qlmps::DirectStateInitMps(mps, stat_labs);
      cout << "Initial mps as direct product state." << endl;
      mps.Dump(sweep_params.mps_path, true);
    }

  }

  if (!has_bond_dimension_parameter) {
    if (world.size() == 1) {
      e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params);
    } else {
      qlmps::FiniteVMPSSweepParams sweep_params(
          params.Sweeps,
          params.Dmin, params.Dmax, params.CutOff,
          qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter)
      );
      e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params, world);
    }
  } else {
    size_t DMRG_time = input_D_set.size();
    std::vector<size_t> MaxLanczIterSet(DMRG_time);
    MaxLanczIterSet.back() = params.MaxLanczIter;
    if (DMRG_time > 1) {
      size_t MaxLanczIterSetSpace;
      MaxLanczIterSet[0] = 3;
      MaxLanczIterSetSpace = (params.MaxLanczIter - 3) / (DMRG_time - 1);
      std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << ", ";
      for (size_t i = 1; i < DMRG_time - 1; i++) {
        MaxLanczIterSet[i] = MaxLanczIterSet[i - 1] + MaxLanczIterSetSpace;
        std::cout << MaxLanczIterSet[i] << ", ";
      }
      std::cout << MaxLanczIterSet.back() << "]" << std::endl;
    } else {
      std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << "]" << std::endl;
    }

    if (world.size() == 1) {
      for (size_t i = 0; i < DMRG_time; i++) {
        size_t D = input_D_set[i];
        std::cout << "D_max = " << D << std::endl;
        qlmps::FiniteVMPSSweepParams sweep_params(
            params.Sweeps,
            D, D, params.CutOff,
            qlmps::LanczosParams(params.LanczErr, MaxLanczIterSet[i])
        );

        e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params);
      }
    } else {
      for (size_t i = 0; i < DMRG_time; i++) {
        size_t D = input_D_set[i];
        std::cout << "D_max = " << D << std::endl;
        qlmps::FiniteVMPSSweepParams sweep_params(
            params.Sweeps,
            D, D, params.CutOff,
            qlmps::LanczosParams(params.LanczErr, MaxLanczIterSet[i])
        );
        e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params, world);
      }
    }
  }
  std::cout << "E0/site: " << e0 / N << std::endl;

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
