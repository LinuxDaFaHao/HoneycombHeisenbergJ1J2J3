#include "DefSpinOne.h"
#include "operators.h"
#include <time.h>
#include <vector>
#include <stdlib.h>     // system
#include "qlmps/qlmps.h"
#include "myutil.h"
#include "two_site_update_noised_finite_vmps_mpi_impl3.h"
#include "params_case.h"

using namespace qlmps;
using namespace qlten;
using namespace std;
using namespace spin_one_model;

int ParserContinueArgs(const int argc, char *argv[],
                       size_t &start_site,
                       char &start_direction) {
  int nOptionIndex = 1;

  string arguement1 = "--start_site=";
  string arguement2 = "--start_direction=";
  bool site_argument_has(false), direction_argument_has(false);
  while (nOptionIndex < argc) {
    if (strncmp(argv[nOptionIndex], arguement1.c_str(), arguement1.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement1.size()];
      start_site = atoi(para_string.c_str());
      site_argument_has = true;
    } else if (strncmp(argv[nOptionIndex], arguement2.c_str(), arguement2.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement2.size()];
      start_direction = *(para_string.c_str());
      direction_argument_has = true;
    }
    nOptionIndex++;
  }

  if (!site_argument_has) {
    std::cout << "Note: no start site argument. exit." << std::endl;
    exit(1);
  }

  if (!direction_argument_has) {
    std::cout << "Note: no thread argument. exit" << std::endl;
    exit(2);
  }

  return 0;

}
int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env(mpi::threading::multiple);
  if (env.thread_level() < mpi::threading::multiple) {
    std::cout << "thread level of env is not right." << std::endl;
    env.abort(-1);
  }
  mpi::communicator world;
  CaseParams params(argv[1]);

  size_t start_site;
  char start_direction;

  ParserContinueArgs(argc, argv, start_site, start_direction);

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

  const SiteVec<TenElemT, U1QN> sites(N, pb_out);
  qlmps::MPOGenerator<TenElemT, U1QN> mpo_gen(sites, qn0);

  qlmps::MPO<Tensor> mpo(N);
  if (IsPathExist(kMpoPath)) {
    for (size_t i = 0; i < mpo.size(); i++) {
      std::string filename = kMpoPath + "/" +
          kMpoTenBaseName + std::to_string(i) + "." + kQLTenFileSuffix;
      mpo.LoadTen(i, filename);
    }
    cout << "FiniteMPO loaded." << endl;
  } else {
    cout << "No mpo directory. exiting" << std::endl;
    exit(0);
  }

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1QN>;
  FiniteMPST mps(sites);

  if (world.rank() == 0) {
    if (params.Threads > 2) {
      qlten::hp_numeric::SetTensorManipulationThreads(params.Threads - 2);
    } else {
      qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);
    }
  } else {
    qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);
  }
  qlmps::FiniteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter),
      params.noise
  );
  if (IsPathExist(kMpsPath)) {//mps only can be load from file
    if (N == GetNumofMps()) {
      cout << "The number of mps files is consistent with mps size." << endl;
      cout << "Directly use mps from files." << endl;
    } else {
      cout << "mps file number do not right" << endl;
      env.abort(-1);
    }
  } else {
    cout << " no mps file" << endl;
    env.abort(-1);
  }

  auto e0 = qlmps::TwoSiteFiniteVMPS2(mps, mpo, sweep_params, world, start_site, start_direction);
  if (world.rank() == 0) {
    std::cout << "E0/site: " << e0 / N << std::endl;
    endTime = clock();
    cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  }
  return 0;

}


