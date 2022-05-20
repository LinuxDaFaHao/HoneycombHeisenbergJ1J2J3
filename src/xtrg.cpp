// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Haoxin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2022/5/13
*
*/

/**
 * @file xtrg.cpp
 * @brief the main file for calculate the thermodynamic of J1-J2-J3 Heisenberg model on Honeycomb lattice
 */


#include "gqmps2/gqmps2.h"
#include "gqten/gqten.h"
#include "../XTRG_common/xtrg.h"

#include "params_case.h"
#include "DefSpinOne.h"
#include "myutil.h"


//TODO: measurement
using namespace std;

int main(int argc, char *argv[]) {


  CaseParams params(argv[1]);
  const size_t Lx = params.Lx;
//  const size_t Ly = params.Ly;
  const size_t N = 2 * Lx * params.Ly;
  std::cout << "The total number of sites: " << N << std::endl;
  double J1 = params.J1;
  double J2 = params.J2;
  double J3 = params.J3;
  std::cout << "Model parameter: J1 = " << J1 << ",\t J2 = " << J2 << ",\t J3 = " << J3 << endl;
  clock_t startTime,endTime;
  startTime = clock();
  gqten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
  gqten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  using namespace spin_one_model;
  const gqmps2::SiteVec<TenElemT, U1QN> sites(N, pb_out);
  gqmps2::MPOGenerator<TenElemT, U1QN> mpo_gen(sites, qn0);

  xtrg::FiniteMPO<TenElemT, U1QN> hamiltonian(N), density_matrix(N);
  if (gqmps2::IsPathExist(kMpoPath)) {
    for(size_t i=0; i < N;i++){
      std::string filename = kMpoPath + "/" +
          kMpoTenBaseName + std::to_string(i) + "." + gqmps2::kGQTenFileSuffix;
      hamiltonian.LoadTen(i,filename);
      density_matrix.alloc(i);
    }
    cout << "FiniteMPO loaded." << endl;
  }else{
    cout << "No mpo directory. exiting" << std::endl;
    exit(0);
  }

  xtrg::XtrgParams xtrg_params(
      params.tau, params.M,
      params.Dmin, params.Dmax, params.CutOff,
      params.TaylorOrder, params.TaylorErr,
      params.Sweeps,params.ConvergeTolerance
  );


  xtrg::Xtrg(hamiltonian, xtrg_params, density_matrix);
  endTime = clock();
  cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}


