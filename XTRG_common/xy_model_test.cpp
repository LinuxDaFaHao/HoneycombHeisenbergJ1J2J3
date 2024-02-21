// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Haoxin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2022/5/12
*
*/

/**
@file xy_model_test.cpp
@brief test for XTRG in 1D XY model
*/


#include "qlmps/qlmps.h"
#include "qlten/qlten.h"

#include "xtrg.h"

using TenElemT = qlten::qlten_Double;
using QNT = qlten::QN<qlten::U1QNVal>;
using QNSctT = qlten::QNSector<QNT>;
using IndexT = qlten::Index<QNT>;
using Tensor = qlten::qltensor<TenElemT, QNT>;
using SiteVecT = qlmps::SiteVec<TenElemT, QNT>;
using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;

//forward declaration
double CalPartitionFunctionOpenXYChain(
    const double J,       //coupling constant
    const double L,       //chain length
    const double beta     // 1/T
);
int main() {
  auto qn0 =  QNT({qlten::QNCard("Sz", qlten::U1QNVal( 0))});
  auto qnup = QNT({qlten::QNCard("Sz", qlten::U1QNVal( 1))});
  auto qndn = QNT({qlten::QNCard("Sz", qlten::U1QNVal(-1))});
  auto pb_out = IndexT(
      {QNSctT(qnup, 1), QNSctT(qndn, 1)},
      qlten::qltenIndexDirType::OUT
  );
  auto pb_in = qlten::InverseIndex(pb_out);
  Tensor sz({pb_in, pb_out});
  Tensor sp({pb_in, pb_out});
  Tensor sm({pb_in, pb_out});
  sz(0, 0) =  0.5;
  sz(1, 1) = -0.5;
  sp(0, 1) = 1;
  sm(1, 0) = 1;
  size_t N = 10;
  double J = 1;
  SiteVecT sites(N, pb_out);
  qlmps::MPOGenerator<TenElemT, QNT> mpo_gen(sites, qn0);
  for (size_t i = 0; i < N-1; ++i) {
    mpo_gen.AddTerm(J/2, sp, i, sm, i+1);
    mpo_gen.AddTerm(J/2, sm, i, sp, i+1);
  }
  auto mpo = mpo_gen.Gen();
  using MPOT = xtrg::FiniteMPO<TenElemT, QNT>;
  MPOT hamiltonian(N), density_matrix(N);
  for(size_t i = 0; i < N; i++) {
    hamiltonian(i) = new Tensor(mpo[i]);
    density_matrix(i) = new Tensor();
  }

  double tau = 0.001;
  double M = 10;

  size_t Dmin = 16;
  size_t Dmax = 48;
  double TruncErr = 1e-7;
  size_t TaylorOrder = 30;
  double TaylorErr = 1e-13;
  size_t Sweeps = 10;
  double ConvergeError = 1e-9;

  xtrg::XtrgParams params(
      tau, M,
      Dmin, Dmax, TruncErr,
      TaylorOrder, TaylorErr,
      Sweeps, ConvergeError
      );

  size_t TenNumThreads = 4;

  qlten::hp_numeric::SetTensorTransposeNumThreads(TenNumThreads);
  qlten::hp_numeric::SetTensorManipulationThreads(TenNumThreads );


//  xtrg::InitializeDensityMatrix(hamiltonian, params, density_matrix);
  xtrg::Xtrg(hamiltonian, params, density_matrix);
  std::cout << "\n" << std::endl;
  double beta = tau;
  for(size_t i = 0; i < M+1; i++) {
    beta = beta*2;
    CalPartitionFunctionOpenXYChain(J, N, beta);
  }
//  std::cout << "E0/site: " << e0 / N << std::endl;    // E0/site: -0.43412
  return 0;
}



double CalPartitionFunctionOpenXYChain(
    const double J,       //coupling constant
    const double L,       //chain length
    const double beta     // 1/T
) {
  double Z = 1.0;
  double Fn = 0;
  for(size_t k = 1; k <=L ; k++) {
    double epsilon_k = -J * std::cos( (k*M_PI)/(L+1) );
    Z *= ( std::exp(-beta*epsilon_k) + 1 );
    Fn = Fn - std::log( std::exp(-beta*epsilon_k) + 1 );
  }
  Fn = Fn / beta;
  std::cout << "beta = " << beta << ", Z = " << Z << ", Fn =" << Fn << std::endl;
  return Z;
}