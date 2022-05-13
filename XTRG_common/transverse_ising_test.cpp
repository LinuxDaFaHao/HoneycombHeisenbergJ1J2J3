// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Haoxin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2022/5/12
*
*/

/**
@file transverse_ising_test.cpp
@brief test for XTRG in 1D transverse ising model
*/


#include "gqmps2/gqmps2.h"
#include "gqten/gqten.h"
using TenElemT = gqten::GQTEN_Double;
using QNT = gqten::QN<gqten::U1QNVal>;
using QNSctT = gqten::QNSector<QNT>;
using IndexT = gqten::Index<QNT>;
using Tensor = gqten::GQTensor<TenElemT, QNT>;
using SiteVecT = gqmps2::SiteVec<TenElemT, QNT>;
using FiniteMPST = gqmps2::FiniteMPS<TenElemT, QNT>;
int main() {
  auto qn0 =  QNT({gqten::QNCard("Sz", gqten::U1QNVal( 0))});
  auto qnup = QNT({gqten::QNCard("Sz", gqten::U1QNVal( 1))});
  auto qndn = QNT({gqten::QNCard("Sz", gqten::U1QNVal(-1))});
  auto pb_out = IndexT(
      {QNSctT(qnup, 1), QNSctT(qndn, 1)},
      gqten::GQTenIndexDirType::OUT
  );
  auto pb_in = gqten::InverseIndex(pb_out);
  Tensor sz({pb_in, pb_out});
  Tensor sp({pb_in, pb_out});
  Tensor sm({pb_in, pb_out});
  sz(0, 0) =  0.5;
  sz(1, 1) = -0.5;
  sp(0, 1) = 1;
  sm(1, 0) = 1;
  size_t N = 20;
  SiteVecT sites(N, pb_out);
  gqmps2::MPOGenerator<TenElemT, QNT> mpo_gen(sites, qn0);
  for (size_t i = 0; i < N-1; ++i) {
    mpo_gen.AddTerm(1.0, sz, i, sz, i+1);
    mpo_gen.AddTerm(0.5, sp, i, sm, i+1);
    mpo_gen.AddTerm(0.5, sm, i, sp, i+1);
  }
  auto mpo = mpo_gen.Gen();
  size_t Sweeps = 10;
  size_t Dmin = 16;
  size_t Dmax = 24;
  double TruncErr = 1e-7;
  double LanczErr = 1e-8;
  size_t MaxLanczIter = 100;
  gqmps2::SweepParams sweep_params(
      Sweeps,
      Dmin, Dmax, TruncErr,
      gqmps2::LanczosParams(LanczErr, MaxLanczIter)
  );
  FiniteMPST mps(sites);
  std::vector<size_t> stat_labs;
  for (size_t i = 0; i < N; ++i) { stat_labs.push_back(i % 2); }
  gqmps2::DirectStateInitMps(mps, stat_labs);
  mps.Dump(sweep_params.mps_path, true);
  size_t TenNumThreads = 4;

  gqten::hp_numeric::SetTensorTransposeNumThreads(TenNumThreads);
  gqten::hp_numeric::SetTensorManipulationThreads(TenNumThreads );

  FiniteMPST mps2 = mps;
  auto e0 = gqmps2::TwoSiteFiniteVMPS(mps, mpo, sweep_params);
  std::cout << "E0/site: " << e0 / N << std::endl;    // E0/site: -0.43412
  return 0;
}