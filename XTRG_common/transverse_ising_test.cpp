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


#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
using TenElemT = qlten::QLTEN_Double;
using QNT = qlten::QN<qlten::U1QNVal>;
using QNSctT = qlten::QNSector<QNT>;
using IndexT = qlten::Index<QNT>;
using Tensor = qlten::QLTensor<TenElemT, QNT>;
using SiteVecT = qlmps::SiteVec<TenElemT, QNT>;
using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;
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
  size_t N = 20;
  SiteVecT sites(N, pb_out);
  qlmps::MPOGenerator<TenElemT, QNT> mpo_gen(sites, qn0);
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
  qlmps::SweepParams sweep_params(
      Sweeps,
      Dmin, Dmax, TruncErr,
      qlmps::LanczosParams(LanczErr, MaxLanczIter)
  );
  FiniteMPST mps(sites);
  std::vector<size_t> stat_labs;
  for (size_t i = 0; i < N; ++i) { stat_labs.push_back(i % 2); }
  qlmps::DirectStateInitMps(mps, stat_labs);
  mps.Dump(sweep_params.mps_path, true);
  size_t TenNumThreads = 4;

  qlten::hp_numeric::SetTensorTransposeNumThreads(TenNumThreads);
  qlten::hp_numeric::SetTensorManipulationThreads(TenNumThreads );

  FiniteMPST mps2 = mps;
  auto e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params);
  std::cout << "E0/site: " << e0 / N << std::endl;    // E0/site: -0.43412
  return 0;
}