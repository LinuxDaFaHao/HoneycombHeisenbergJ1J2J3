// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Haoxin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2022/5/12
*
* Description: GraceQ/XTRG project.
*/


/**
@file xtrg_params.h
@brief The struct XtrgParams
*/



#ifndef HONEYCOMBHEISENBERG_XTRG_COMMON_XTRG_PARAMS_H
#define HONEYCOMBHEISENBERG_XTRG_COMMON_XTRG_PARAMS_H

#include "qlmps/consts.h"

namespace xtrg {
using namespace qlmps;

const size_t kMaxTaylorExpansionOrder = 30;
const double kToleranceTaylorExpansionError = 1e-13;

const size_t kMaxVariationSweeps = 10;
const double kVariationConvergeTolerance = 1e-13;

const std::string kMpoPathPreFix = "mpo";


// MPO variation optimize params
struct MpoVOptimizeParams {
  MpoVOptimizeParams(
      const std::string& initial_mpo_path,
      const size_t Dmin, const size_t Dmax,
      const double trunc_err,
      const size_t sweeps = kMaxVariationSweeps,
      const double converge_tolerance = kVariationConvergeTolerance,
      const std::string& temp_path = kRuntimeTempPath
      ) :
      Dmin(Dmin), Dmax(Dmax), trunc_err(trunc_err),
      sweeps(sweeps), converge_tolerance(converge_tolerance),
      initial_mpo_path(initial_mpo_path),
      temp_path(temp_path) {}

  size_t Dmin;
  size_t Dmax;
  double trunc_err;

  size_t sweeps;
  double converge_tolerance;

  std::string initial_mpo_path;
  std::string temp_path;
};



struct XtrgParams {
  XtrgParams(
      const double tau,
      const size_t M,
      const size_t dmin, const size_t dmax, const double trunc_err,
      const size_t taylor_expansion_order = kMaxTaylorExpansionOrder,
      const double tolerace_taylor_expansion = kToleranceTaylorExpansionError,
      const size_t sweeps_variation = kMaxVariationSweeps,
      const double converge_tolerance_variation = kVariationConvergeTolerance,
      const std::string& mpo_path_prefix = kMpoPathPreFix,
      const std::string& temp_path = kRuntimeTempPath
  ) :
      tau(tau), M(M),
      Dmin(dmin), Dmax(dmax), trunc_err(trunc_err),
      taylor_expansion_order(taylor_expansion_order),
      tolerace_taylor_expansion(tolerace_taylor_expansion),
      sweeps_variation(sweeps_variation),
      converge_tolerance_variation(converge_tolerance_variation),
      mpo_path_prefix(mpo_path_prefix),
      temp_path(temp_path) {}

  // physical
  double tau;     // the initial value of reciprocal temperature
  size_t M;       // the times of double. when M = 0, only density matrix of beta = tau is found.

  // precision
  // global
  size_t Dmin;
  size_t Dmax;      // bond dimension for constructing density matrix and doubling
  double trunc_err; // truncation error globally

  // local
  size_t taylor_expansion_order; // maximal taylor expansion order
  double tolerace_taylor_expansion; // convergence tolerance of taylor expansion, i.e. the upper bond of remainder terms

  size_t sweeps_variation;  //maximal sweep time when variational optimize FiniteMPO
  double converge_tolerance_variation; // convergence tolerance when variational optimize FiniteMPO

  // Advanced parameters
  /// Runtime temporary files directory path
  std::string mpo_path_prefix;
  std::string temp_path;
};
} /* xtrg */





#endif //HONEYCOMBHEISENBERG_XTRG_COMMON_XTRG_PARAMS_H
