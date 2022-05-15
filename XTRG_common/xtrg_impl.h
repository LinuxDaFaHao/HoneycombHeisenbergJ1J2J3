// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Haoxin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2022/5/10
*
* Description: GraceQ/XTRG project.
*/

/**
@file xtrg_impl.h
@brief The entrance of the xtrg algorithm
*/

#ifndef HONEYCOMBHEISENBERG_XTRG_COMMON_XTRG_IMPL_H
#define HONEYCOMBHEISENBERG_XTRG_COMMON_XTRG_IMPL_H

#include "gqten/gqten.h"
#include "gqmps2/gqmps2.h"
#include "xtrg_params.h"

namespace xtrg {
using namespace gqten;
using namespace gqmps2;


std::string GenDensityMatrixMpoPath(
    const std::string& mpo_path_prefix,
    const double tau
);

void DumpThermodynamicQuantity(
    const std::vector<double> beta,
    const std::vector<double> Fn_set,
    const std::string &basename
);

/**
 * * The entrance of the xtrg algorithm.
 * Up to now, only bosonic operators can be measured.
 * Up to now, only one point functions and a group/type of correlations can be measured.
 *
 * The results of partition function, free energy, and correlations will be write to disk at finial.
 *
 * @tparam TenElemT
 * @tparam QNT
 * @param Hamiltonian
 * @param xtrg_params
 * @param density_matrix  output, must be presetted with length N
 */
template <typename TenElemT, typename QNT>
void Xtrg(
    const FiniteMPO<TenElemT, QNT>& Hamiltonian,
    const XtrgParams xtrg_params,
    FiniteMPO<TenElemT, QNT>& density_matrix
    ) {
  assert( Hamiltonian.size() == density_matrix.size() );
  double beta(xtrg_params.tau);
  double norm2 = InitializeDensityMatrix(Hamiltonian, xtrg_params, density_matrix);
  std::string mpo_rho_path = GenDensityMatrixMpoPath(xtrg_params.mpo_path_prefix, beta);
  density_matrix.Dump(mpo_rho_path);

  std::cout << "========== Doubling the density matrices ==========="
            << std::endl;

  const size_t M = xtrg_params.M;
  std::vector<double> beta_set(M+1), Fn_set(M+1);
  for(size_t step = 0; step < xtrg_params.M + 1; step++) {
    beta = beta*2.0;
    std::cout << "step = " << step
              << ",  beta = " << beta;
    if(step == 0) {
      Fn_set[step] = - std::log(norm2) / xtrg_params.tau;
    } else {
      Fn_set[step] = Fn_set[step - 1] - 2.0 * std::log(norm2) /beta;
    }
    beta_set[step] = beta;
    std::cout << ",\tF(beta = " << beta <<") =" << Fn_set[step] <<"\n" <<std::endl;
    if(step < M) {
      Timer double_rho_timer("double_density_matrix");
      mpo_rho_path = GenDensityMatrixMpoPath(xtrg_params.mpo_path_prefix, beta);
      MpoVOptimizeParams mpo_v_optimize_params(
          mpo_rho_path,
          xtrg_params.Dmin, xtrg_params.Dmax,
          xtrg_params.trunc_err, xtrg_params.sweeps_variation,
          xtrg_params.converge_tolerance_variation,
          xtrg_params.temp_path
      );
      norm2 =  density_matrix.SquareAndNormlize(mpo_v_optimize_params);
      density_matrix.Dump(mpo_rho_path);
      double_rho_timer.PrintElapsed();
    }
  }
  DumpThermodynamicQuantity(beta_set,  Fn_set, "free_energy");
  return;
}


std::string GenDensityMatrixMpoPath(
    const std::string& mpo_path_prefix,
    const double tau
) {
  return mpo_path_prefix + "_rho"  + "_tau" + std::to_string(tau);
}



void DumpThermodynamicQuantity(
    const std::vector<double> beta,
    const std::vector<double> Fn_set,
    const std::string &basename
) {
  auto file = basename + ".json";
  std::ofstream ofs(file);

  ofs << "[\n";

  for(size_t i = 0; i < beta.size(); i++) {
    ofs << "  [";

    DumpAvgVal(ofs, beta[i]);
    ofs << ", ";
    DumpAvgVal(ofs, Fn_set[i]);
    if (i == beta.size()-1) {
      ofs << "]\n";
    } else {
      ofs << "],\n";
    }
  }
  ofs << "]";
  ofs.close();
}





} //xtrg









#endif //HONEYCOMBHEISENBERG_XTRG_COMMON_XTRG_IMPL_H
