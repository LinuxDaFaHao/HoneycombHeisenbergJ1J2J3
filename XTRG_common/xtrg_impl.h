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


void DumpThermodynamicQuantity(
    const std::vector<double> beta,
    const std::vector<double> Z_set,
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


  MpoVOptimizeParams mpo_v_optimize_params(
      xtrg_params.Dmin, xtrg_params.Dmax,
      xtrg_params.trunc_err, xtrg_params.sweeps_variation,
      xtrg_params.converge_tolerance_variation,
      xtrg_params.temp_path
      );
  std::cout << "========== Doubling the density matrices ==========="
            << std::endl;

  const size_t M = xtrg_params.M;
  std::vector<double> beta_set(M+1), Z_set(M+1), Fn_set(M+1);
  for(size_t step = 0; step < xtrg_params.M; step++) {
    beta = beta*2.0;
    std::cout << "step = " << step
              << ",  beta = " << beta;
    double Z = norm2 * norm2; // partition function of square of current density matrix;
    double F = -std::log(Z) / (beta);  // partition function of 2*beta
    std::cout << "F(beta = " << beta <<") =" << F <<"\n" <<std::endl;
    beta_set[step] = beta;
    Z_set[step] = Z;
    Fn_set[step] = F;

    Timer double_rho_timer("double_density_matrix");
    norm2 = Z * density_matrix.SquareAndNormlize(mpo_v_optimize_params);
    double_rho_timer.PrintElapsed();
  }
  beta_set.back() = beta*2.0;
  std::cout << "step = end"
            << ",  beta = " << 2 * beta;
  Z_set.back() = norm2 * norm2;
  Fn_set.back() = -std::log(Z_set.back()) / (beta_set.back());
  std::cout << "F(beta = " << beta_set.back() <<") =" << Fn_set.back() <<"\n" <<std::endl;
  DumpThermodynamicQuantity(beta_set, Z_set, Fn_set, "free_energy");
  return;
}


void DumpThermodynamicQuantity(
    const std::vector<double> beta,
    const std::vector<double> Z_set,
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
    DumpAvgVal(ofs, Z_set[i]);
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
