// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Haoxin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2022/5/10
*
* Description: GraceQ/XTRG project.
*/

/**
@file density_matrix_initializer.h
@brief initialize the density matrix with beta=tau (a small value) by Taylor expansion
*/


#ifndef HONEYCOMBHEISENBERG_XTRG_COMMON_DENSITY_MATRIX_INITIALIZER_H
#define HONEYCOMBHEISENBERG_XTRG_COMMON_DENSITY_MATRIX_INITIALIZER_H



#include "qlten/qlten.h"
#include "qlmps/qlmps.h"
#include "./mpopp.h"
#include "./mpo_utility.h"
#include "./xtrg_params.h"
namespace xtrg {
using namespace qlten;
using namespace qlmps;


//forward declaration

std::string GenHnMpoPath(
    const std::string& mpo_path_prefix,
    const std::size_t n, //the power
    const double tau
);
/* FYI, lanczos figure
 * |----0                       0-----
 * |          2         2            |
 * |          |         |            |
 * |----1 0-------3 0--------3  1-----
 * |          |        |             |
 * |          1       1 2            |
 * |          |        |             |
 * |----2 0-------------------3 2----|
 */
/**
 * Find the MPO of density matrix with beta=params.tau by Taylor expansion.
 * The output density_matrix is a normlized density matrix (by 2 norm, not the usual trace 1).
 *
 *
 * @tparam TenElemT
 * @tparam QNT
 * @param hamiltonian
 * @param params
 * @param density_matrix
 * @return the normalization factor
 */
template <typename TenElemT, typename QNT>
double InitializeDensityMatrix(
    const FiniteMPO<TenElemT, QNT>& hamiltonian,
    const XtrgParams params,
    FiniteMPO<TenElemT, QNT>& density_matrix //output
    ) {
  assert( hamiltonian.size() == density_matrix.size() );
  const QLTEN_Double tau = params.tau;
  const QLTEN_Double trunc_err = params.trunc_err;
  const size_t Dmin = params.Dmin;
  const size_t Dmax = params.Dmax;
  const size_t sweep_time_max = params.sweeps_variation;
  const std::string temp_path = params.temp_path;
  const size_t max_taylor_expansion_order = params.taylor_expansion_order;
  const QLTEN_Double tolerance_taylor_expansion_error = params.tolerace_taylor_expansion;

  std::cout << "Constructing the density matrix with tau = " << tau << " by Taylor expansion" << std::endl;
  Timer taylor_expansion_timer("taylor_expansion");
//  using Tensor = QLTensor<TenElemT, QNT>;
  using MPOT = FiniteMPO<TenElemT, QNT>;
  const size_t N = hamiltonian.size();
  const size_t bond_dimension_of_H = hamiltonian.GetMaxBondDimension();

  std::cout << "Maximal bond dimension of hamiltonian = " << bond_dimension_of_H << std::endl;
  MPOT m_tau_h(hamiltonian); // minus tau times hamiltonian
  m_tau_h.Scale(-tau);
  GenerateIndentiyMPO(hamiltonian, density_matrix);
  TenElemT t0 = density_matrix.Trace();
  std::cout << "power i = 0, Identity, Trace(0) = " << std::setprecision(7) << std::scientific << t0 << std::endl;
  MPOT A(m_tau_h);  //temp mpo's
  // A =  (-tau * hamiltonian) ^ i / ( i! )
  // temp_multiply_mpo -> A


  if (!IsPathExist(params.temp_path)) {
    CreatPath(params.temp_path);
  }
  double norm2(1.0);
  for(size_t i = 1; i < max_taylor_expansion_order; i++){
    std::cout << "power i = " << i << std::endl;
    Timer power_n_timer("power_n");
    MPOT temp_multiply_mpo(N);
    std::string mpo_Hn_path = xtrg::GenHnMpoPath(params.mpo_path_prefix, i, tau);
    if(i==1) {
      //suppose bond_dimension_of_H <= Dmax.
    } else if( std::pow(bond_dimension_of_H, i) <= Dmax) {
      MpoProduct(A, m_tau_h, temp_multiply_mpo);
    } else { //need compress
      temp_multiply_mpo.Load(mpo_Hn_path);

      MpoProduct(A, m_tau_h, temp_multiply_mpo,
                 Dmin, Dmax, trunc_err,
                 sweep_time_max, params.converge_tolerance_variation,
                 temp_path);
    }
    if(i > 1){
      temp_multiply_mpo.Scale((1.0)/i);
      // TODO: write a health operator=(&&)
      // The natural operator=(&&) will std::move(rhs.ten_cano_type_) and release
      // its memory. We don't want it. We just want set it to NONE;
      A = std::move(temp_multiply_mpo);
    }
    density_matrix += A;
    A.Dump(mpo_Hn_path);

    TenElemT t = A.Trace();
    size_t bond_dimension_of_A = A.GetMaxBondDimension();
    if( density_matrix.GetMaxBondDimension() > params.Dmax ) {
      norm2 *= density_matrix.Truncate(1e-16, params.Dmax, params.Dmax);
    }
    double power_n_elapsed_time = power_n_timer.Elapsed();
    std::cout << "power i = " << i ;
    std::cout << " D(H^"<< i <<") = " << std::setw(5) << bond_dimension_of_A
              << " Time = " << std::setw(8) << power_n_elapsed_time
              << " Trace("<<i<<") = " << std::setprecision(7) << std::scientific << t;
    std::cout << std::scientific << std::endl;
    if( i%2 == 0) {
      assert(t > 0.0);
      if( std::abs(t / t0) < tolerance_taylor_expansion_error) {
        break;
      }
    }
  }

  density_matrix.Centralize(0);
  norm2 *= density_matrix.Normlize();
  std::cout << "The series expansion was finished. The bond dimension of rho(tau) = "
            << density_matrix.GetMaxBondDimension() << std::endl;

  taylor_expansion_timer.PrintElapsed();
  return norm2;
}

std::string GenHnMpoPath(
    const std::string& mpo_path_prefix,
    const std::size_t n, //the power
    const double tau
) {
  return mpo_path_prefix + "_H" + std::to_string(n) + "_tau" + std::to_string(tau);
}

} //xtrg





#endif //HONEYCOMBHEISENBERG_XTRG_COMMON_DENSITY_MATRIX_INITIALIZER_H



