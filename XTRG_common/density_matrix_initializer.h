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



#include "gqten/gqten.h"
#include "gqmps2/gqmps2.h"
#include "./mpopp.h"
#include "./mpo_utility.h"
#include "./xtrg_params.h"
namespace xtrg {
using namespace gqten;
using namespace gqmps2;


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
  const GQTEN_Double tau = params.tau;
  const GQTEN_Double trunc_err = params.trunc_err;
  const size_t Dmin = params.Dmin;
  const size_t Dmax = params.Dmax;
  const size_t sweep_time_max = params.sweeps_variation;
  const std::string temp_path = params.temp_path;
  const size_t max_taylor_expansion_order = params.taylor_expansion_order;
  const GQTEN_Double tolerance_taylor_expansion_error = params.tolerace_taylor_expansion;

  std::cout << "Constructing the density matrix with tau = " << tau << " by Taylor expansion" << std::endl;
  Timer taylor_expansion_timer("taylor_expansion");
//  using Tensor = GQTensor<TenElemT, QNT>;
  using MPOT = FiniteMPO<TenElemT, QNT>;
  const size_t N = hamiltonian.size();
  const size_t bond_dimension_of_H = hamiltonian.GetMaxBondDimension();

  std::cout << "Maximal bond dimension of hamiltonian = " << bond_dimension_of_H << std::endl;
  MPOT m_tau_h(hamiltonian); // minus tau times hamiltonian
  m_tau_h.Scale(-tau);
  GenerateIndentiyMPO(hamiltonian, density_matrix);
  MPOT temp_sum_mpo(N), A(m_tau_h), temp_multiply_mpo(N);  //temp mpo's
  // A =  (-tau * hamiltonian) ^ i / ( i! )
  // temp_multiply_mpo -> A


  if (!IsPathExist(params.temp_path)) {
    CreatPath(params.temp_path);
  }

  for(size_t i = 1; i < max_taylor_expansion_order; i++){
    std::cout << "power i = " << i ;
    Timer power_n_timer("power_n");
    if(i==1) {
      //suppose bond_dimension_of_H <= Dmax.
    } else if( std::pow(bond_dimension_of_H, i) <= Dmax) {
      MpoProduct(A, m_tau_h, temp_multiply_mpo);
    } else { //need compress
      MpoProduct(A, m_tau_h, temp_multiply_mpo,
                 Dmin, Dmax, trunc_err,
                 sweep_time_max, params.converge_tolerance_variation,
                 temp_path);
      std::cout << "power i = " << i ;
    }
    if(i > 1){
      temp_multiply_mpo.Scale((1.0)/i);
      // TODO: write test, debug and fix operator= (&, &&)
//      A = std::move(temp_multiply_mpo);
      for(size_t i = 0; i < N; i++) {
        delete A(i);
        A(i) = temp_multiply_mpo(i);
        temp_multiply_mpo(i) = nullptr;
      }
    }
    density_matrix += A;

    TenElemT t = A.Trace();
    size_t bond_dimension_of_A = A.GetMaxBondDimension();
    size_t power_n_elapsed_time = power_n_timer.Elapsed();
    std::cout << " D(H^"<< i <<") = " << std::setw(5) << bond_dimension_of_A
              << " Time = " << std::setw(8) << power_n_elapsed_time
              << " Trace("<<i<<") = " << std::setprecision(7) << std::scientific << t;
    std::cout << std::scientific << std::endl;
    if( i%2 == 0) {
      assert(t > 0.0);
      if( std::abs(t) < tolerance_taylor_expansion_error) {
        break;
      }
    }
  }

  density_matrix.Centralize(0);
  double norm2 = density_matrix.Normlize();
  std::cout << "The series expansion was finished. The bond dimension of rho(tau) = "
            << density_matrix.GetMaxBondDimension() << std::endl;
//TODO: test the results of here if need
  taylor_expansion_timer.PrintElapsed();
  return norm2;
}





} //xtrg





#endif //HONEYCOMBHEISENBERG_XTRG_COMMON_DENSITY_MATRIX_INITIALIZER_H



