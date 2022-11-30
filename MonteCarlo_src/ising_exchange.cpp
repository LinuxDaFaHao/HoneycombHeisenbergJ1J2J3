/**
 * Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
 * Created Data: 2022/10/16.
 *
 * Description: Calculation for Ising model via Exchange Monte Carlo + cluster update.
 */

#include <assert.h>
#include "ExchangeMCExecutor.h"
#include <iostream>
#include "MC_params_case.h"

const size_t DimDof = 1; //ising model
const size_t NumOfCouplingType = 4; // D, J1, J2, J3;
const double ising_square_beta_c = 0.4406867935097714683578828953614;
const double ising_honeycomb_beta_c = 0.65847894846240828670147493539844;

int main(int argc, char **argv) {
  clock_t startTime, endTime;
  startTime = clock();
  random_engine.seed((int) startTime + 35);

  CaseParams params(argv[1]);
  size_t sweeps = params.Sweeps;
  std::string geometry = params.Geometry;
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;
  if (argc < 5) {
    std::cout << "lack arguments" << std::endl;
    exit(1);
  }
  double beta1 = std::atof(argv[2]);
  double beta2 = std::atof(argv[3]);
  size_t thread = std::atoi(argv[4]);
  double beta_max = std::max(beta1, beta2);
  double beta_min = std::min(beta1, beta2);
  std::cout << "Read argument beta_max = " << beta_max << std::endl;
  std::cout << "              beta_min = " << beta_min << std::endl;
  std::cout << "              thread_num = " << thread << std::endl;

  double IdentityMatrix[DimDof][DimDof], ZeroMatrix[DimDof][DimDof];
  for (size_t i = 0; i < DimDof; i++) {
    for (size_t j = 0; j < DimDof; j++) {
      if (i == j) {
        IdentityMatrix[i][j] = 1;
      } else {
        IdentityMatrix[i][j] = 0;
      }
      ZeroMatrix[i][j] = 0;
    }
  }
  auto I = CouplingStructure<DimDof>(IdentityMatrix);
  auto J1 = I;
  J1 *= params.J1[0];
  auto J2 = I;
  J2 *= params.J2[0];
  auto J3 = I;
  J3 *= params.J3[0];
  auto D = I;
  D *= params.D[0];
  std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures = {D, J1, J2, J3};

  PhysParams phys_params(
      geometry,
      Lx,
      Ly,
      beta_min,
      beta_max,
      coupling_structures,
      LocalDOF<DimDof>({1.0}),
      LocalDOF<DimDof>({1.0}),
      LocalDOF<DimDof>({1.0})
  );
  ExchangeMCParams mc_params;
  mc_params.thread_num = thread;
  mc_params.adjust_temperature_times = 10;
  mc_params.adjust_temperature_samples = sweeps / 100;
  mc_params.sweeps = sweeps;
  mc_params.sample_interval = params.SampleInterval;
  mc_params.warmup_sample_num = params.Warmup_Samples;
  mc_params.exchange_interval = params.ExchangeInterval;
  mc_params.print_interval = 100;
  mc_params.filename_postfix = "ising" + params.Geometry
      + "J1" + std::to_string(params.J1[0])
      + "J2" + std::to_string(params.J2[0])
      + "J3" + std::to_string(params.J3[0])
      + "beta_max" + std::to_string(beta_max)
      + "beta_min" + std::to_string(beta_min)
      + "L" + std::to_string(Lx);

  ExchangeMCExecutor executor(mc_params, phys_params);
  executor.Execute();
  endTime = clock();
  std::cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  return 0;
}