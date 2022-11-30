/*
* Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2022/10/3.
*
* Description: Simulation for the clock model with random field which doesn't break time reversal symmetry
*
*/


#include <assert.h>
#include "ClockMCExecutor.h"
#include <iostream>
#include "mpi.h"
#include "clock_params_case.h"


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  clock_t startTime, endTime;
  startTime = clock();
  int world_size(1), world_rank(0);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  random_engine.seed((int) startTime + world_rank);

  CaseParams params(argv[1]);
  size_t sweeps = params.Sweeps;
  std::string geometry = params.Geometry;
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;
  double beta = params.beta;


  if (argc > 2) {
    beta = std::atof(argv[2]);
  }

  double J_matrix[1][1] = {{params.J}};
  double Zero_matrix[1][1] = {{0}};
  auto coupling_J = CouplingStructure<clock_dim_dof>(J_matrix);
  auto D = CouplingStructure<clock_dim_dof>(Zero_matrix);

  MCParams mc_params;
  mc_params.sweeps = sweeps;
  mc_params.cluster_radius = params.ClusterRadius;
  mc_params.sample_interval = params.SampleInterval;
  mc_params.print_interval = 100;
  mc_params.filename_postfix = "clock-rank" + std::to_string(world_rank)
      + params.Geometry
      + "J" + std::to_string(params.J)
      + "hp" + std::to_string(params.hp)
      + "beta" + std::to_string(beta)
      + "L" + std::to_string(Lx);



  std::array<CouplingStructure<clock_dim_dof>, clock_interaction_range> coupling_structures = {D, coupling_J};
  ClockPhysParams phys_params(
      geometry,
      Lx, Ly,
      beta,
      coupling_structures,
      params.hp
      );

  ClockMCExecutor<6> executor(mc_params, phys_params, world_rank);
  executor.Execute();




  endTime = clock();
  std::cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  MPI_Finalize();
  return 0;
}