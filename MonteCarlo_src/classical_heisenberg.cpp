/*
* Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2022/5/22.
*
* Description: Simulation for the classical heisenberg model.
*
*/

/**
@file classical_heisenberg.cpp
@brief Simulation for the classical heisenberg model,
       support the square/honeycomb lattice, with D(single ion anistropy), J1, J2, J3 interactions.
       Hamiltonian: H = -\sum_i S_i*D*S_i - \sum_<i,j> S_i*J1*S_j - J2 terms - J3 terms
@usage mpirun -n ${SLURM_NTASKS} ./heisenberg params.json
*/

#include <assert.h>
#include "WolfMCExecutor.h"
#include <iostream>
#include "mpi.h"
#include "MC_params_case.h"

const size_t DimDof = 3; //heisenberg model
const size_t NumOfCouplingType = 4; // D, J1, J2, J3

void Vec2Mat(
    std::vector<double> j_vec,
    double j_mat[3][3]
    ) {
  for(size_t i = 0; i < 3; i++ ) {
    for(size_t j = 0; j < 3; j++ ) {
      j_mat[i][j] = j_vec[3*i + j];
    }
  }
}


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  clock_t startTime,endTime;
  startTime = clock();
  int world_size, world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  random_engine.seed((int) startTime + world_rank);

  CaseParams params(argv[1]);
  size_t sweeps = params.Sweeps;
  std::string geometry = params.Geometry;
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;
  double beta = params.beta;

  if(argc > 2) {
    beta = std::atof(argv[2]);
  }

  double IdentityMatrix[DimDof][DimDof], ZeroMatrix[DimDof][DimDof];
  for(size_t i = 0; i < DimDof; i++){
    for(size_t j = 0; j < DimDof; j++) {
      if(i == j) {
        IdentityMatrix[i][j] = 1;
      } else{
        IdentityMatrix[i][j] = 0;
      }
      ZeroMatrix[i][j] = 0;
    }
  }
  double J1Matrix[3][3], J2Matrix[3][3], J3Matrix[3][3], DMatrix[3][3];
  Vec2Mat(params.J1, J1Matrix);
  auto J1 = CouplingStructure<DimDof>(J1Matrix);
  Vec2Mat(params.J2, J2Matrix);
  auto J2 = CouplingStructure<DimDof>(J2Matrix);
  Vec2Mat(params.J3, J3Matrix);
  auto J3 = CouplingStructure<DimDof>(J3Matrix);
  Vec2Mat(params.D, DMatrix);
  auto D = CouplingStructure<DimDof>(DMatrix);
  std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures = {D, J1,J2,J3};

  PhysParams phys_params(
      geometry,
      Lx,
      Ly,
      beta,
      coupling_structures
      );
  MCParams mc_params;
  mc_params.sweeps = sweeps;
  mc_params.print_interval = 100;
  mc_params.filename_postfix = "hei-rank" + std::to_string(world_rank)
      + params.Geometry
      + "J1zz" + std::to_string(J1Matrix[2][2])
      + "J2zz" + std::to_string(J2Matrix[2][2])
      + "J3zz" + std::to_string(J3Matrix[2][2])
      + "Dzz" + std::to_string(DMatrix[2][2])
      + "beta" + std::to_string(beta)
      + "L" + std::to_string(Lx);

  WolfMCExecutor executor(mc_params, phys_params);
  executor.Execute();
  endTime = clock();
  std::cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  MPI_Finalize();
  return 0;
}