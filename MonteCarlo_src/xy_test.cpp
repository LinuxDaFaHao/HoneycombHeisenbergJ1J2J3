//
// Created by Hao-Xin on 2022/5/22.
//



#include <assert.h>
#include "WolfMCExecutor.h"
#include <iostream>
#include <fstream>
#include "mpi.h"
#include "MC_params_case.h"

const size_t DimDof = 2; //XY model
const size_t NumOfCouplingType = 2; // only D, J1;
const double xy_square_beta_c = 1.12;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  clock_t startTime,endTime;
  startTime = clock();
  int world_size, world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  random_engine.seed((int) startTime + world_rank);
//  random_engine.seed(0);

  CaseParams params(argv[1]);
  size_t sweeps = params.Sweeps;
  std::string geometry = "Square";
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;
  double beta = params.beta;
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
  double Matrix1[3][3] = {{0,0,0},{0,0,0},{0,0,1}};
  double Matrixz[3][3] = {{0,0,0},{0,0,0},{0,0,1}};
  auto J1 = CouplingStructure<DimDof>(IdentityMatrix);
  auto D = CouplingStructure<DimDof>(ZeroMatrix);
  auto J3 = D;
  auto J2 = CouplingStructure<DimDof>(ZeroMatrix);
  std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures = {D, J1};

  PhysParams phys_params(
      geometry,
      Lx,
      Ly,
      beta,
      coupling_structures
      );
  MCParams mc_params;
  mc_params.sweeps = sweeps;
  mc_params.sample_interval = 10;
  mc_params.print_interval = 100;
  mc_params.filename_postfix = "xy" + std::to_string(world_rank)
      + "beta" + std::to_string(beta)
      + "L" + std::to_string(Lx);

  WolfMCExecutor executor(mc_params, phys_params);
  executor.Execute();
  endTime = clock();
  std::cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  MPI_Finalize();
  return 0;
}