//
// Created by Hao-Xin on 2022/5/22.
//



#include <assert.h>
#include "WolfMCExecutor.h"
#include <iostream>
#include <fstream>
#include "mpi.h"
#include "MC_params_case.h"

const size_t DimDof = 1; //ising model
const size_t NumOfCouplingType = 4; // D, J1, J2, J3;
const double ising_square_beta_c = 0.4406867935097714683578828953614;
const double ising_honeycomb_beta_c = 0.65847894846240828670147493539844;

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
  std::string geometry = params.Geometry;
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;
  double beta = params.beta;
  if(beta < 0.0) {
    if(geometry == "Square"){
      beta = ising_square_beta_c;
    } else if (geometry == "Honeycomb") {
      beta = ising_honeycomb_beta_c;
    }
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
  auto I = CouplingStructure<DimDof>(IdentityMatrix);
  auto J1 = I; J1 *= params.J1[0];
  auto J2 = I; J2 *= params.J2[0];
  auto J3 = I; J3 *= params.J3[0];
  auto D = I; D *= params.D[0];
  std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures = {D, J1, J2, J3};

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
  mc_params.filename_postfix = "ising-rank" + std::to_string(world_rank)
      + params.Geometry
      + "J1" + std::to_string(params.J1[0])
      + "J2" + std::to_string(params.J2[0])
      + "J3" + std::to_string(params.J3[0])
      + "D" + std::to_string(params.D[0])
      + "beta" + std::to_string(beta)
      + "L" + std::to_string(Lx);

  WolfMCExecutor executor(mc_params, phys_params);
  executor.Execute();
  endTime = clock();
  std::cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  MPI_Finalize();
  return 0;
}