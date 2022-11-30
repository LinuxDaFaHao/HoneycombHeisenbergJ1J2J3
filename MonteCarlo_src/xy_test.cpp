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
const size_t NumOfCouplingType = 4; // only D, J1;
const double xy_square_beta_c = 1.12;

void Vec2Mat(
    std::vector<double> j_vec,
    double j_mat[DimDof][DimDof]
) {
  for(size_t i = 0; i < DimDof; i++ ) {
    for(size_t j = 0; j < DimDof; j++ ) {
      j_mat[i][j] = j_vec[DimDof*i + j];
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
//  random_engine.seed(0);

  CaseParams params(argv[1]);
  size_t sweeps = params.Sweeps;
  std::string geometry = params.Geometry;
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
//  double J1Matrix[2][2], J2Matrix[2][2], J3Matrix[2][2], DMatrix[2][2];
//  Vec2Mat(params.J1, J1Matrix);
//  auto J11 = CouplingStructure<DimDof>(J1Matrix);
//  Vec2Mat(params.J12, J1Matrix);
//  auto J12 = CouplingStructure<DimDof>(J1Matrix);
//  Vec2Mat(params.J13, J1Matrix);
//  auto J13 = CouplingStructure<DimDof>(J1Matrix);
//
//  Vec2Mat(params.J2, J2Matrix);
//  auto J21 = CouplingStructure<DimDof>(J2Matrix);
//  Vec2Mat(params.J22, J2Matrix);
//  auto J22 = CouplingStructure<DimDof>(J2Matrix);
//  Vec2Mat(params.J23, J2Matrix);
//  auto J23 = CouplingStructure<DimDof>(J2Matrix);
//
//  Vec2Mat(params.J3, J3Matrix);
//  auto J31 = CouplingStructure<DimDof>(J3Matrix);
//  Vec2Mat(params.J32, J3Matrix);
//  auto J32 = CouplingStructure<DimDof>(J3Matrix);
//  Vec2Mat(params.J33, J3Matrix);
//  auto J33 = CouplingStructure<DimDof>(J3Matrix);

  auto J1 = CouplingStructure<DimDof>(ZeroMatrix);
  auto D = CouplingStructure<DimDof>(ZeroMatrix);

  auto J2 = CouplingStructure<DimDof>(ZeroMatrix);
  auto J3 = CouplingStructure<DimDof>(IdentityMatrix);
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