//
// Created by Hao-Xin on 2022/5/22.
//



#include <assert.h>
#include "XYMCExecutor.h"
#include <iostream>
#include <fstream>
#include "mpi.h"
#include "MC_params_case.h"

const size_t DimDof = 2; //XY model
const size_t NumOfCouplingType = 2; // only D, J1;
//const double xy_square_beta_c = 1.12;

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

  CaseParams params(argv[1]);
  size_t sweeps = params.Sweeps;
  std::string geometry = params.Geometry;
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;
  double beta = params.beta;
  if (argc > 2) {
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
  double J1Matrix[2][2];
  Vec2Mat(params.J1, J1Matrix);
  auto J1 = CouplingStructure<DimDof>(J1Matrix);
  auto D = CouplingStructure<DimDof>(ZeroMatrix);
  double hp = params.J2[0];

  std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures = {D, J1};

  XYPhysParams phys_params(
      geometry,
      Lx,
      Ly,
      beta,
      coupling_structures,
      hp
      );
  MCParams mc_params;
  mc_params.sweeps = sweeps;
  mc_params.sample_interval = params.SampleInterval;
  mc_params.print_interval = 100;
  mc_params.filename_postfix = "xy-rank" + std::to_string(world_rank)
      + params.Geometry
      + "J1" + std::to_string(J1Matrix[0][0])
      + "Hp" + std::to_string(hp)
      + "beta" + std::to_string(beta)
      + "L" + std::to_string(Lx);

  XYMCExecutor executor(mc_params, phys_params);
  executor.Execute();
  endTime = clock();
  std::cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  MPI_Finalize();
  return 0;
}