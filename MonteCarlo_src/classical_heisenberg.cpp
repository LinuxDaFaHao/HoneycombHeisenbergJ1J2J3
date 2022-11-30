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


void Vec2Mat(
    std::vector<double> j_vec,
    double j_mat[3][3]
) {
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      j_mat[i][j] = j_vec[3 * i + j];
    }
  }
}

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

  double J1Matrix[3][3], J2Matrix[3][3], J3Matrix[3][3], DMatrix[3][3];
  double RMatrix[3][3] = {{cos(2 * M_PI / 3), -sin(2 * M_PI / 3), 0.0},
                          {sin(2 * M_PI / 3), cos(2 * M_PI / 3), 0.0},
                          {0.0, 0.0, 1.0}};
  double RinvMatrix[3][3] = {{cos(-2 * M_PI / 3), -sin(-2 * M_PI / 3), 0.0},
                             {sin(-2 * M_PI / 3), cos(-2 * M_PI / 3), 0.0},
                             {0.0, 0.0, 1.0}};
  auto R = CouplingStructure<DimDof>(RMatrix);
  auto Rinv = CouplingStructure<DimDof>(RinvMatrix);

  Vec2Mat(params.J1, J1Matrix);
  auto J11 = CouplingStructure<DimDof>(J1Matrix);

  Vec2Mat(params.J2, J2Matrix);
  auto J21 = CouplingStructure<DimDof>(J2Matrix);

  Vec2Mat(params.J3, J3Matrix);
  auto J31 = CouplingStructure<DimDof>(J3Matrix);

  bool Interaction_SU2 = true;
  Interaction_SU2 = Interaction_SU2 && J11.IsIsometry();
  Interaction_SU2 = Interaction_SU2 && J21.IsIsometry();
  Interaction_SU2 = Interaction_SU2 && J31.IsIsometry();

  Vec2Mat(params.D, DMatrix);
  auto D = CouplingStructure<DimDof>(DMatrix);

  MCParams mc_params;
  mc_params.sweeps = sweeps;
  mc_params.sample_interval = params.SampleInterval;
  mc_params.print_interval = 100;
  mc_params.filename_postfix = "hei-rank" + std::to_string(world_rank)
      + params.Geometry
      + "J1zz" + std::to_string(J1Matrix[2][2])
      + "J2zz" + std::to_string(J2Matrix[2][2])
      + "J3zz" + std::to_string(J3Matrix[2][2])
      + "Dzz" + std::to_string(DMatrix[2][2])
      + "beta" + std::to_string(beta)
      + "L" + std::to_string(Lx);

  if (geometry == "Honeycomb" && Interaction_SU2) {
    const size_t NumOfCouplingType = 4; // D, J1, J2, J3
    std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures = {D, J11, J21, J31};
    PhysParams phys_params(
        geometry,
        Lx,
        Ly,
        beta,
        coupling_structures
    );
    WolfMCExecutor executor(mc_params, phys_params, world_rank);
    executor.Execute();
  } else if(geometry == "Honeycomb3" && Interaction_SU2) {
    const size_t NumOfCouplingType = 10;
    std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures = {D, J11, J11, J11,
                                                                                    J21, J21, J21,
                                                                                    J31, J31, J31};
    PhysParams phys_params(
        geometry,
        Lx,
        Ly,
        beta,
        coupling_structures
    );
    WolfMCExecutor executor(mc_params, phys_params, world_rank);
    executor.Execute();
  } else if(geometry == "Honeycomb3" && J21.IsIsometry()){
    const size_t NumOfCouplingType = 10;
    CouplingStructure<DimDof> J32 = R * (J31 * Rinv);
    CouplingStructure<DimDof> J33 = Rinv * (J31 * R);
    CouplingStructure<DimDof> J12 = Rinv * J11 * R;
    CouplingStructure<DimDof> J13 = R * J11 * Rinv;

    std::array<CouplingStructure<DimDof>, NumOfCouplingType> coupling_structures = {D, J11, J12, J13,
                                                                                    J21, J21, J21,
                                                                                    J31, J32, J33};
    std::array<double, DimDof> gs1_data, gs2_data, gs3_data;
    for(size_t i =0; i < DimDof; i++) {
      gs1_data[i] = params.GS1[i];
      gs2_data[i] = params.GS2[i];
      gs3_data[i] = params.GS3[i];
    }

    PhysParams phys_params(
        geometry,
        Lx,
        Ly,
        beta,
        coupling_structures,
        LocalDOF<DimDof>(gs1_data),
        LocalDOF<DimDof>(gs2_data),
        LocalDOF<DimDof>(gs3_data)
    );
    WolfMCExecutor executor(mc_params, phys_params, world_rank);
    executor.Execute();
  } else {
    std::cout << "Unexpected parameter case." << std::endl;
    exit(1);
  }

  endTime = clock();
  std::cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  MPI_Finalize();
  return 0;
}