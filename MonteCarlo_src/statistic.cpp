/*
* Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2022/8/22.
*
* Description: Honeycomb lattice
*
*/


#include <assert.h>
#include "WolfMCExecutor.h"
#include <iostream>
#include "mpi.h"
#include "MC_params_case.h"

const size_t DimDof = 3; //heisenberg model

template<typename DataType>
std::vector<DataType> LoadData(
    const std::string &filename,
    size_t &data_size
) {
  std::ifstream ifs(filename, std::ifstream::binary);
  if (!ifs.good()) {
    std::cout << "open file " << filename << " failed." << std::endl;
    exit(1);
  }
  std::vector<DataType> data(data_size);
  for (size_t i = 0; i < data_size; i++) {
#ifndef NDEBUG
    if (ifs.eof()) {
      std::cout << "warning: data number in file " << filename << "is only " << i + 1 << std::endl;
      data_size = i + 1;
      break;
    }
#endif
    ifs >> data[i];
  }
  ifs.close();
  return data;
}

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
  int world_size, world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  random_engine.seed((int) startTime + world_rank);

  CaseParams params(argv[1]);
  size_t sweeps = params.Sweeps;
  std::string geometry = params.Geometry;
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;
  size_t N;
  if (geometry == "Honeycomb" || geometry == "Honeycomb3") {
    N = 2 * Lx * Ly;
  } else if(geometry == "Square") {
    N = Lx * Ly;
  } else {
    std::cout << " unknown geometry " << geometry << std::endl;
    exit(1);
  }
  double beta = params.beta;

  if (argc > 2) {
    beta = std::atof(argv[2]);
  }

  double J1Matrix[3][3], J2Matrix[3][3], J3Matrix[3][3], DMatrix[3][3];
  Vec2Mat(params.J1, J1Matrix);
  Vec2Mat(params.J2, J2Matrix);
  Vec2Mat(params.J3, J3Matrix);
  Vec2Mat(params.D, DMatrix);


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

  gqten::Timer load_data_timer("load_data");
  std::vector<double> res;
  auto &res_ = res;
  size_t data_size = mc_params.sweeps;
  auto energy_data = LoadData<double>("energy" + mc_params.filename_postfix, data_size);
//  auto stiffness_data = LoadData<double>("stiffness" + mc_params.filename_postfix, data_size);
  std::array<std::vector<double>, DimDof> sum_spin_data;
  for (size_t i = 0; i < DimDof; i++) {
    sum_spin_data[i] = LoadData<double>("sum_spin" + std::to_string(i) + mc_params.filename_postfix, data_size);
  }
  std::array<std::vector<double>, DimDof> zigzag_af_magnetization_data1, zigzag_af_magnetization_data2,
      zigzag_af_magnetization_data3;
  for (size_t i = 0; i < DimDof; i++) {
    zigzag_af_magnetization_data1[i] =
        LoadData<double>("magnetization" + std::to_string(i) + mc_params.filename_postfix, data_size);
    zigzag_af_magnetization_data2[i] =
        LoadData<double>("magnetization2" + std::to_string(i) + mc_params.filename_postfix, data_size);
    zigzag_af_magnetization_data3[i] =
        LoadData<double>("magnetization3" + std::to_string(i) + mc_params.filename_postfix, data_size);
  }

  std::vector<double> complex_xy_order_parameter_real_, complex_xy_order_parameter_imag_;
  complex_xy_order_parameter_real_ =
      LoadData<double>("xy_order_parameter_real" + mc_params.filename_postfix, data_size);
  complex_xy_order_parameter_imag_ =
      LoadData<double>("xy_order_parameter_imag" + mc_params.filename_postfix, data_size);

  load_data_timer.PrintElapsed();

  gqten::Timer statistic_timer("statistic");

  auto half_data_of_energy = std::vector(energy_data.begin() + data_size / 2, energy_data.begin() + data_size);
  res.push_back(Mean(half_data_of_energy) / N);
  // specific heat
  res.push_back(
      Variance(half_data_of_energy, res[0] * N) * beta * beta / N);
  // magnetic susceptibility
//  for (size_t i = 0; i < DimDof; i++) {
//    auto half_data_of_sum_spin =
//        std::vector(sum_spin_data[i].begin() + data_size / 2, sum_spin_data[i].begin() + data_size);
//    res.push_back(Variance(half_data_of_sum_spin) / phys_params.N);
//  }


// anti-ferromagnetic susceptibility
  for (size_t i = 0; i < DimDof; i++) {
    std::vector<double> valid_data_of_magnetization_abs1(mc_params.sweeps / 2, 0.0),
        valid_data_of_magnetization_abs2(mc_params.sweeps / 2, 0.0),
        valid_data_of_magnetization_abs3(mc_params.sweeps / 2, 0.0);
    for (size_t j = 0; j < valid_data_of_magnetization_abs1.size(); j++) {
      valid_data_of_magnetization_abs1[j] = abs(zigzag_af_magnetization_data1[i][j + sweeps / 2]);
      valid_data_of_magnetization_abs2[j] = abs(zigzag_af_magnetization_data2[i][j + sweeps / 2]);
      valid_data_of_magnetization_abs3[j] = abs(zigzag_af_magnetization_data3[i][j + sweeps / 2]);
    }
    res_.push_back(Variance((valid_data_of_magnetization_abs1)));
    res_.push_back(Variance((valid_data_of_magnetization_abs2)));
    res_.push_back(Variance((valid_data_of_magnetization_abs3)));
  }



  // stiffness
//  res.push_back(Mean(std::vector(stiffness_data.begin() + data_size / 2, stiffness_data.begin() + data_size)));

  //stiffness correction
//  res.push_back(0.0);

  // zig-zag anti-ferromagnetic magnetization
  std::array<std::vector<double>, DimDof> valid_data_of_magnetization_squareQ1;
  std::array<std::vector<double>, DimDof> valid_data_of_magnetization_squareQ2;
  std::array<std::vector<double>, DimDof> valid_data_of_magnetization_squareQ3;
  // Q1
  double M2Q1 = 0.0;
  for (size_t dof = 0; dof < DimDof; dof++) {
    valid_data_of_magnetization_squareQ1[dof] = std::vector<double>(sweeps / 2, 0.0);
    for (size_t i = 0; i < valid_data_of_magnetization_squareQ1[dof].size(); i++) {
      valid_data_of_magnetization_squareQ1[dof][i] =
          zigzag_af_magnetization_data1[dof][i + sweeps / 2] * zigzag_af_magnetization_data1[dof][i + sweeps / 2];
    }
    M2Q1 += Mean(valid_data_of_magnetization_squareQ1[dof]);
    res_.push_back(Mean(valid_data_of_magnetization_squareQ1[dof]) / N / N);
  }
  // Q2
  double M2Q2 = 0.0;
  for (size_t dof = 0; dof < DimDof; dof++) {
    valid_data_of_magnetization_squareQ2[dof] = std::vector<double>(sweeps / 2, 0.0);
    for (size_t i = 0; i < valid_data_of_magnetization_squareQ2[dof].size(); i++) {
      valid_data_of_magnetization_squareQ2[dof][i] =
          zigzag_af_magnetization_data2[dof][i + sweeps / 2] * zigzag_af_magnetization_data2[dof][i + sweeps / 2];
    }
    M2Q2 += Mean(valid_data_of_magnetization_squareQ2[dof]);
    res_.push_back(Mean(valid_data_of_magnetization_squareQ2[dof]) / N / N);
  }
  // Q3
  double M2Q3 = 0.0;
  for (size_t dof = 0; dof < DimDof; dof++) {
    valid_data_of_magnetization_squareQ3[dof] = std::vector<double>(sweeps / 2, 0.0);
    for (size_t i = 0; i < valid_data_of_magnetization_squareQ3[dof].size(); i++) {
      valid_data_of_magnetization_squareQ3[dof][i] =
          zigzag_af_magnetization_data3[dof][i + sweeps / 2] * zigzag_af_magnetization_data3[dof][i + sweeps / 2];
    }
    M2Q3 += Mean(valid_data_of_magnetization_squareQ3[dof]);
    res_.push_back(Mean(valid_data_of_magnetization_squareQ3[dof]) / N / N);
  }


  //Binder ratio
  std::vector<double> valid_data_of_M2_Q1(sweeps / 2, 0.0), valid_data_of_M2_Q2(sweeps / 2, 0.0),
      valid_data_of_M2_Q3(sweeps / 2, 0.0);
  std::vector<double> valid_data_of_magnetization_quarticQ1(sweeps / 2, 0.0),
      valid_data_of_magnetization_quarticQ2(sweeps / 2, 0.0),
      valid_data_of_magnetization_quarticQ3(sweeps / 2, 0.0);

  double binder_ratio1, binder_ratio2, binder_ratio3;
  //Q1
  for (size_t i = 0; i < valid_data_of_magnetization_quarticQ1.size(); i++) {
    double sum(0.0);
    for (size_t dof = 0; dof < DimDof; dof++) {
      sum += valid_data_of_magnetization_squareQ1[dof][i];
    }
    valid_data_of_M2_Q1[i] = sum;
    valid_data_of_magnetization_quarticQ1[i] = sum * sum;
  }
  double M4Q1 = Mean(valid_data_of_magnetization_quarticQ1);
  binder_ratio1 = M4Q1 / M2Q1 / M2Q1;
  //Q2
  for (size_t i = 0; i < valid_data_of_magnetization_quarticQ2.size(); i++) {
    double sum(0.0);
    for (size_t dof = 0; dof < DimDof; dof++) {
      sum += valid_data_of_magnetization_squareQ2[dof][i];
    }
    valid_data_of_M2_Q2[i] = sum;
    valid_data_of_magnetization_quarticQ2[i] = sum * sum;
  }
  double M4Q2 = Mean(valid_data_of_magnetization_quarticQ2);
  binder_ratio2 = M4Q2 / M2Q2 / M2Q2;
  //Q3
  for (size_t i = 0; i < valid_data_of_magnetization_quarticQ3.size(); i++) {
    double sum(0.0);
    for (size_t dof = 0; dof < DimDof; dof++) {
      sum += valid_data_of_magnetization_squareQ3[dof][i];
    }
    valid_data_of_M2_Q3[i] = sum;
    valid_data_of_magnetization_quarticQ3[i] = sum * sum;
  }
  double M4Q3 = Mean(valid_data_of_magnetization_quarticQ3);
  binder_ratio3 = M4Q3 / M2Q3 / M2Q3;

  std::vector<double> valid_data_of_M2(sweeps / 2, 0.0); //m^2 = m_1^2 + m_2^2 + m_3^2
  std::vector<double> valid_data_of_magnetization_quartic(sweeps / 2, 0.0); // m^4
  for (size_t i = 0; i < valid_data_of_M2.size(); i++) {
    valid_data_of_M2[i] = valid_data_of_M2_Q1[i] + valid_data_of_M2_Q2[i] + valid_data_of_M2_Q3[i];
    valid_data_of_magnetization_quartic[i] = valid_data_of_M2[i] * valid_data_of_M2[i];
  }
  double M4 = Mean(valid_data_of_magnetization_quartic);
  double M2 = Mean(valid_data_of_M2);

  res_.push_back(binder_ratio1);
  res_.push_back(binder_ratio2);
  res_.push_back(binder_ratio3);
  res_.push_back(M4 / M2 / M2);


  //C3 symmetry breaking
  std::vector<double> data_of_c3_order_parameter(sweeps / 2, 0.0);
  std::vector<double> data_of_c3_order_parameter_abs(sweeps / 2, 0.0);
  for (size_t i = 0; i < data_of_c3_order_parameter.size(); i++) {
    data_of_c3_order_parameter[i] = valid_data_of_M2_Q1[i] + valid_data_of_M2_Q2[i] - 2 * valid_data_of_M2_Q3[i];
    data_of_c3_order_parameter_abs[i] = abs(data_of_c3_order_parameter[i]);
  }
  res_.push_back(Variance(data_of_c3_order_parameter)); // <b^2>, assume b = 0.0;
  res_.push_back(Variance(data_of_c3_order_parameter_abs)); // susceptibility of C3 order parameter b


  double msquare, mquartic;
  bool has_gs_data = true;
  if (has_gs_data) {
    std::vector<double> m2_data(sweeps / 2);
    std::vector<double> m4_data(sweeps / 2);
    for (size_t i = 0; i < sweeps / 2; i++) {
      auto &real = complex_xy_order_parameter_real_[i + sweeps / 2];
      auto &imag = complex_xy_order_parameter_imag_[i + sweeps / 2];
      m2_data[i] = real * real + imag * imag;
      m4_data[i] = m2_data[i] * m2_data[i];
    }
    msquare = Mean(m2_data);
    mquartic = Mean(m4_data);
    res_.push_back(msquare);

    double binder_ratio_m = mquartic / msquare / msquare;
    res_.push_back(binder_ratio_m);
  }

  for (size_t i = 0; i < res.size(); i++) {
    std::cout << "res[" << i << "] = " << res[i] << std::endl;
  }
  DumpData("summary" + mc_params.filename_postfix, res);
  statistic_timer.PrintElapsed();

  endTime = clock();
  std::cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  MPI_Finalize();
  return 0;
}

