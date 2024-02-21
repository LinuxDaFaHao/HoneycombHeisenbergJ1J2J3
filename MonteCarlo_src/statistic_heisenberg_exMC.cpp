/*
* Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2022/12/02.
*
* Description: Honeycomb lattice, U1 symmetric, Heisenberg model, exchange monte carlo
*
*/


#include <assert.h>
#include "ExchangeMCExecutor.h"
#include <iostream>
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
  clock_t startTime, endTime;
  startTime = clock();

  CaseParams params(argv[1]);
  size_t sweeps = params.Sweeps;
  std::string geometry = params.Geometry;
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;
  size_t N = 2 * Lx * Ly;
  if (argc < 6) {
    std::cout << "lack arguments" << std::endl;
    exit(1);
  }

  double beta1 = std::atof(argv[2]);
  double beta2 = std::atof(argv[3]);
  size_t thread = std::atoi(argv[4]);
  std::string path_to_old_res = std::string(argv[5]);

  double beta_max = std::max(beta1, beta2);
  double beta_min = std::min(beta1, beta2);

  std::cout << "Read argument beta_max = " << beta_max << std::endl;
  std::cout << "              beta_min = " << beta_min << std::endl;
  std::cout << "              thread_num = " << thread << std::endl;

  double J1Matrix[3][3], J2Matrix[3][3], J3Matrix[3][3], DMatrix[3][3];
  Vec2Mat(params.J1, J1Matrix);
  auto J11 = CouplingStructure<DimDof>(J1Matrix);

  Vec2Mat(params.J2, J2Matrix);
  auto J21 = CouplingStructure<DimDof>(J2Matrix);

  Vec2Mat(params.J3, J3Matrix);
  auto J31 = CouplingStructure<DimDof>(J3Matrix);

  Vec2Mat(params.D, DMatrix);
  bool Interaction_SU2 = true;
  Interaction_SU2 = Interaction_SU2 && J11.IsIsometry();
  Interaction_SU2 = Interaction_SU2 && J21.IsIsometry();
  Interaction_SU2 = Interaction_SU2 && J31.IsIsometry();

  if (!Interaction_SU2) {
    std::cout << "coupling has anistropy." << std::endl;
    exit(1);
  }

  if (geometry != "Honeycomb") {
    std::cout << "beyond assumption geometry." << std::endl;
    exit(1);
  }

  ExchangeMCParams mc_params;
  mc_params.thread_num = thread;
  mc_params.adjust_temperature_times = 5;
  mc_params.adjust_temperature_samples = sweeps / 100;
  mc_params.sweeps = sweeps;
  mc_params.sample_interval = params.SampleInterval;
  mc_params.warmup_sample_num = params.Warmup_Samples;
  mc_params.exchange_interval = params.ExchangeInterval;
  mc_params.print_interval = 100;
  mc_params.filename_postfix = "hei" + params.Geometry
      + "J1zz" + std::to_string(J1Matrix[2][2])
      + "J2zz" + std::to_string(J2Matrix[2][2])
      + "J3zz" + std::to_string(J3Matrix[2][2])
      + "Dzz" + std::to_string(DMatrix[2][2])
      + "beta_max" + std::to_string(beta_max)
      + "beta_min" + std::to_string(beta_min)
      + "L" + std::to_string(Lx);

  qlten::Timer load_data_timer("load_data");
  StatisticResults<DimDof> results_(thread);
  results_.LoadData(path_to_old_res + "results" + mc_params.filename_postfix);
  size_t data_size = mc_params.sweeps;

  std::vector<double> &beta_set_ = results_.beta;
  std::vector<std::array<std::vector<double>, DimDof>> zigzag_af_magnetization1_set_(thread);
  std::vector<std::array<std::vector<double>, DimDof>> zigzag_af_magnetization2_set_(thread);
  std::vector<std::array<std::vector<double>, DimDof>> zigzag_af_magnetization3_set_(thread);

#pragma omp parallel for default(none) \
                shared(zigzag_af_magnetization1_set_, zigzag_af_magnetization2_set_, \
                        zigzag_af_magnetization3_set_, mc_params, beta_set_, N,  data_size, results_)\
                num_threads(mc_params.thread_num)\
                schedule(static)
  for (size_t replica_id = 0; replica_id < mc_params.thread_num; replica_id++) {
    for (size_t i = 0; i < DimDof; i++) {
      zigzag_af_magnetization1_set_[replica_id][i] =
          LoadData<double>("magnetization1" + std::to_string(beta_set_[replica_id]) + std::to_string(i)
                               + mc_params.filename_postfix, data_size);
      zigzag_af_magnetization2_set_[replica_id][i] =
          LoadData<double>("magnetization2" + std::to_string(beta_set_[replica_id]) + std::to_string(i)
                               + mc_params.filename_postfix, data_size);
      zigzag_af_magnetization3_set_[replica_id][i] =
          LoadData<double>("magnetization3" + std::to_string(beta_set_[replica_id]) + std::to_string(i)
                               + mc_params.filename_postfix, data_size);

    }
  }

  load_data_timer.PrintElapsed();

  qlten::Timer statistic_timer("statistic");

  size_t warmup_sample_num = mc_params.warmup_sample_num;
  size_t N_square = N * N;
  size_t N_quartic = N_square * N_square;
#pragma omp parallel for default(none) \
                shared(zigzag_af_magnetization1_set_, zigzag_af_magnetization2_set_, \
                        zigzag_af_magnetization3_set_, mc_params, beta_set_, N, N_square, N_quartic, sweeps, warmup_sample_num, results_)\
                num_threads(mc_params.thread_num)\
                schedule(static)
  for (size_t replica_id = 0; replica_id < mc_params.thread_num; replica_id++) {
    std::array<std::vector<double>, DimDof> thermalized_magnetization_component_square_Q1_data,
        thermalized_magnetization_component_square_Q2_data, thermalized_magnetization_component_square_Q3_data;

    for (size_t dof = 0; dof < DimDof; dof++) {
      thermalized_magnetization_component_square_Q1_data[dof] =
          std::vector<double>(sweeps - warmup_sample_num, 0.0);
      for (size_t i = 0; i < thermalized_magnetization_component_square_Q1_data[dof].size(); i++) {
        thermalized_magnetization_component_square_Q1_data[dof][i] =
            zigzag_af_magnetization1_set_[replica_id][dof][i + warmup_sample_num]
                * zigzag_af_magnetization1_set_[replica_id][dof][i + warmup_sample_num];
      }
    }
    for (size_t dof = 0; dof < DimDof; dof++) {
      thermalized_magnetization_component_square_Q2_data[dof] =
          std::vector<double>(sweeps - warmup_sample_num, 0.0);
      for (size_t i = 0; i < thermalized_magnetization_component_square_Q2_data[dof].size(); i++) {
        thermalized_magnetization_component_square_Q2_data[dof][i] =
            zigzag_af_magnetization2_set_[replica_id][dof][i + warmup_sample_num]
                * zigzag_af_magnetization2_set_[replica_id][dof][i + warmup_sample_num];
      }
    }

    for (size_t dof = 0; dof < DimDof; dof++) {
      thermalized_magnetization_component_square_Q3_data[dof] =
          std::vector<double>(sweeps - warmup_sample_num, 0.0);
      for (size_t i = 0; i < thermalized_magnetization_component_square_Q3_data[dof].size(); i++) {
        thermalized_magnetization_component_square_Q3_data[dof][i] =
            zigzag_af_magnetization3_set_[replica_id][dof][i + warmup_sample_num]
                * zigzag_af_magnetization3_set_[replica_id][dof][i + warmup_sample_num];
      }
    }

    std::vector<double> thermalized_magnetization_square_q1_data(sweeps - warmup_sample_num, 0.0),
        thermalized_magnetization_square_q2_data(sweeps - warmup_sample_num, 0.0),
        thermalized_magnetization_square_q3_data(sweeps - warmup_sample_num, 0.0);

    for (size_t i = 0; i < sweeps - warmup_sample_num; i++) {
      double sum(0.0);
      for (size_t dof = 0; dof < DimDof; dof++) {
        sum += thermalized_magnetization_component_square_Q1_data[dof][i];
      }
      thermalized_magnetization_square_q1_data[i] = sum;
    }

    for (size_t i = 0; i < sweeps - warmup_sample_num; i++) {
      double sum(0.0);
      for (size_t dof = 0; dof < DimDof; dof++) {
        sum += thermalized_magnetization_component_square_Q2_data[dof][i];
      }
      thermalized_magnetization_square_q2_data[i] = sum;
    }

    for (size_t i = 0; i < sweeps - warmup_sample_num; i++) {
      double sum(0.0);
      for (size_t dof = 0; dof < DimDof; dof++) {
        sum += thermalized_magnetization_component_square_Q3_data[dof][i];
      }
      thermalized_magnetization_square_q3_data[i] = sum;
    }
    //C3 symmetry order parameter and its susceptibility
    using complexT = std::complex<double>;
    std::vector<complexT> thermalized_data_of_c3_order_parameter(sweeps - warmup_sample_num, 0.0);
    std::vector<double> thermalized_data_of_c3_order_parameter_abs(sweeps - warmup_sample_num, 0.0);
    std::vector<double> thermalized_data_of_c3_order_parameter_square(sweeps - warmup_sample_num, 0.0);
    for (size_t i = 0; i < sweeps - warmup_sample_num; i++) {
      thermalized_data_of_c3_order_parameter[i] = thermalized_magnetization_square_q1_data[i]
          + std::exp(complexT(0, 2 * M_PI / 3)) * thermalized_magnetization_square_q2_data[i]
          + std::exp(complexT(0, 4 * M_PI / 3)) * thermalized_magnetization_square_q3_data[i];
      thermalized_data_of_c3_order_parameter_abs[i] = abs(thermalized_data_of_c3_order_parameter[i]);
      thermalized_data_of_c3_order_parameter_square[i] =
          thermalized_data_of_c3_order_parameter_abs[i] * thermalized_data_of_c3_order_parameter_abs[i];
    }
    double B_c3_square = Mean(thermalized_data_of_c3_order_parameter_square);
    double B_c3_abs = Mean(thermalized_data_of_c3_order_parameter_abs);
    results_.c3_order_parameter_square[replica_id] = B_c3_square / N_quartic;
    results_.c3_order_parameter_susceptibility[replica_id] =
        (B_c3_square - B_c3_abs * B_c3_abs) / N_square * beta_set_[replica_id];
    results_.binder_ratio_c3[replica_id] = B_c3_square / B_c3_abs / B_c3_abs;
  }
  statistic_timer.PrintElapsed();

  results_.DumpData("results" + mc_params.filename_postfix);
  endTime = clock();
  std::cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  return 0;
}

