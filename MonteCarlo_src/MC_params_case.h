//
// Created by Hao-Xin on 2022/5/28.
//

#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_MC_PARAMS_CASE_H_
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_MC_PARAMS_CASE_H_


#include "gqmps2/case_params_parser.h"
using gqmps2::CaseParamsParserBasic;

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    Geometry = ParseStr("Geometry");
    Ly = ParseInt("L");
    Lx = ParseInt("L");
    J1 = ParseDoubleVec("J1");
    J2 = ParseDoubleVec("J2");
    J3 = ParseDoubleVec("J3");
    D = ParseDoubleVec("D");
    Sweeps = ParseInt("Sweeps");
    SampleInterval = ParseInt("SampleInterval");
    beta = ParseDouble("beta");
  }

  std::string Geometry;
  size_t Ly;
  size_t Lx;
  std::vector<double> J1;
  std::vector<double> J2;
  std::vector<double> J3;
  std::vector<double> D;
  size_t Sweeps;
  size_t SampleInterval;
  double beta;
};


#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_MC_PARAMS_CASE_H_
