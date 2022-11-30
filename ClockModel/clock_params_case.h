//
// Created by Hao-Xin on 2022/5/28.
//

#ifndef HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCK_PARAMS_CASE_H
#define HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCK_PARAMS_CASE_H


#include "gqmps2/case_params_parser.h"
using gqmps2::CaseParamsParserBasic;

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    Geometry = ParseStr("Geometry");
    Ly = ParseInt("L");
    Lx = ParseInt("L");
    J = ParseDouble("J");
    hp = ParseDouble("hp");
    Sweeps = ParseInt("Sweeps");
    ClusterRadius = ParseInt("ClusterRadius");
    SampleInterval = ParseInt("SampleInterval");
    beta = ParseDouble("beta");
  }

  std::string Geometry;
  size_t Ly;
  size_t Lx;
  double J;
  double hp;
  size_t Sweeps;
  size_t ClusterRadius;
  size_t SampleInterval;
  double beta;
};


#endif //HONEYCOMBHEISENBERG_CLOCKMODEL_CLOCK_PARAMS_CASE_H
