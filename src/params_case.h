/*
 * File Name: params_case.h
 * Description: Declare CaseParams class used set parameters by users
 * Created by Hao-Xin on 2022/04/09.
 *
 */


#ifndef HONEYCOMBHEISENBERGJ1J2J3_SRC_PARAMS_CASE_H
#define HONEYCOMBHEISENBERGJ1J2J3_SRC_PARAMS_CASE_H

#include "gqmps2/case_params_parser.h"
using gqmps2::CaseParamsParserBasic;

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    Geometry = ParseStr("Geometry");
    Ly = ParseInt("Ly");
    Lx = ParseInt("Lx");
    J1 = ParseDouble("J1");
    J2 = ParseDouble("J2");
    J3 = ParseDouble("J3");
    Sweeps = ParseInt("Sweeps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
    tau = ParseDouble("tau");
    M = ParseInt("M");
    ConvergeTolerance = ParseDouble("ConvergeTolerance");
    TaylorOrder = ParseInt("TaylorOrder");
    TaylorErr = ParseDouble("TaylorErr");
    Threads = ParseInt("Threads");
    Perturbation=ParseDouble("Perturbation");
    wavelength = ParseInt("wavelength");
    noise = ParseDoubleVec("noise");
    SymmetryMode = ParseInt("SymmetryMode");
  }

  std::string Geometry;
  size_t Ly;
  size_t Lx;
  double J1;
  double J2;
  double J3;
  size_t Sweeps;
  size_t Dmin;
  size_t Dmax;
  double CutOff;
  double LanczErr;
  size_t MaxLanczIter;
  double tau;
  size_t M;
  double ConvergeTolerance;
  size_t TaylorOrder;
  double TaylorErr;
  size_t Threads;
  double Perturbation;
  size_t wavelength;
  std::vector<double> noise;
  size_t SymmetryMode;//useless upto now
};

#endif //HONEYCOMBHEISENBERGJ1J2J3_SRC_PARAMS_CASE_H
