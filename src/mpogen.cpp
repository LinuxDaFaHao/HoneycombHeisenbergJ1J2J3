/*
 * File Name: mpogen.cpp
 * Description: Generate and dump mpo for the spin model
 * Created by Hao-Xin on 2022/04/09.
 *
 */

#include "DefSpinOne.h"
#include "operators.h"
#include "params_case.h"
#include "gqmps2/gqmps2.h"
#include "honeycomb_latiice.h"

#include <ctime>

using namespace std;
using namespace gqten;
using namespace gqmps2;


int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);
  const size_t Lx = params.Lx;
  const size_t Ly = params.Ly;
  const size_t N = 2 * Lx * params.Ly;
  cout << "The total number of sites: " << N << endl;
  double J1 = params.J1;
  double J2 = params.J2;
  double J3 = params.J3;
  cout << "Model parameter: J1 = " << J1 << ",\t J2 = " << J2 << ",\t J3 = " << J3 << endl;
  clock_t startTime,endTime;
  startTime = clock();


  if (!IsPathExist(kMpoPath)) {
    CreatPath(kMpoPath);
  }

  using namespace spin_one_model;
  OperatorInitial();
  const SiteVec<TenElemT, U1QN> sites=SiteVec<TenElemT, U1QN>(N, pb_out);
  gqmps2::MPOGenerator<TenElemT, U1QN> mpo_gen(sites, qn0);

  HoneycombLattice* plattice = nullptr;
  if(params.Geometry == "XC"){
    plattice = new HoneycombXCCylinder(Ly, params.Lx);
  } else if(params.Geometry == "YC"){
    plattice = new HoneycombYCCylinder(Ly, params.Lx);
  } else {
    std::cout << "Do not support geometry now. Exit program." << std::endl;
    exit(1);
  }

  for(auto& link : plattice->nearest_neighbor_links) {
    size_t site1 = std::get<0>(link),
        site2 = std::get<1>(link);
    if(site1 > site2) {
      std::swap(site1, site2);
    }
    mpo_gen.AddTerm(J1, sz, site1, sz, site2);
    mpo_gen.AddTerm(J1/2, sp, site1, sm, site2);
    mpo_gen.AddTerm(J1/2, sm, site1, sp, site2);
    std::cout << "add Spin-Spin interaction for bond (" << site1 << ", " << site2 << "); \n";
  }
  std::cout << "Number of NN-bond = " << plattice->nearest_neighbor_links.size() << std::endl;

  for(auto& link : plattice->next_nearest_neighbor_links) {
    size_t site1 = std::get<0>(link),
        site2 = std::get<1>(link);
    if(site1 > site2) {
      std::swap(site1, site2);
    }
    mpo_gen.AddTerm(J2, sz, site1, sz, site2);
    mpo_gen.AddTerm(J2/2, sp, site1, sm, site2);
    mpo_gen.AddTerm(J2/2, sm, site1, sp, site2);
    std::cout << "add Spin-Spin interaction for NNN (" << site1 << ", " << site2 << "); \n";
  }

  std::cout << "Number of NNN-pairs = " << plattice->next_nearest_neighbor_links.size() << std::endl;

  for(auto& link : plattice->third_nearest_neighbor_links) {
    size_t site1 = std::get<0>(link),
        site2 = std::get<1>(link);
    if(site1 > site2) {
      std::swap(site1, site2);
    }
    mpo_gen.AddTerm(J3, sz, site1, sz, site2);
    mpo_gen.AddTerm(J3/2, sp, site1, sm, site2);
    mpo_gen.AddTerm(J3/2, sm, site1, sp, site2);
    std::cout << "add Spin-Spin interaction for NNN (" << site1 << ", " << site2 << "); \n";
  }

  std::cout << "Number of TNN-pairs = " << plattice->third_nearest_neighbor_links.size() << std::endl;

  auto mpo = mpo_gen.Gen();
  cout << "FiniteMPO generated." << endl;

  for(size_t i=0; i<mpo.size();i++){
    std::string filename = kMpoPath + "/" +
        kMpoTenBaseName + std::to_string(i)
        + "." + kGQTenFileSuffix;
    mpo.DumpTen(i,filename);
  }

  endTime = clock();
  cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
