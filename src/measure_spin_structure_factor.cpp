#include "gqmps2/gqmps2.h"
#include "gqten/gqten.h"
#include <ctime>
#include "DefSpinOne.h"
#include "operators.h"
#include "params_case.h"
#include "myutil.h"
#include "my_measure.h"       //MeasureTwoSiteOp

using namespace gqmps2;
using namespace gqten;
using namespace std;

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;

  CaseParams params(argv[1]);
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;
  size_t N = 2 * Lx * Ly;

  clock_t startTime, endTime;
  startTime = clock();

  gqten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
  gqten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  using namespace spin_one_model;
  OperatorInitial();
  const SiteVec<TenElemT, U1QN> sites = SiteVec<TenElemT, U1QN>(N, pb_out);

  using FiniteMPST = gqmps2::FiniteMPS<TenElemT, U1QN>;
  FiniteMPST mps(sites);
//  mps.Load();
//  cout << "bond dimension of middle mps = ";
//  cout << mps[N / 2].GetShape()[0] << endl;

  /******** define the measure_tasks ********/
  std::vector<MeasureGroupTask> measure_tasks;
  measure_tasks.reserve(N);
  size_t begin_x = Lx / 2;
  size_t end_x = 3 * Lx / 2 + 1;
  for (size_t i = begin_x * Ly; i < end_x * Ly; i++) {
    const size_t site1 = i;
    std::vector<size_t> site2;
    site2.reserve(N - i);
    for (size_t j = i + 1; j < end_x * Ly; j++) {
      site2.push_back(j);
    }
    measure_tasks.push_back(MeasureGroupTask(site1, site2));
  }

  Timer one_site_timer("measure spin_ structure factors");
  MeasureTwoSiteOp(mps, kMpsPath, sz, sz, measure_tasks, "zzsf", world);
  MeasureTwoSiteOp(mps, kMpsPath, sp, sm, measure_tasks, "pmsf", world);
  MeasureTwoSiteOp(mps, kMpsPath, sm, sp, measure_tasks, "mpsf", world);
  one_site_timer.PrintElapsed();

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}