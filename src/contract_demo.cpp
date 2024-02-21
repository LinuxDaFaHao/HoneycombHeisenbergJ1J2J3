//
// Created by Hao-Xin on 2023/2/17.
//


#include "qlmps/qlmps.h"
#include <iostream>
#include <vector>
#include "myutil.h"

using std::string;
using std::vector;
using namespace qlmps;
using namespace qlten;
using namespace std;


using TenElemT = qlten::qlten_Complex;

using U1QN = qlten::special_qn::U1QN; // Sz
using QNSctT = qlten::QNSector<U1QN>;
using IndexT = qlten::Index<U1QN>;
using Tensor = QLTensor<TenElemT, U1QN>;
const U1QN qn0 = U1QN(0);
int main(int argc, char *argv[]){
  if(argc != 7) {
    std::cout << "argument number is wrong!" << std::endl;
    exit(-1);
  }
  size_t dima_0 = std::atoi(argv[1]);
  size_t dima_1 = std::atoi(argv[2]);
  size_t dima_2 = std::atoi(argv[3]);

  size_t dimb_0 = std::atoi(argv[4]);
  size_t dimb_1 = std::atoi(argv[5]);

  size_t thread = std::atoi(argv[6]);
  std::cout << "tensor A size = ("  << dima_0
                                    <<", " << dima_1
                                    <<", " << dima_2
                                    <<")" << std::endl;
  std::cout << "tensor B size = ("  << dimb_0
            <<", " << dimb_1
            <<")" << std::endl;
  std::cout << "thread number = " << thread << std::endl;


  qlten::hp_numeric::SetTensorManipulationThreads(thread);

  qlten::Timer generator_tensor_timer("random generate tensors");
  const IndexT pba0 = IndexT({QNSctT(U1QN(0), dima_0),},
                               qlten::TenIndexDirType::IN
  );
  const IndexT pba1 = IndexT({QNSctT(U1QN(0), dima_1),},
                            qlten::TenIndexDirType::IN
  );
  const IndexT pba2 = IndexT({QNSctT(U1QN(0), dima_2),},
                            qlten::TenIndexDirType::OUT
  );
  Tensor tensor_a({pba0, pba1,pba2});
  tensor_a.Random(qn0);
  const IndexT pbb0 = IndexT({QNSctT(U1QN(0), dimb_0),},
                             qlten::TenIndexDirType::IN
  );
  const IndexT pbb1 = IndexT({QNSctT(U1QN(0), dimb_1),},
                             qlten::TenIndexDirType::OUT
  );

  Tensor tensor_b({pbb0, pbb1});
  tensor_b.Random(qn0);
  Tensor res;
  generator_tensor_timer.PrintElapsed();
  qlten::Timer contract_timer("contract");
  Contract(&tensor_a, &tensor_b, {{0},{0}},&res);
  contract_timer.PrintElapsed();
  return 0;
}

