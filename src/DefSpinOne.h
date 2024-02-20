/*
 * File Name: DefSpinOne.h
 * Description: Define the types, constant and spin_-1 Hilbert Space in a spin_-1 model
 * Created by Hao-Xin on 2022/04/09.
 *
 */

#ifndef HONEYCOMBHEISENBERGJ1J2J3_SRC_SPIN_ONE_H
#define HONEYCOMBHEISENBERGJ1J2J3_SRC_SPIN_ONE_H

#include "qlten/qlten.h"


using TenElemT = qlten::QLTEN_Double;

using U1QN = qlten::special_qn::U1QN; // Sz

using qlten::QLTensor;

const std::string kMpoPath = "mpo";
const std::string kMpoTenBaseName = "mpo_ten";

namespace spin_one_model {
using QNSctT = qlten::QNSector<U1QN>;
using IndexT = qlten::Index<U1QN>;
using Tensor = QLTensor<TenElemT, U1QN>;
const U1QN qn0 = U1QN(0);
const IndexT pb_out = IndexT({   //QNSctT( U1QN(Sz * 3), degeneracy )
                                  QNSctT(U1QN(1), 1),
                                  QNSctT(U1QN(0), 1),
                                  QNSctT(U1QN(-1), 1),
                              },
                              qlten::TenIndexDirType::OUT
);
const IndexT pb_in = qlten::InverseIndex(pb_out);

}//spin_one_model


#endif //HONEYCOMBHEISENBERGJ1J2J3_SRC_SPIN_ONE_H
