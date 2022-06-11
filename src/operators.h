/*
 * File Name: operators.h
 * Description: Declare the operators in spin_-1/2 model
 * Created by Hao-Xin on 2022/04/09.
 *
 */

#ifndef HONEYCOMBHEISENBERGJ1J2J3_SRC_SPIN_OPERATORS_H
#define HONEYCOMBHEISENBERGJ1J2J3_SRC_SPIN_OPERATORS_H

#include "DefSpinOne.h"

namespace spin_one_model {
//Spin-1 operators
extern Tensor sz, sp, sm, id, sz_square;
void OperatorInitial();
}

#endif //HONEYCOMBHEISENBERGJ1J2J3_SRC_SPIN_OPERATORS_H
