/*
 * File Name: operators.cpp
 * Description: Define the operators in spin-1 model
 * Created by Hao-Xin on 2022/04/09.
 *
 */

#include "DefSpinOne.h"
#include <cmath>
using std::sqrt;

namespace spin_one_model {
Tensor sz = Tensor({pb_in, pb_out});
Tensor sp = Tensor({pb_in, pb_out});
Tensor sm = Tensor({pb_in, pb_out});
Tensor id = Tensor({pb_in, pb_out});


void OperatorInitial(){
  static bool initialized = false;
  if(!initialized){
    //suppose operator multiply on the mps, first index is the column index.
    sz({0, 0}) = 1.0;
    sz({1, 1}) = 0.0;
    sz({2, 2}) = -1.0;

    sp({1, 0}) = sqrt(2.0);
    sp({2, 1}) = sqrt(2.0);

    sm({0, 1}) = sqrt(2.0);
    sm({1, 2}) = sqrt(2.0);

    id({0, 0}) = 1;
    id({1, 1}) = 1;
    id({2, 2}) = 1;

    initialized=true;
  }
}

}//namespace spin_model







