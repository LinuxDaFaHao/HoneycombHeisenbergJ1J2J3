//
// Created by Hao-Xin on 2022/4/9.
//

#ifndef HONEYCOMBHEISENBERGJ1J2J3_SRC_HONEYCOMB_LATIICE_H
#define HONEYCOMBHEISENBERGJ1J2J3_SRC_HONEYCOMB_LATIICE_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <assert.h>

using std::vector;
using std::tuple;
using std::min;
using std::max;

using Link = tuple<size_t, size_t>;
//using OritationedLink = tuple<size_t, size_t>;

class HoneycombLattice {
 public:
  HoneycombLattice(void) = default;
  HoneycombLattice(const size_t Ly, const size_t Lx) : Ly(Ly), Lx(Lx), N(2*Lx*Ly) {}

  const HoneycombLattice& operator=(const HoneycombLattice& lhs) {
    Ly = lhs.Ly;
    Lx = lhs.Lx;
    N = lhs.N;
    nearest_neighbor_links = lhs.nearest_neighbor_links;
    next_nearest_neighbor_links = lhs.next_nearest_neighbor_links;
    third_nearest_neighbor_links = lhs.third_nearest_neighbor_links;
    return *this;
  }

  inline void Print(void) {
    std::cout << "System size (Lx, Ly) = ( "
              << Lx << ", "
              << Ly << ")" << "\n";
  }


  size_t Ly;
  size_t Lx;
  size_t N;
  vector<Link> nearest_neighbor_links;
  vector<Link> next_nearest_neighbor_links;
  vector<Link> third_nearest_neighbor_links;
};

// For the definition of geometry please find in https://online.kitp.ucsb.edu/online/fragnets12/honeydmrg/pdf1/Zhu_Honeycomb_Fragnets12_KITP.pdf
class HoneycombYCCylinder : public HoneycombLattice {
 public:
  HoneycombYCCylinder() = default;
  HoneycombYCCylinder(const size_t Ly, const size_t Lx);
};

HoneycombYCCylinder::HoneycombYCCylinder(const size_t Ly, const size_t Lx) :
HoneycombLattice(Ly, Lx) {
  if(Ly % 2 != 0) {
    std::cout << "Ly = " << Ly << std::endl;
    std::cout << "Unexpected odd Ly!" << std::endl;
    exit(1);
  }
  nearest_neighbor_links.reserve( 3 * N / 2 );
  next_nearest_neighbor_links.reserve( 3 * N);
  third_nearest_neighbor_links.reserve( 3 * N / 2);
  for (size_t i = 0; i < N; ++i) {
    const size_t y = i % Ly; //y coordinate of site i
    const size_t x = i / Ly; //x coordinate of site i
    const size_t Tx = y + Ly * ((x + 1) % (2 * Lx)); //x-directional translation site of site i
    const size_t Txx = y + Ly * ((x + 2) % (2 * Lx));
    const size_t Ty = (y + 1) % Ly + Ly * x;   //y-directional translation site of site i
    const size_t Txy = (y + 1) % Ly + Ly * ((x + 1) % (2 * Lx)); //right-up-directional translation site of site i
    const size_t Txmy = (y - 1 + Ly) % Ly + Ly * ((x + 1) % (2 * Lx));
    const size_t Txxy = (y + 1) % Ly  + Ly * ((x + 2) % (2 * Lx));
    const size_t Txxmy = (y - 1 + Ly) % Ly  + Ly * ((x + 2) % (2 * Lx));


    //J1
    if( x < 2 * Lx - 1) {
      if( (x + y) % 2 == 0  ) {
        nearest_neighbor_links.push_back( Link{i, Ty});
        nearest_neighbor_links.push_back( Link{i, Tx});
      } else {
        nearest_neighbor_links.push_back( Link{i, Tx});
      }
    } else {
      if ( y % 2 == 1) {
        nearest_neighbor_links.push_back( Link{i, Ty});
      }
    }


    //J2
    if( x < 2 * Lx - 2) {
      next_nearest_neighbor_links.push_back( Link{i, Txx});
      next_nearest_neighbor_links.push_back( Link{i, Txy});
      next_nearest_neighbor_links.push_back( Link{i, Txmy});
    } else if ( x == 2 * Lx - 2) {
      next_nearest_neighbor_links.push_back( Link{i, Txy});
      next_nearest_neighbor_links.push_back( Link{i, Txmy});
    } else {
      //do nothing
    }

    //J3
    if( x < 2 * Lx - 2) {
      if( (x + y) % 2 == 1  ) {
        third_nearest_neighbor_links.push_back( Link{i, Ty});
        third_nearest_neighbor_links.push_back( Link{i, Txxmy});
      } else {
        third_nearest_neighbor_links.push_back( Link{i, Txxy});
      }
    } else if ( (x + y) % 2 == 1) {
        third_nearest_neighbor_links.push_back( Link{i, Ty});
    } else {
      //do nothing
    }

  }
}

// For the definition of geometry please find in https://online.kitp.ucsb.edu/online/fragnets12/honeydmrg/pdf1/Zhu_Honeycomb_Fragnets12_KITP.pdf
// XC-Ly Cylinder
class HoneycombXCCylinder : public HoneycombLattice {
 public:
  HoneycombXCCylinder() = default;
  HoneycombXCCylinder(const size_t Ly, const size_t Lx);
};

HoneycombXCCylinder::HoneycombXCCylinder(const size_t Ly, const size_t Lx) :
    HoneycombLattice(Ly, Lx) {
  if(Ly % 2 != 0) {
    std::cout << "Ly = " << Ly << std::endl;
    std::cout << "Unexpected odd Ly!" << std::endl;
    exit(1);
  }
  nearest_neighbor_links.reserve( 3 * N / 2 );
  next_nearest_neighbor_links.reserve( 3 * N);
  third_nearest_neighbor_links.reserve( 3 * N / 2);
  for (size_t i = 0; i < N; ++i) {
    const size_t y = i % Ly; //y coordinate of site i
    const size_t x = i / Ly; //x coordinate of site i
    const size_t Tx = y + Ly * ((x + 1) % (2*Lx) ); //x-directional translation site of site i
    const size_t Ty = (y + 1) % Ly + Ly * x;   //y-directional translation site of site i
    const size_t Tyy = (y + 2) % Ly + Ly * x;
    const size_t Txy = (y + 1) % Ly + Ly * ((x + 1) % (2*Lx) ); //right-up-directional translation site of site i
    const size_t Txyy = (y + 2) % Ly + Ly * ((x + 1) % (2*Lx) );
    const size_t Txmyy = (y - 2 + Ly) % Ly + Ly * ((x + 1) % (2*Lx) );
    const size_t Txmy = (y - 1 + Ly) % Ly + Ly * ((x + 1) % (2*Lx) );

    // J1
    if( x == 0) {
      nearest_neighbor_links.push_back( Link{i, Ty});
      if( y % 2 == 0) {
        nearest_neighbor_links.push_back( Link{i, Tx});
      }
    } else if( x == 2 * Lx - 1) {
      nearest_neighbor_links.push_back( Link{i, Ty});
    } else if( (x + y) % 2 == 0  ) {
      nearest_neighbor_links.push_back( Link{i, Ty});
      nearest_neighbor_links.push_back( Link{i, Tx});
    } else if( (x + y) % 2 == 1  ) { //else all case
      nearest_neighbor_links.push_back( Link{i, Ty});
    }

    //J2
    if( x < 2 * Lx - 1) {
      next_nearest_neighbor_links.push_back( Link{i, Tyy});
      next_nearest_neighbor_links.push_back( Link{i, Txy});
      next_nearest_neighbor_links.push_back( Link{i, Txmy});
    } else {
      next_nearest_neighbor_links.push_back( Link{i, Tyy});
    }


    //J3
    if( x < 2 * Lx - 1) {
      if( (x + y) % 2 == 0 ) {
        third_nearest_neighbor_links.push_back( Link{i, Txyy});
        third_nearest_neighbor_links.push_back( Link{i, Txmyy});
      } else {
        third_nearest_neighbor_links.push_back( Link{i, Tx});
      }
    }
  }
}

#endif //HONEYCOMBHEISENBERGJ1J2J3_SRC_HONEYCOMB_LATIICE_H
