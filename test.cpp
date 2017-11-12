//
//  test.cpp
//
//
//  Created by Federico on 12/11/17.
//

#include <iostream>
#include "Lattice.h"

using namespace std;

void test(){
  Lattice uni(5,1,2);
  Lattice bi(3,2,4);
  Lattice tri(3,3,6);



  uni.printLattice();
  cout << endl;
  bi.printLattice();
  cout << endl;
  tri.printLattice();
  cout << endl;

  cout << uni.energy() << endl;
  cout << bi.energy() << endl;
  cout << tri.energy() << endl;
}
