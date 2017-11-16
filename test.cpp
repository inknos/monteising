//
//  test.cpp
//
//
//  Created by Federico on 12/11/17.
//

#include <iostream>
#include <fstream>

#include "Lattice.h"
#include "TRandom3.h"
#include "TStopwatch.h"

using namespace std;

void test(){
  gRandom -> SetSeed(357);
  
  TStopwatch timer;
  Lattice uni(10,1,2);
  Lattice bi(10,2,4);
  Lattice tri(10,3,6);

  ofstream file;

  uni.printLattice("uno.csv");
  bi.printLattice("due.csv");
  tri.printLattice("tre.csv");
  
  cout << uni << endl;
  cout << bi << endl;
  cout << tri << endl;

  
  timer.Start();
  cout << uni.energy() << endl;
  cout << bi.energy() << endl;
  cout << tri.energy() << endl;
  timer.Stop();
  
  cout << "seed" << gRandom->GetSeed() << endl;
  timer.Print();
}











