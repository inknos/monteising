//
//  test.cpp
//
//
//  Created by Federico on 12/11/17.
//

#include <iostream>
#include <fstream>
#include <ctime>


#include "Lattice.h"
#include "TRandom3.h"
#include "TStopwatch.h"

using namespace std;

void test(){
  gRandom -> SetSeed(time(0));
  
  TStopwatch timer;
  Lattice uni(1000,1,2);
  Lattice bi(1000,2,4);
  Lattice tri(100,3,6);

  ofstream file;

  //uni.printLattice("uno.csv");
  //bi.printLattice("due.csv");
  //tri.printLattice("tre.csv");
  
  //cout << uni << endl;
  //cout << bi << endl;
  //cout << tri << endl;

  
  timer.Start();
  cout << uni.energy() << endl;
  cout << bi.energy() << endl;
  cout << tri.energy() << endl;
  timer.Stop();
  timer.Print();

  
  timer.Start();
  cout << uni.energy2() << endl;
  cout << bi.energy2() << endl;
  cout << tri.energy2() << endl;
  timer.Stop();
  timer.Print();
  
  timer.Start();
  cout << uni.energyFede() << endl;
  cout << bi.energyFede() << endl;
  cout << tri.energyFede() << endl;
  timer.Stop();
  timer.Print();
  
  cout << "seed" << gRandom->GetSeed() << endl;
}











