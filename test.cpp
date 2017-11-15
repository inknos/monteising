//
//  test.cpp
//
//
//  Created by Federico on 12/11/17.
//

#include <iostream>
#include "Lattice.h"
#include "TRandom3.h"
#include "TStopwatch.h"

using namespace std;

void test(){
  gRandom -> SetSeed(4357);
  
  TStopwatch timer;
  Lattice uni(5,1,2);
  Lattice bi(3,2,4);
  Lattice tri(3,3,6);
  Lattice uui(5,1,2);
  Lattice bbi(3,2,4);
  Lattice ttri(100,3,6);

  cout << uni << endl;
  cout << uui << endl;
  cout << bi << endl;
  cout << bbi  << endl;
  cout << "tre\n";
  cout << tri << endl;
  cout << "tre bis\n";
  //cout << ttri<< endl;
  
  timer.Start();
  cout << uni.energy() << endl;
  cout << uui.energy() << endl;
  cout << bi.energy() << endl;
  cout << bbi.energy()<< endl;
  cout << tri.energy() << endl;
  cout <<ttri.energy() << endl;
  timer.Stop();
  
  cout << "seed" << gRandom->GetSeed() << endl;
  timer.Print();
}











