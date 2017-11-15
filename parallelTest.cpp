
#include <iostream>
#include "Lattice.h"
#include "TRandom3.h"
#include "TStopwatch.h"

using namespace std;

void parallelTest(){

  gRandom -> SetSeed(4357);

  int nn = 0;
  int tt = 4;
  cout << "Lattice N: ";
  cin >> nn;
  cout << "Threads: ";
  cin >> tt;
  cout << endl;
  cout << endl;

  TStopwatch timer;
  Lattice a(nn,1,2);
  Lattice b(nn,2,4);
  Lattice c(nn,3,6);
  
  //Lattice aP(nn,1,2);
  //Lattice bP(nn,2,4);
  //Lattice cP(nn,3,6);

  cout << "1-D, " << nn << " entries, seq.\n";
  timer.Start();
  cout << a.energy() << endl;
  timer.Stop();
  timer.Print();
  cout << endl;
  
  cout << "2-D, " << nn << " entries, seq.\n";
  timer.Start();
  cout << b.energy() << endl;
  timer.Stop();
  timer.Print();
  cout << endl;

  cout << "3-D, " << nn << " entries, seq.\n";
  timer.Start();
  cout << c.energy() << endl;
  timer.Stop();
  timer.Print();
  cout << endl;

  cout << "1-D, " << nn << " entries, "<< tt << " threads\n";
  timer.Start();
  cout << a.energyParallel(tt)<< endl;
  timer.Stop();
  timer.Print();
  cout << endl;

  cout << "2-D, " << nn << " entries, " << tt << " threads\n";
  timer.Start();
  cout << b.energyParallel(tt) << endl;
  timer.Stop();
  timer.Print();
  cout << endl;

  cout << "3-D, " << nn << " entries, "<< tt << " threads\n";
  timer.Start();
  cout << c.energyParallel(tt) << endl;
  timer.Stop(); 
  timer.Print();
  cout << endl;

  cout << "seed" << gRandom->GetSeed() << endl;
}
