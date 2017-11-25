//
//  test_cooling.cpp
//
//
//  Created by Federico on 25/11/17.
//

#include <stdio.h>
#include <string>
#include "SimulationLattice.h"
#include "Block.h"
#include "TGraph.h"
#include "TFile.h"
#include "TString.h"
#include "TStopwatch.h"
#define STEPS 1000

void test_cooling(){
  TStopwatch timer;
  SimulationLattice s(10, 2, 1, "test_cooling.root", STEPS, 2, 5, 1);
  timer.Start();
  s.run();
  timer.Stop();
  timer.Print();
  
  TFile file("test_cooling.root");

  double* x = new double[STEPS];
  double* y = new double[STEPS];

  Block * bb = *(Block**) file.Get("Block");
  
  timer.Start();
  for(unsigned int i = 0; i < STEPS + 1; i++){
    cout << i << '\n' << flush; 
    y[i] = bb[i].E;
    x[i] = i;
  }
  timer.Stop();
  timer.Print();
  TGraph* gr = new TGraph(STEPS, x, y);
  gr->Draw();

  delete[] x;
  delete[] y;
}


