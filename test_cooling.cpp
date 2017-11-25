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


void test_cooling(){
  SimulationLattice s(100, 2, 1, "test_cooling.root", 10000, 0, 5, 1);
  s.run();
  TFile file("test_cooling.root");
  double* x = new double[10000];
  double* y = new double[10000];
  
  for(unsigned int i = 0; i < 10000; i++){
    y[i] = ((Block*) file.Get( (TString)(string("0;") + std::to_string(i+2)) ))->E;
    x[i] = i;
  }
  TGraph* gr = new TGraph(10000, x, y);
  gr->Draw();
  
  delete[] x;
  delete[] y;
}
