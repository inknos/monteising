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
#define STEPS 10000

void test_cooling(){
  TStopwatch timer;
  SimulationLattice s(10, 2, 1, "test_cooling.root", 1500, STEPS, 2., 5., 5);
  
  timer.Start();
  s.run();
  timer.Stop();
  timer.Print();
  
  TFile file("test_cooling.root");

  TClonesArray *hits = new TClonesArray("Block", STEPS + 1);

  TTree *tree = (TTree*) file.Get("tree");

  TBranchClones* branch = (TBranchClones*) ( tree -> GetBranch(TString("Lattice_") + TString(std::to_string(0).c_str()))  );

  branch -> SetAddress(&hits);
 
  tree -> GetEntry(0);
  
  double* x = new double[STEPS + 1];
  double* y = new double[STEPS + 1];
  
  for(uint  j = 0; j < STEPS +1; j++){
    
    y[j] = ((Block*) hits->At(j)) -> E;
    x[j] = j;
  }
  
  TGraph* gr = new TGraph(STEPS + 1, x, y);
  gr->Draw();

  delete[] x;
  delete[] y;
}


