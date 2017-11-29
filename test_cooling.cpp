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

void test_cooling(int n = 10, int dim = 2, int dim_v = 3, const char* nn = "test_cooling.root",
                  int i0 = 100, int iter = 1000, double t1 = 0., double t2 = 2., int step = 2 ){
  TStopwatch timer;
  
  SimulationLattice s(n, dim, dim_v, nn, i0, iter, t1, t2, step);

  timer.Start();
  s.run();
  timer.Stop();
  timer.Print();
  
  TFile file("test_cooling.root");

  TClonesArray *hits = new TClonesArray("Block", iter + 1);

  TTree *tree = (TTree*) file.Get("Lattice_0");

  TBranchClones* branch = (TBranchClones*) ( tree -> GetBranch(TString("T_") + TString(std::to_string(1).c_str()))  );

  branch -> SetAddress(&hits);

  tree -> GetEntry(0);

  double* x = new double[iter + 1];
  double* y = new double[iter + 1];

  for(uint  j = 0; j < iter +1; j++){

    y[j] = ((Block*) hits->At(j)) -> M;
    x[j] = j;
  }
  
  timer.Start();
  TGraph* gr = new TGraph(iter + 1, x, y);
  gr->Draw();
  timer.Stop();
  timer.Print();

  delete[] x;
  delete[] y;

  
  AnalysisLattice a(nn, "out.root");
  timer.Start();
  a.run();
  timer.Stop();
  timer.Print();
  
}
