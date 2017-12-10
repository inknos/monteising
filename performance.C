//
//  performance.C
//  
//
//  Created by Federico on 10/12/17.
//

#include <stdio.h>
#include "AnalysisLattice.h"
#include "SimulationLattice.h"
#include "TStopwatch.h"
#include "TString.h"

#define unsigned int uint;

void performanceSim(uint size, uint N, uint I0, uint I, uint tempSteps);
void performance(uint par);

void performanceSim(uint size=10, uint N=1, uint I0=1000, uint I=1000, uint tempSteps=4){
  TStopwatch timer;
  timer.Start();
  SimulationLattice s(size, 2, N, "performance.root", I0, I, 0., 3.5, tempSteps);
  s.run();
  timer.Stop();
  timer.Print();
  
  timer.Start();
  AnalysisLattice a("performance.root", "performanceOut.root");
  a.run();
  timer.Stop();
  timer.Print();
  
}

void performance(uint par){
  std::cout << "\t \t \t \t" << "PAR : " << par << endl;
  if(par == 1){
    performanceSim();
    performanceSim(2*10);
    performanceSim(4*10);
  }
  if(par == 2){
    performanceSim();
    performanceSim(10, 2*1);
    performanceSim(10, 4*1);
  }
  if(par == 3){
    performanceSim();
    performanceSim(10, 1, 2*1000);
    performanceSim(10, 1, 4*1000);
  }
  if(par == 4){
    performanceSim();
    performanceSim(10, 1, 1000, 2*1000);
    performanceSim(10, 1, 1000, 4*1000);
  }
  if(par == 5){
    performanceSim();
    performanceSim(10, 1, 1000, 1000, 2*4);
    performanceSim(10, 1, 1000, 1000, 4*4);
  }
}
























