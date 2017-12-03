//
//  sap.C
//  
//
//  Created by Federico on 03/12/17.
//

#include <stdio.h>
#include "SimulationLattice.h"
#include "AnalysisLattice.h"
#include "TF1.h"
#include "TMath.h"


{
/*
  SimulationLattice s(10, 2, 10, "provaM.root", 1000, 10000, 1.5, 3., 50);
  s.run();
  
  AnalysisLattice a("provaM.root", "outM.root");
  a.run();
  
  TGraphErrors * g = a.drawLattice(5, TEMPERATURE, MAGNETIZATION, true);
  g -> Draw();
  
  TF1 * f = new TF1("f", "[0] * TMath::Power( TMath::Abs( x - 2.27 ) , [1] )", 1.5, 2.27);
  //f -> Draw("same");
  g -> Fit(f);



  s = SimulationLattice(10, 2, 10, "provaM.root", 1000, 10000, 1.5, 3., 50);
  s.run();
  
  a = AnalysisLattice("provaM.root", "outM.root");
  a.run();
*/
  g = a.drawLattice(5, TEMPERATURE, MAGNETIZATION, true);
  g -> Draw();
  
  f = new TF1("f", "TMath::Power( TMath::Abs( x - 2.27 ) , [0] )", 1.5, 2.27);
  g -> Fit(f);
  //f -> Draw("same");
}
