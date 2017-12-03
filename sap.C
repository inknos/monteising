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

/*
{
  SimulationLattice s(10, 2, 10, "prova1.root", 1000, 10000, 1.5, 3., 20);
  s.run();
  
  AnalysisLattice a("prova1.root", "out1.root");
  a.run();
  
  TGraphErrors * g = a.drawLattice(3, TEMPERATURE, MAGNETIZATION, true);
  g -> Draw();
  
  TF1 * f = new TF1("f", "TMath::Power( TMath::Abs( x - 2.27 ) , 0.125 )", 1.5, 2.27);
  f -> Draw("same");
  
}
*/

{
  s = SimulationLattice(10, 2, 10, "prova1.root", 1000, 10000, 1.5, 3., 13);
  s.run();
  
  a = AnalysisLattice("prova1.root", "out1.root");
  a.run();
  
  g = a.drawLattice(5, TEMPERATURE, MAGNETIZATION, true);
  g -> Draw();
  
  f = new TF1("f", "[0]*TMath::Power( TMath::Abs( x - 2.27 ) , [1] )", 1.5, 2.27);
  g -> Fit(f);
  f -> Draw("same");
  
  
}

