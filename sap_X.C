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
  //SimulationLattice s(10, 2, 10, "provaX.root", 1000, 10000, 1.5, 3., 20);
  //s.run();
  
  a = AnalysisLattice("prova.root", "out.root");
  a.run();
  
  g = a.drawLatticeMean(TEMPERATURE, MAGNETIZATION);
  g -> Draw();
  

  f = new TF1("f", "TMath::Power( TMath::Abs( x - 2.30 ) , [0] )", 1.5, 2.3);
  //f -> Draw("same");
  g -> Fit(f);




/*
  s = SimulationLattice(10, 2, 10, "provaX.root", 1000, 10000, 1.5, 3., 50);
  s.run();
  
  a = AnalysisLattice("provaX.root", "outX.root");
  a.run();

  g = a.drawLattice(3, TEMPERATURE, MAGNETIZATION, true);
  g -> Draw();
  
  f = new TF1("f", "TMath::Power( TMath::Abs( x - 2.27 ) , [0] )", 1.5, 3.5);
  g -> Fit(f);
  //f -> Draw("same");
*/
}
