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


void sap_X(){
  SimulationLattice s(10, 2, 10, "prova.root", 1000, 10000, 1.5, 3., 20);
  s.run();
  
  AnalysisLattice a("prova.root", "out.root");
  a.run();

    
  TGraphErrors *g = a.drawLattice(3 ,  TEMPERATURE , MAGNETIZATION , true);
  //TGraphErrors *g = a.drawLatticeMean(TEMPERATURE, MAGNETIZATION);
  g -> Draw();
  

  TF1 * f = new TF1("f", "TMath::Power( TMath::Abs( x - 2.27 ) , [0] )", 1.5, 2.27);
  f -> Draw("same");
  g -> Fit(f);

/*

  g = a.drawLattice(3, TEMPERATURE, MAGNETIZATION, true);
  g = a.drawLatticeMean(TEMPERATURE, MAGNETIZATION);
  g -> Draw();
  
  f1 = new TF1("f", "TMath::Power( TMath::Abs( x - 2.27 ) , [0] )", 1.5, 3.5);

  g -> Fit(f);
  f -> Draw("same");

*/
  
}
