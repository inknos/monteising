#include "AnalysisLattice.h"

#include "TString.h"
#include "TCanvas.h"


void exportGraph(const TString& filename){
  AnalysisLattice a("in.root", filename);

  new TCanvas;
  a.drawLattice(0, TEMPERATURE, MAGNETIZATION) -> Draw("ALP");
  new TCanvas;
  a.fitLattice(false, TEMPERATURE, MAGNETIZATION) -> Draw();

  new TCanvas;

  new TCanvas;

  
}
