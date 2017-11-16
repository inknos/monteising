#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"

#include "DrawLattice.h"
#include "Lattice.h"

ClassImp(DrawLattice)

typedef unsigned int uint;

DrawLattice::DrawLattice() : lattice(), N(1), dim(1), num_spin(1),
			     fname(""), lname(""), cname("Default Canvas"),
			     ctitle("cv"), gname("Default Ising") {} 

/* Public Functions */
Lattice DrawLattice::getLattice() const { return lattice; }

uint DrawLattice::getN() const { return N; }

uint DrawLattice::getDim() const { return dim; }

uint DrawLattice::getNumSpin() const { return num_spin; }

void DrawLattice::readFile(const TString& fname, const TString& lname){
  TFile f(fname);
  lattice = * (Lattice*) f.Get(lname);
  f.Close();
  setN();
  setDim();
  setNumSpin();
}

/* Private Functions */
void DrawLattice::setN(){ N = lattice.getN(); }

void DrawLattice::setN(uint n){ N = n; }

void DrawLattice::setDim(){ dim = lattice.getDim(); }

void DrawLattice::setDim(uint d){ dim = d; }

void DrawLattice::setNumSpin(){ num_spin = lattice.getNumSpin(); }

void DrawLattice::setNumSpin(uint ns){ num_spin = ns; }

void DrawLattice::draw2D(){
  uint index = 0;
  uint zeros = 0;
  uint ones  = 0;

  for(uint i = 0; i < num_spin; i++){
    lattice.getSpin(i) ? ones++ : zeros++;
  }

  TCanvas* c1 = new TCanvas("c1","Ising Canvas",200,10,700,500);

  TMultiGraph *mg = new TMultiGraph();
   mg->SetTitle("Ising");
  
  TGraph * g0 = new TGraph(zeros);
  for(uint i = 0; i < num_spin; i++){
    if(!lattice.getSpin(i)){
      g0->SetPoint(index,
                   (int) 1 + (i % N),
                   (int) (N - i / N)
                   );
      index++;
    }
  }
  g0->SetMarkerStyle(21);
  g0->SetMarkerColor(kBlue);
  //g0->GetXaxis()->SetRangeUser(0., 100.);
  //g0->GetYaxis()->SetRangeUser(0., 100.);

  index = 0;
  TGraph * g1 = new TGraph(ones);
  for(uint i = 0; i < num_spin; i++){
    if(lattice.getSpin(i)){
      g1->SetPoint(index,
                   (int) 1 + (i % N),
                   (int) (N - i / N)
                   );
      index++;
    }
  }
  g1->SetMarkerStyle(21);
  g1->SetMarkerColor(kRed);
  //g1->GetXaxis()->SetRangeUser(0., 100.);
  //g1->GetYaxis()->SetRangeUser(0., 100.);

  mg->Add(g0);
  mg->Add(g1);
  mg->Draw("AP");
}

void DrawLattice::draw3D(){
  int index = 0;
  int po = (int) pow(N, 2);

  TCanvas* c1 = new TCanvas("c1","Ising Canvas",200,10,700,500);

  TGraph2D * g0 = new TGraph2D();
  for(uint i = 0; i < num_spin; i++){
    if(!lattice.getSpin(i)){
      g0->SetPoint(index,
                   (int) 1 + (i % N),
                   (int) N - (i % po) / N,
                   (int) N - (i / po)
                   );
      index++;
    }
  }
  g0->SetMarkerStyle(20);
  g0->SetMarkerSize(0.6);
  g0->SetMarkerColor(kBlue);
  g0->Draw("AP");


  index = 0;
  TGraph2D * g1 = new TGraph2D();
  for(uint i = 0; i < num_spin; i++){
    if(lattice.getSpin(i)){
      g1->SetPoint(index,
                   (int) 1 + (i % N),
                   (int) N - (i % po) / N,
                   (int) N - (i / po)
                   );
      index++;
    }
  }
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.6);
  g1->SetMarkerColor(kRed);
  g1->Draw("AP same");
}

void DrawLattice::draw(){
  void (DrawLattice::*func)();
  switch(dim){
  case 1: return;
  case 2: {
    func = &DrawLattice::draw2D;
    break;
  }
  case 3: {
    func = &DrawLattice::draw3D;
    break;  
  }
  default: return;
  }
  (*this.*func)();
}
