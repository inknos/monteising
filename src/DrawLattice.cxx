#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include "TLegend.h"

#include "DrawLattice.h"

#define COLOR_ZEROS kBlue
#define COLOR_ONES  kRed

#define MARKER_STYLE_2D 21

#define MARKER_STYLE_3D 20

#define MARKER_SIZE_3D 0.01

ClassImp(DrawLattice)

//Default
DrawLattice::DrawLattice() : 
           lattice(), 
           N(1), dim(1), num_spin(1),
			     fname(""), lname(""), 
			     cname("cv"), ctitle("Default Canvas"), 
			     gname("gr"), gtitle("Ising") {} 

//Standard1 - Draw from Lattice
DrawLattice::DrawLattice(const Lattice& lat) :
           lattice(lat), 
           N(lat.getN()), dim(lat.getDim()), num_spin(lat.getNumSpin()),
           fname(""), lname(""), 
           cname("cv"), ctitle("default canvas"), 
           gname("gr"), gtitle("Ising") {}

//Standard2 - Draw from File
DrawLattice::DrawLattice(const TString& _fname, const TString& _lname) :
           fname(_fname), lname(_lname), 
           cname("cv"), ctitle("default canvas"),
           gname("gr"), gtitle("Ising") {readFile(fname, lname);}

// readFile for Standard2
void DrawLattice::readFile(const TString& fname, const TString& lname){
  TFile f(fname);
  lattice = * (Lattice*) f.Get(lname);
  f.Close();
  setN();
  setDim();
  setNumSpin();
}


/* DRAW FUNCTION */ 

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

/* GETTERS */ 

Lattice DrawLattice::getLattice() const { return lattice; }

uint DrawLattice::getN() const { return N; }

uint DrawLattice::getDim() const { return dim; }

uint DrawLattice::getNumSpin() const { return num_spin; }


/* Private Functions */

/* DRAW FUNCTIONS */

void DrawLattice::draw2D(){
  uint index = 0;
  uint zeros = 0;
  uint ones  = 0;

  for(uint i = 0; i < num_spin; i++){
    lattice.getSpin(i) ? ones++ : zeros++;
  }

  TCanvas* c1 = new TCanvas(cname,ctitle,200,10,700,500);

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(gtitle);
  
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
  g0->SetMarkerStyle(MARKER_STYLE_2D);
  g0->SetMarkerColor(COLOR_ZEROS);
  
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
  g1->SetMarkerStyle(MARKER_STYLE_2D);
  g1->SetMarkerColor(COLOR_ONES);
  
  mg->Add(g0);
  mg->Add(g1);
  mg->SetTitle("2D Ising lattice");
  mg->Draw("AP");
}

void DrawLattice::draw3D(){
  int index = 0;
  int po = (int) pow(N, 2);

  TCanvas* c1 = new TCanvas(cname, ctitle,200,10,700,500);

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
  g0->SetTitle("3D Ising lattice");
  g0->SetMarkerStyle(MARKER_STYLE_3D);
  g0->SetMarkerSize(MARKER_SIZE_3D);
  g0->SetMarkerColor(COLOR_ZEROS);
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
  g1->SetTitle("3D Ising lattice");
  g1->SetMarkerStyle(MARKER_STYLE_3D);
  g1->SetMarkerSize(MARKER_SIZE_3D);
  g1->SetMarkerColor(COLOR_ONES);
  g1->Draw("AP same");
}

/* SETTERS */

void DrawLattice::setN(){ N = lattice.getN(); }

void DrawLattice::setN(uint _N){ N = _N; }

void DrawLattice::setDim(){ dim = lattice.getDim(); }

void DrawLattice::setDim(uint _dim){ dim = _dim; }

void DrawLattice::setNumSpin(){ num_spin = lattice.getNumSpin(); }

void DrawLattice::setNumSpin(uint _num_spin){ num_spin = _num_spin; }

