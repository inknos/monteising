#ifndef DRAW_LATTICE_H
#define DRAW_LATTICE_H

#include "Lattice.h"

#include "TString.h"
#include "TFile.h"

class DrawLattice {
 private:
  Lattice lattice;
  unsigned int N;
  unsigned int dim;
  unsigned int num_spin;

  TString fname;
  TString lname;
  
  TString cname;
  TString ctitle;

  TString gname;
  TString gtitle;

  TFile   file;
  
  void setN();

  void setN(unsigned int);
  
  void setDim();

  void setDim(unsigned int);
  
  void setNumSpin();

  void setNumSpin(unsigned int);

  void draw2D();

  void draw3D();
  
 public:
  DrawLattice();

  DrawLattice(const Lattice&);

  DrawLattice(const TFile&, const TString&);
  
  Lattice getLattice() const;

  unsigned int getN() const;

  unsigned int getDim() const;

  unsigned int getNumSpin() const;
  
  void readFile(const TString&, const TString&);

  void draw();
  
  ClassDef(DrawLattice, 1)
  
};

#endif
