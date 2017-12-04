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
  
  void setN();

  void setN(unsigned int _N);
  
  void setDim();

  void setDim(unsigned int _dim);
  
  void setNumSpin();

  void setNumSpin(unsigned int _num_spin);

  void draw2D();

  void draw3D();
  
 public:
  DrawLattice();

  DrawLattice(const Lattice& lat);

  DrawLattice(const TString& fi, const TString& nf);
  
  Lattice getLattice() const;

  unsigned int getN() const;

  unsigned int getDim() const;

  unsigned int getNumSpin() const;
  
  void readFile(const TString& fname, const TString& lname);

  void draw();
  
  ClassDef(DrawLattice, 1)
  
};

#endif
