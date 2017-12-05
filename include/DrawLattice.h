#ifndef DRAW_LATTICE_H
#define DRAW_LATTICE_H

#include "Lattice.h"

#include "TString.h"
#include "TFile.h"

class DrawLattice {
 
 private:
 
  Lattice lattice;            //Lattice object    
  uint N , dim , num_spin;    //N , dim , num_spin of lattice
  
  TString fname;              //name of root file 
  TString lname;              //name of Lattice saved in file fname.root
  
  TString cname , ctitle;     //canvas name and canvas title
  TString gname , gtitle;     //graph name and graph title
  
  /* DRAW FUNCTIONS */

  void draw2D();

  void draw3D();
  
  /* SETTERS */
   
  void setN();

  void setN(uint _N);
  
  void setDim();

  void setDim(uint _dim);
  
  void setNumSpin();

  void setNumSpin(uint _num_spin);
  
  
  
 public:
 
  /* CONSTRUCTORS */
  
  DrawLattice();

  DrawLattice(const Lattice& lat);

  DrawLattice(const TString& _fname, const TString& _lname);
  
  void readFile(const TString& _fname, const TString& _lname);
  
  /* DRAW FUNCTION */
  
  void draw();
  
  /* GETTERS */
  
  Lattice getLattice() const;

  uint getN() const;

  uint getDim() const;

  uint getNumSpin() const;
  
  
  
  
  
  ClassDef(DrawLattice, 1)
  
};

#endif
