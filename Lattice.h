#ifndef LATTICE_H
#define LATTICE_H

#include "TRandom3.h"

#include <math.h>
#include <iostream>


class Lattice {
 private:
  unsigned int dim;
  unsigned int q;
  unsigned int N;
  bool * lattice;

 public:
  Lattice();

  Lattice(const unsigned int& _N, const unsigned int& _dim, unsigned const int& _q);

  Lattice(const Lattice &obj);
    
  ~Lattice();
  
  int getDim();

  int getN();
  
  int getSpin(int i);

  int getQ();
  
  void printLattice();
    

};

#endif
