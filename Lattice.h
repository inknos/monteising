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
  unsigned int num_spin;
  bool * lattice;

 public:
  Lattice();

  Lattice(const unsigned int& _N, const unsigned int& _dim, unsigned const int& _q);

  Lattice(const Lattice &obj);
    
  ~Lattice();
  
  Lattice& operator=(const Lattice& obj);
  
  int getDim();

  int getN();
  
  int getSpin(int i);

  int getQ();
  
  int getNumSpin();
  
  void printLattice();
  
  double energy(bool pauli=false);

};

#endif
