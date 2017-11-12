#ifndef LATTICE_H
#define LATTICE_H

#include "TRandom3.h"

#include <math.h>
#include <iostream>


class Lattice {
 private:
  unsigned int dim;       // dimension of the lattice
  unsigned int q;         // 
  unsigned int N;         // number of spins in one dimension
  unsigned int num_spin;  // pow(N, dim) total number of spins
  bool * lattice;         // [num_spin] lattice

 public:
  /* Public Constructors */
  Lattice();

  Lattice(const unsigned int& _N, const unsigned int& _dim, unsigned const int& _q);

  Lattice(const Lattice &obj);

  ~Lattice();

  /* Public Operators */
  Lattice& operator=(const Lattice& obj);

  friend std::ostream &operator<<(std::ostream &out, const Lattice &lat);

  /* Getters */
  int getDim();

  int getN();

  int getSpin(unsigned int i);

  int getQ();

  int getNumSpin();

  /* Print - deprecated*/
  //void printLattice();

  /* Physical functions */
  int energy(bool pauli=false);

  ClassDef(Lattice,1)
};

#endif
