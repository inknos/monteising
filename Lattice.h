#ifndef LATTICE_H
#define LATTICE_H

#include "TRandom3.h"

#include <math.h>
#include <iostream>


class Lattice {
 private:
  const unsigned int dim;       // dimension of the lattice
  const unsigned int q;         // 
  const unsigned int N;         // number of spins in one dimension
  const unsigned int num_spin;  // pow(N, dim) total number of spins
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
  int getDim() const;

  int getN() const;

  int getSpin(const unsigned int & i) const;

  int getQ() const;

  int getNumSpin() const;

  /* Print - deprecated*/
  //void printLattice();

  /* Physical functions */
  int energy() const;

  int energyParallel(int) const;
  
  ClassDef(Lattice,1)
};

#endif
