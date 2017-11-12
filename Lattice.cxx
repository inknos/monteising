#include "TRandom3.h"

#include "Lattice.h"

#include <math.h>
#include <iostream>

typedef unsigned int uint;

Lattice::Lattice() : q(2), N(1), dim(1) {
  lattice = new bool[1];
  lattice[0] = gRandom->Rndm() > 0.5 ? 0 : 1;
}

Lattice::Lattice(const uint& _N, const uint& _dim, const uint& _q) :
  N(_N) , dim(_dim), q(_q) {
  lattice = new bool[ (int) pow(N, dim) ];
  for(uint i = 0; i < (uint) pow(N, dim); i++){
    lattice[i] = gRandom->Rndm() > 0.5 ? 0 : 1;
  }
}

Lattice::~Lattice(){
  delete []lattice;
}

int Lattice::getDim(){ return dim; }

int Lattice::getN(){ return N; }

int Lattice::getSpin(int i){
  if(i < 0) return i;
  if((uint) i > dim) return i;
  return lattice[i];
}

int Lattice::getQ(){ return q; }

void Lattice::printLattice(){
  for(uint i = 0; i < (uint) pow(N, dim); i++){
    std::cout << lattice[i] << " " << std::flush;
    //if( (i + 1 ) % N == 0) std::cout << std::endl << std::flush;
    for(uint j = 1; j < dim; j++){
      if( (i + 1 ) % (uint) pow(N, dim - j) == 0 ) std::cout << std::endl << std::flush;
    }
  }
}
