#include "TRandom3.h"

#include "Lattice.h"

#include <math.h>
#include <iostream>

typedef unsigned int uint;

Lattice::Lattice() : q(2), N(1), dim(1) , num_spin(1){
  lattice = new bool[1];
  lattice[0] = gRandom->Rndm() > 0.5 ? 0 : 1;
}

Lattice::Lattice(const uint& _N, const uint& _dim, const uint& _q) :
  N(_N) , dim(_dim) , q(_q) , num_spin(pow(_N, _dim)){
  lattice = new bool[ num_spin ];
  for(uint i = 0; i < num_spin; i++){
    lattice[i] = gRandom->Rndm() > 0.5 ? 0 : 1;
  }
}

Lattice::Lattice(const Lattice &obj) :
N(obj.N) , dim(obj.dim) , q(obj.q) , num_spin(obj.num_spin){
  lattice = new bool[obj.dim];
  for(uint i = 0; i < obj.dim ; i++){
    lattice[i] = obj.lattice[i];
  }
}

Lattice::~Lattice(){
  delete []lattice;
}

Lattice& Lattice::operator=(const Lattice& obj){
  if(&obj == this) return *this;
  this -> ~Lattice();
  new(this) Lattice(obj);
  return *this;
}

int Lattice::getDim(){ return dim; }

int Lattice::getN(){ return N; }

int Lattice::getSpin(int i){
  if(i < 0) return i;
  if((uint) i > dim) return i;
  return lattice[i];
}

int Lattice::getQ(){ return q; }

int Lattice::getNumSpin(){ return num_spin; }

void Lattice::printLattice(){
  for(uint i = 0; i < num_spin; i++){
    std::cout << lattice[i] << " " << std::flush;
    //if( (i + 1 ) % N == 0) std::cout << std::endl << std::flush;
    for(uint j = 1; j < dim; j++){
      if( (i + 1 ) % (uint) pow(N, dim - j) == 0 ) std::cout << std::endl << std::flush;
    }
  }
}

double Lattice::energy(bool){
  double E_tmp = 0;
  double pow_tmp;
  for(uint i = 0; i < num_spin; i++){
    for(uint d = 0; d < dim; d++){
      pow_tmp = pow(N, d);
      lattice[i] ^ lattice[(int) (i + pow_tmp)            %num_spin] ? E_tmp -= 1 : E_tmp += 1;
      lattice[i] ^ lattice[(int) (i - pow_tmp + num_spin) %num_spin] ? E_tmp -= 1 : E_tmp += 1;
    }
  }
  E_tmp *= 0.5;
  return E_tmp;
}
























