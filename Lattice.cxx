#include "TRandom3.h"

#include "Lattice.h"

#include <math.h>
#include <iostream>

#include <omp.h>

typedef unsigned int uint;

ClassImp(Lattice)

Lattice::Lattice() : q(2), N(1), dim(1) , num_spin(1){
  lattice = new bool[1];
  lattice[0] = gRandom->Rndm() > 0.5 ? 0 : 1;
}

Lattice::Lattice(const uint& _N, const uint& _dim, const uint& _q = 2) :
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

std::ostream& operator<<(std::ostream &out, const Lattice& lat) {
  for(uint i = 0; i < lat.num_spin; i++){
    out << lat.lattice[i] << " " << std::flush;
    for(uint j = 1; j < lat.dim; j++){
      if( (i + 1 ) % (uint) pow(lat.N, lat.dim - j) == 0 ) out << std::endl << std::flush;
    }
  }
  return out;
}

int Lattice::getDim() const { return dim; }

int Lattice::getN() const{ return N; }

int Lattice::getSpin(const uint & i) const{
  //if(i < 0) return i;
  if(i > dim) return i;
  return lattice[i];
}

int Lattice::getQ() const { return q; }

int Lattice::getNumSpin() const { return num_spin; }

/*
  void Lattice::printLattice(){
  for(uint i = 0; i < num_spin; i++){
  std::cout << lattice[i] << " " << std::flush;
  //if( (i + 1 ) % N == 0) std::cout << std::endl << std::flush;
  for(uint j = 1; j < dim; j++){
  if( (i + 1 ) % (uint) pow(N, dim - j) == 0 ) std::cout << std::endl << std::flush;
  }
  }
  }
*/

int Lattice::energy() const {
  int E_tmp = 0;
  uint pow_tmp1;
  uint pow_tmp2;
  uint i_tmp;
  for(uint i = 0; i < num_spin; i++){
    for(uint d = 0; d < dim; d++){
      pow_tmp1 = (uint) pow(N, d);
      pow_tmp2= (uint) pow(N, d + 1);
      i_tmp = (int) (i / pow_tmp2);
      lattice[i] ^ lattice[(int) ( i_tmp * pow_tmp2 + (i + pow_tmp1)            % pow_tmp2 ) %num_spin] ? E_tmp -= 1 : E_tmp += 1;
      lattice[i] ^ lattice[(int) ( i_tmp * pow_tmp2 + (i - pow_tmp1 + num_spin) % pow_tmp2 ) %num_spin] ? E_tmp -= 1 : E_tmp += 1;
    }
  }
  E_tmp *= 0.5;
  return E_tmp;
}

int Lattice::energyParallel(int nt = 4) const {
  int E_tmp = 0;
  uint i(0);
  uint d(0);
  uint pow_tmp1(0);
  uint pow_tmp2(0);
  uint i_tmp(0);

#pragma omp parallel num_threads(nt)
  uint num_spin_t = num_spin;
#pragma omp for private(pow_tmp1, pow_tmp2, i_tmp, d)
  for(i = 0; i < num_spin; i++){
    uint e = 0;
    for(d = 0; d < dim; d++){
      pow_tmp1 = (uint) pow(N, d);
      pow_tmp2= (uint) pow(N, d + 1);
      i_tmp = (int) (i / pow_tmp2);
      lattice[i] ^ lattice[(int) ( i_tmp * pow_tmp2 + (i + pow_tmp1)            % pow_tmp2 ) %num_spin] ? e -= 1 : e += 1;
      lattice[i] ^ lattice[(int) ( i_tmp * pow_tmp2 + (i - pow_tmp1 + num_spin) % pow_tmp2 ) %num_spin] ? e -= 1 : e += 1;
    }
    E_tmp += e;
    // if collapse no code here
  }
  return E_tmp * 0.5;

}
