#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"

#include "Lattice.h"

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

#define H 1
#define K 1
#define J 1

//#include <omp.h>

typedef unsigned int uint;

ClassImp(Lattice)

Lattice::Lattice() : TObject(), q(2), N(1), dim(1) , num_spin(1){
  lattice = new bool[1];
  lattice[0] = gRandom->Rndm() > 0.5 ? 0 : 1;
}

Lattice::Lattice(const uint& _N, const uint& _dim, const uint& _q = 2) :
  TObject(), N(_N) , dim(_dim) , q(_q) , num_spin(pow(_N, _dim)){
  lattice = new bool[ num_spin ];
  for(uint i = 0; i < num_spin; i++){
    lattice[i] = gRandom->Rndm() > 0.5 ? 0 : 1;
  }
}

/* Copy-constructor */
Lattice::Lattice(const Lattice &obj) :
  TObject(), N(obj.N) , dim(obj.dim) , q(obj.q) , num_spin(obj.num_spin){
  lattice = new bool[obj.num_spin];
  for(uint i = 0; i < obj.num_spin; i++){
    lattice[i] = obj.lattice[i];
  }
}

Lattice::~Lattice(){
  delete []lattice;
}

/* Overload */
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

bool Lattice::operator==(const Lattice& obj){
  if(N != obj.N || dim != obj.dim) return false;
  for(uint i = 0; i < num_spin; i++){
    if(lattice[i] != obj.lattice[i]) return false;
  }
  return true;
}

/* Getters */
int Lattice::getDim() const { return dim; }

int Lattice::getN() const{ return N; }

bool Lattice::getSpin(const uint & i) const{
  //if(i < 0) return i;
  if(i > num_spin) return i;
  return lattice[i];
}

int Lattice::getQ() const { return q; }

int Lattice::getNumSpin() const { return num_spin; }

bool Lattice::flipSpin(const uint& n){
  if(n > num_spin) return false;
  //std::cout << "I'm fliping" << std::endl << std::flush;
  lattice[n] = !lattice[n];
  return true;
}

void Lattice::printLatticeCSV(const TString& name) const {
  std::ofstream file;
  file.open(name + TString(".csv"), std::ofstream::out | std::ofstream::trunc);
  file <<"#," << "N:" << N  << ", " << "dim:" << dim << ", " << "\n";
  for(uint i = 0; i < num_spin; i++){
    if(lattice[i]) file << 1 << ", " << std::flush;
    else file << -1 << ", " << std::flush;
    for(uint j = 1; j < dim; j++){
      if( (i + 1 ) % (uint) pow(N, dim - j) == 0 ) file << std::endl << std::flush;
    }
  }
  file.close();
}

void Lattice::printLatticeROOT(const TString& name, const TString& ln = "lat") const {
  TFile f(name + TString(".root"), "RECREATE");
  this -> Write(ln);
  f.Write();
  f.Close();
}

int Lattice::energy(const bool& p = false) const {
  int E_tmp = 0;
  uint pow_tmp1;
  uint pow_tmp2;
  uint i_tmp;
  for(uint i = 0; i < num_spin; i++){
    for(uint d = 0; d < dim; d++){
      //
      if(d == 0){ // skips first cycle pow and prevents one pow in every future cycle
        pow_tmp1 = 1;
        pow_tmp2 = N;
      }
      else{
        pow_tmp1 = pow_tmp2;             // this one is prevented by if statement
        pow_tmp2 = (uint) pow(N, d + 1); // one single pow for each loop
      }
      //
      if(d == dim - 1){ // last loop => biggest torus => no need to control position
        lattice[i] ^ lattice[ (int) ( (i + pow_tmp1)            % pow_tmp2 ) ] ? E_tmp -= 1 : E_tmp += 1;
        lattice[i] ^ lattice[ (int) ( (i - pow_tmp1 + pow_tmp2) % pow_tmp2 ) ] ? E_tmp -= 1 : E_tmp += 1;
      }
      else{
        i_tmp = ( (int) (i / pow_tmp2) ) * pow_tmp2;
        lattice[i] ^ lattice[ (int) ( i_tmp + (i + pow_tmp1)            % pow_tmp2 ) ] ? E_tmp -= 1 : E_tmp += 1;
        lattice[i] ^ lattice[ (int) ( i_tmp + (i - pow_tmp1 + pow_tmp2) % pow_tmp2 ) ] ? E_tmp -= 1 : E_tmp += 1;
      }
      /* Pauli term */
      if(p){ lattice[i] ? E_tmp += H : E_tmp -= H; }
      //
    }
  }
  return - E_tmp * 0.5;
}

/*
  int Lattice::energyParallel(int nt = 4) const {
  int E_tmp = 0;
  uint i(0);
  uint d(0);
  uint pow_tmp1(0);
  uint pow_tmp2(0);
  uint i_tmp(0);

  uint num_spin_t = num_spin;
  uint dim_t = dim;
  uint N_t = N;
  uint e = 0;

  #pragma omp parallel for private(pow_tmp1, pow_tmp2, i_tmp, d) shared(E_tmp) firstprivate(num_spin_t, dim_t, N_t, e)
  for(i = 0; i < num_spin_t; i++){
  e = 0;
  for(d = 0; d < dim_t; d++){
  pow_tmp1 = (uint) pow(N_t, d);
  pow_tmp2= (uint) pow(N_t, d + 1);
  i_tmp = ( (int) (i / pow_tmp2) - 1 ) * pow_tmp2;
  lattice[i] ^ lattice[(int) ( i_tmp + (i + pow_tmp1)            % pow_tmp2 )] ? e -= 1 : e += 1;
  lattice[i] ^ lattice[(int) ( i_tmp + (i - pow_tmp1 + pow_tmp2) % pow_tmp2 )] ? e -= 1 : e += 1;
  }
  E_tmp += e;
  // if collapse no code here
  }
  return E_tmp * 0.5;
  }
*/


int Lattice::dE(const uint& spin, const bool& p = false) const {
  int dE_tmp = 0;
  if(spin > num_spin){ return 0; }
  uint pow_tmp1;
  uint pow_tmp2;
  uint i_tmp;
  for(uint d = 0; d < dim; d++){
    //
    if(d == 0){ // skips first cycle pow and prevents one pow in every future cycle
      pow_tmp1 = 1;
      pow_tmp2 = N;
    }
    else{
      pow_tmp1 = pow_tmp2;             // this one is prevented by if statement
      pow_tmp2 = (uint) pow(N, d + 1); // one single pow for each loop
    }
    //
    if(d == dim - 1){ // last loop => biggest torus => no need to control position
      lattice[spin] ^ lattice[ (int) ( (spin + pow_tmp1)            % pow_tmp2 ) ] ? dE_tmp -= 1 : dE_tmp += 1;
      lattice[spin] ^ lattice[ (int) ( (spin - pow_tmp1 + pow_tmp2) % pow_tmp2 ) ] ? dE_tmp -= 1 : dE_tmp += 1;
    }
    else{
      i_tmp = ( (int) (spin / pow_tmp2) ) * pow_tmp2;
      lattice[spin] ^ lattice[ (int) ( i_tmp + (spin + pow_tmp1)            % pow_tmp2 ) ] ? dE_tmp -= 1 : dE_tmp += 1;
      lattice[spin] ^ lattice[ (int) ( i_tmp + (spin - pow_tmp1 + pow_tmp2) % pow_tmp2 ) ] ? dE_tmp -= 1 : dE_tmp += 1;
    }
    //if(p){ lattice[spin] ? E_tmp += H : E_tmp -= H; }
    //
  }
  return dE_tmp * 2;
}

float Lattice::magnetization() const {
  float M_temp = 0;
  for(uint i = 0; i < num_spin; i ++) lattice[i] ? M_temp += 1: M_temp -=1;
  return M_temp / num_spin;
}

/* Cooling */
void Lattice::cooling(const float& t = 0.5){
  double T = gRandom -> Rndm() * 1e-5 + t;
  uint spin = (uint) ( gRandom -> Rndm() * num_spin );
  int tmp_spin = dE(spin);
  if( tmp_spin < 0 ){
    flipSpin(spin);
  }
  else{
    if(gRandom -> Rndm() < TMath::Exp( - tmp_spin / T)){
      flipSpin(spin);
    }
  }
}
