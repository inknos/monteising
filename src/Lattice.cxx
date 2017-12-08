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
#define K 1e-23
#define J 1
#define NOISE 1e-5

ClassImp(Lattice)

// setting static member Temperature
double Lattice::T = 0;

/* CONSTRUCTORS */

//Default
Lattice::Lattice() : TObject(), N(1), dim(1) , num_spin(1){
  lattice = new bool[1];
  lattice[0] = gRandom->Rndm() > 0.5 ? 0 : 1;
}

//Standard
Lattice::Lattice(const uint& _N, const uint& _dim) :
  TObject(), N(_N) , dim(_dim) , num_spin(pow(_N, _dim)){
  lattice = new bool[ num_spin ];
  for(uint i = 0; i < num_spin; i++){
    lattice[i] = gRandom->Rndm() > 0.5 ? 0 : 1;
  }
}

//Copy
Lattice::Lattice(const Lattice &obj) :
  TObject(), N(obj.N) , dim(obj.dim) , num_spin(obj.num_spin){
  lattice = new bool[obj.num_spin];
  for(uint i = 0; i < obj.num_spin; i++){
    lattice[i] = obj.lattice[i];
  }
}

/* DESTRUCTOR */

Lattice::~Lattice(){
  delete []lattice;
}



/* PHYSICAL AND NUMERICAL FUNCTIONS */

bool Lattice::flipSpin(const uint& n){
  if(n > num_spin) return false;
  lattice[n] = !lattice[n];
  return true;
}

int Lattice::dE(const uint& spin ) const {
  int dE_tmp = 0;
  if(spin > num_spin){ return 0; }
  uint pow_tmp1;
  uint pow_tmp2;
  uint i_tmp;
  for(uint d = 0; d < dim; d++){
    if(d == 0){
      pow_tmp1 = 1;
      pow_tmp2 = N;
    }
    else{
      pow_tmp1 = pow_tmp2;             
      pow_tmp2 = (uint) pow(N, d + 1); 
    }
    if(d == dim - 1){
      lattice[spin] ^ lattice[ (int) ( (spin + pow_tmp1)            % pow_tmp2 ) ] ? dE_tmp -= 1 : dE_tmp += 1;
      lattice[spin] ^ lattice[ (int) ( (spin - pow_tmp1 + pow_tmp2) % pow_tmp2 ) ] ? dE_tmp -= 1 : dE_tmp += 1;
    }
    else{
      i_tmp = ( (int) (spin / pow_tmp2) ) * pow_tmp2;
      lattice[spin] ^ lattice[ (int) ( i_tmp + (spin + pow_tmp1)            % pow_tmp2 ) ] ? dE_tmp -= 1 : dE_tmp += 1;
      lattice[spin] ^ lattice[ (int) ( i_tmp + (spin - pow_tmp1 + pow_tmp2) % pow_tmp2 ) ] ? dE_tmp -= 1 : dE_tmp += 1;
    }
  }
  return dE_tmp * 2;
}

int Lattice::energy() const {
  int E_tmp = 0;
  uint pow_tmp1 = 1;
  uint pow_tmp2 = N;
  uint i_tmp;
  for(uint d = 0; d < dim; d++) {
    for(uint i = 0; i < num_spin; i++) {
      if(d == dim - 1) {
        lattice[i] ^ lattice[ (int) ( (i + pow_tmp1) % pow_tmp2 ) ] ? E_tmp -= 1 : E_tmp += 1;
      }
      else{
        i_tmp = ( (int) (i / pow_tmp2) ) * pow_tmp2;
        lattice[i] ^ lattice[ (int) (i_tmp + (i + pow_tmp1) % pow_tmp2 ) ] ? E_tmp -= 1 : E_tmp += 1;
      }
    }
    if(d == dim - 1) continue;
    pow_tmp1 = pow_tmp2;             
    pow_tmp2 = (uint) pow(N, d + 2); 
  }
  return - E_tmp;
}

float Lattice::magnetization() const {
  float M_temp = 0;
  for(uint i = 0; i < num_spin; i ++) lattice[i] ? M_temp += 1: M_temp -=1;
  return M_temp / num_spin;
}

void Lattice::cooling(){
  double _T = gRandom -> Rndm() * NOISE + Lattice::T;
  uint spin = (uint) ( gRandom -> Rndm() * num_spin );
  int tmp_spin = dE(spin);
  if( tmp_spin < 0 ){
    flipSpin(spin);
  }
  else{
    if(gRandom -> Rndm() < TMath::Exp( - tmp_spin / _T)){
      flipSpin(spin);
    }
  }
}

void Lattice::cooling(const uint& iter) {
  double _T;
  uint spin;
  int tmp_spin;

  for(uint i = 0; i < iter; i++){
    _T = gRandom -> Rndm() * NOISE + Lattice::T;
    spin = (uint) ( gRandom -> Rndm() * num_spin );
    tmp_spin = dE(spin);
    if( tmp_spin < 0 ){
      flipSpin(spin);
    }
    else{
      if(gRandom -> Rndm() < TMath::Exp( - tmp_spin / _T)){
        flipSpin(spin);
      }
    }
  }
}

double * Lattice::coolingPar(){
  double * arr = new double[4];
  arr[0] = gRandom -> Rndm() * NOISE + Lattice::T;
  arr[1] = 0;
  arr[2] = 0;
  arr[3] = 0;
  uint spin = (uint) ( gRandom -> Rndm() * num_spin );
  int tmp_spin = dE(spin);
  if( tmp_spin < 0 ){
    flipSpin(spin);
    arr[1] = (double) tmp_spin;
    getSpin(spin) ? arr[2] = 2. / num_spin : arr[2] = - ( 2. / num_spin );
    arr[3] = (double) tmp_spin / num_spin;
  }
  else{
    if(gRandom -> Rndm() < TMath::Exp( - tmp_spin / arr[0])){
      flipSpin(spin);
      arr[1] = (double) tmp_spin;
      getSpin(spin) ? arr[2] = 2. / num_spin : arr[2] = - ( 2. / num_spin );
      arr[3] = (double) tmp_spin / num_spin;
    }
  }
  return arr;
}

/* OVERLOADED OPERATORS */

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

/* GETTERS AND SETTERS */

uint Lattice::getN() const{ return N; }

uint Lattice::getDim() const { return dim; }

uint Lattice::getNumSpin() const { return num_spin; }

bool Lattice::getSpin(const uint & i) const{
  if(i > num_spin) return i;
  return lattice[i];
}

double Lattice::getT() { return Lattice::T; }

void Lattice::setT(const double& _T){ Lattice::T = _T; }

/* OTHERS */  

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

void Lattice::printLatticeROOT(const TString& name, const TString& ln) const {
  TFile f(name + TString(".root"), "RECREATE");
  this -> Write(ln);
  f.Write();
  f.Close();
}


