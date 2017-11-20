#include "SimulationLattice.h"

#include "TFile.h"
#include "TString.h"

#include <string>
#include <iostream>

ClassImp(SimulationLattice)

SimulationLattice::SimulationLattice() : dim_vector(1){
  lattice_vector = new Lattice[1];
}

SimulationLattice::SimulationLattice(const Lattice& lat , const uint & _dim_vector ) : dim_vector(_dim_vector){
  lattice_vector = new Lattice[_dim_vector];
  for(uint i = 0 ; i < _dim_vector ; i++) lattice_vector[i] = lat;
}

SimulationLattice::SimulationLattice(const uint& _N , const uint& _dim , const uint& _dim_vector) : dim_vector(_dim_vector){
  lattice_vector = new Lattice[_dim_vector];
  for(uint i = 0 ; i < _dim_vector ; i++) lattice_vector[i] = Lattice(_N , _dim );
}

SimulationLattice::SimulationLattice(const SimulationLattice& obj) : dim_vector(obj.dim_vector){
  lattice_vector = new Lattice[obj.dim_vector];
  for(uint i = 0 ; i < obj.dim_vector ; i++) lattice_vector[i] = obj.lattice_vector[i];
}

SimulationLattice& SimulationLattice::operator=(const SimulationLattice& obj){
  if(&obj == this) return *this;
  this -> ~SimulationLattice();
  new(this) SimulationLattice(obj);
  return *this;
}

SimulationLattice::~SimulationLattice(){
  delete []lattice_vector;
}

Lattice SimulationLattice::getLattice(const uint& i) const { return lattice_vector[i]; }

uint SimulationLattice::getDimVector() const { return dim_vector; }

uint SimulationLattice::getT() { return Lattice::getT(); }

void SimulationLattice::setT(const double& _T) { Lattice::setT(_T); }

void SimulationLattice::simulation(const TString& fname,
                                   const uint& iter,
                                   const double& tempmin = 0,
                                   const double& tempmax = 5,
                                   const uint& tempstep= 10){
  /* Ntot steps = t * dim_vector * iter */
  
  TFile f(fname, "RECREATE");
  std::string lv("lattice_vector");
  for(double t = tempmin; t <= tempmax; t += (tempmax - tempmin) / tempstep){
    Lattice::setT(t);
    //std::cout << "Cooling at T = " << t << std::endl << std::flush;
    for(uint i = 0; i < dim_vector; i++){
      lattice_vector[i].cooling(iter);
      lattice_vector[i].Write( (lv + std::to_string(i) + "_" + std::to_string(t)).c_str() );
      //std::cout << i << "/" << dim_vector << std::endl << std::flush;
    }
  }
  f.Write();
  f.Close();
}
