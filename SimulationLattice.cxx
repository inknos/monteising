#include "SimulationLattice.h"



ClassImp(AnalysisLattice)


SimulationLattice::SimulationLattice() : dim_vector(1){
  lattice_vector = new Lattice[1]; 
}



SimulationLattice::SimulationLattice(const Lattice& lat , const uint & _dim_vector ) : dim_vector(_dim_vector){
  lattice_vector = new Lattice[_dim_vector]; 
  for(int i = 0 ; i < _dim_vector ; i++) lattice_vector[i] = lat; 
}

SimulationLattice::SimulationLattice(const uint& _N , const uint& _dim , const uint& _dim_vector) : dim_vector(_dim_vector){
  lattice_vector = new Lattice[_dim_vector]; 
  for(int i = 0 ; i < _dim_vector ; i++) lattice_vector[i] = Lattice(_N , _dim );
} 

SimulationLattice::SimulationLattice(const SimulationLattice& obj) : dim_vector(obj.dim_vector){
  lattice_vector = new Lattice[obj.dim_vector]; 
  for(int i = 0 ; i < obj.dim_vector ; i++) lattice_vector[i] = obj.lattice_vector[i];
}

SimulationLattice& SImulationLattice::operator=(const SimulationLattice& obj){
  if(&obj == this) return *this;
  this -> ~SimulationLattice();
  new(this) SimulationLattice(obj);
  return *this;
}






SimulationLattice::~SimulationLattice(){
  delete []lattice_vector;
}      





