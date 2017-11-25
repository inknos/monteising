#include "SimulationLattice.h"

#include "TFile.h"
#include "TString.h"
#include "Lattice.h" 
#include <string>
#include <iostream>

ClassImp(SimulationLattice)

SimulationLattice::SimulationLattice() :
N(1) , dim_vector(1) , dim(1) , file(""), iter(1), tempmin(0.), tempmax(1.), tempstep(1.) {
  lattice_vector = new Lattice[1];
}

SimulationLattice::SimulationLattice(const uint& _N ,
                                     const uint& _dim ,
                                     const uint& _dim_vector,
                                     const TString& _file,
                                     const uint& _iter,
                                     const double& _tempmin,
                                     const double& _tempmax,
                                     const uint& _tempstep
                                     ) :
  N(_N) , dim(_dim) , dim_vector(_dim_vector), file(_file),     
  iter(_iter), tempmin(_tempmin), tempmax(_tempmax),
  tempstep(_tempstep)
{
  lattice_vector = new Lattice[_dim_vector];
  for(uint i = 0 ; i < _dim_vector ; i++) lattice_vector[i] = Lattice(_N , _dim );
}

SimulationLattice::SimulationLattice(const Lattice& _lat ,
                                     const uint & _dim_vector,
                                     const TString& _file,
                                     const uint& _iter,
                                     const double& _tempmin,
                                     const double& _tempmax,
                                     const uint& _tempstep) :
N(_lat.getN()) , dim(_lat.getDim()) , dim_vector(_dim_vector)
{
  SimulationLattice(_lat.getN(), _lat.getDim(), _dim_vector,
                    _file, _iter, _tempmin, _tempmax, _tempstep);
}

SimulationLattice::SimulationLattice(const SimulationLattice& obj) :
  N(obj.N) , dim(obj.dim) , dim_vector(obj.dim_vector), file(obj.file), iter(obj.iter),
  tempmin(obj.tempmin), tempmax(obj.tempmax), tempstep(obj.tempstep)
{
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

uint SimulationLattice::getN() const { return N; }

uint SimulationLattice::getDim() const { return dim; }

TString SimulationLattice::getFile() const { return file; }

uint SimulationLattice::getIter() const { return iter; }

double SimulationLattice::getTempMin() const { return tempmin; }

double SimulationLattice::getTempMax() const{ return tempmax; }

double SimulationLattice::getTempStep() const {return tempstep; }

void SimulationLattice::setFile(const TString& f) { file =f; }

void SimulationLattice::setIter(const uint& i) { iter = i; }

void SimulationLattice::setTempMin(const double& t) { tempmin = t; }

void SimulationLattice::setTempMax(const double& t) { tempmax = t; }

void SimulationLattice::setTempStep(const uint& t) { tempstep = t; }

uint SimulationLattice::getT() { return Lattice::getT(); }

void SimulationLattice::setT(const double& _T) { Lattice::setT(_T); }



void SimulationLattice::run(){
  
  TFile f(file, "recreate");
  Block* block_vector = new Block[iter];
  uint I_0 = 100;
  int E_tmp;
  double M_tmp;
  double S_tmp;
  double T_tmp;
  setT(tempmin);
  double* data = new double[4];
 
 
  for(uint i = 0; i < dim_vector; i++){
    
    lattice_vector[i].cooling(I_0);
    E_tmp = lattice_vector[i].energy(false);
    M_tmp = lattice_vector[i].magnetization();
    S_tmp = ( (double) E_tmp )/ lattice_vector[i].getNumSpin();
    T_tmp = Lattice::getT();
    
    block_vector[0] = Block(i, T_tmp, E_tmp, M_tmp, S_tmp, I_0);
    block_vector[0].Write( (TString)std::to_string(i) );
    
    for(uint j = 0; j < iter; j++){
      
      data = lattice_vector[i].coolingPar();
      E_tmp += data[1];
      M_tmp += data[2];
      S_tmp += data[3];
      T_tmp = data[0];
      
      block_vector[j] = Block(i, T_tmp, E_tmp, M_tmp, S_tmp, j);
    }
    
    for(uint k = 0; k < iter ; k++){
      block_vector[k].Write( std::to_string(i).c_str() );
    }
    std::cout << "Ho scritto il reticolo " << i << std::endl;
  }
  f.Write();
  f.Close();
}



/*
void SimulationLattice::simulation(const TString& fname,
                                   const uint& iter,
                                   const double& tempmin = 0,
                                   const double& tempmax = 5,
                                   const uint& tempstep= 10){
  // Ntot steps = t * dim_vector * iter
  
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
*/
