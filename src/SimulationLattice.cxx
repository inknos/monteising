#include "SimulationLattice.h"

#include "TTree.h"
#include "TBranch.h"
#include "TBranchClones.h"
#include "TFile.h"
#include "TString.h"
#include "Lattice.h"
#include "TClonesArray.h"
#include "TDatime.h"

//#include "omp.h"
#include <string>
#include <iostream>

ClassImp(SimulationLattice)

SimulationLattice::SimulationLattice() :
N(1) , dim_vector(1) , dim(1) , file(""), I0(1), iter(1), tempmin(0.), tempmax(1.), tempstep(1.) {
  lattice_vector = new Lattice[1];
}

SimulationLattice::SimulationLattice(const uint& _N ,
                                     const uint& _dim ,
                                     const uint& _dim_vector
                                     ) :
  N(_N), dim(_dim), dim_vector(_dim_vector) {
  SimulationLattice(_N, dim, _dim_vector, "", 1, 1, 0., 1., 1.);
}

SimulationLattice::SimulationLattice(const uint& _N ,
                                     const uint& _dim ,
                                     const uint& _dim_vector,
                                     const TString& _file,
                                     const uint& _i0,
                                     const uint& _iter,
                                     const double& _tempmin,
                                     const double& _tempmax,
                                     const uint& _tempstep
                                     ) :
N(_N) , dim(_dim) , dim_vector(_dim_vector), file(_file), I0(_i0),     
  iter(_iter), tempmin(_tempmin), tempmax(_tempmax),
  tempstep(_tempstep)
{
  lattice_vector = new Lattice[_dim_vector];
  for(uint i = 0 ; i < _dim_vector ; i++) lattice_vector[i] = Lattice(_N , _dim );
}

SimulationLattice::SimulationLattice(const Lattice& _lat ,
                                     const uint & _dim_vector,
                                     const TString& _file,
                                     const uint& _i0,
                                     const uint& _iter,
                                     const double& _tempmin,
                                     const double& _tempmax,
                                     const uint& _tempstep) :
  N(_lat.getN()) , dim(_lat.getDim()) , dim_vector(_dim_vector)
{
  SimulationLattice(_lat.getN(), _lat.getDim(), _dim_vector,
                    _file, _i0, _iter, _tempmin, _tempmax, _tempstep);
}

SimulationLattice::SimulationLattice(const SimulationLattice& obj) :
  N(obj.N) , dim(obj.dim) , dim_vector(obj.dim_vector), file(obj.file), I0(obj.I0),
  iter(obj.iter), tempmin(obj.tempmin), tempmax(obj.tempmax), tempstep(obj.tempstep)
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

uint SimulationLattice::getI0() const { return I0; }

void SimulationLattice::setFile(const TString& f) { file =f; }

void SimulationLattice::setIter(const uint& i) { iter = i; }

void SimulationLattice::setTempMin(const double& t) { tempmin = t; }

void SimulationLattice::setTempMax(const double& t) { tempmax = t; }

void SimulationLattice::setTempStep(const uint& t) { tempstep = t; }

void SimulationLattice::setI0(const uint& _I0) { I0 = _I0; }

uint SimulationLattice::getT() { return Lattice::getT(); }

void SimulationLattice::setT(const double& _T) { Lattice::setT(_T); }

void SimulationLattice::run(){

  TDatime datime;
  unsigned int datime_t = datime.Get();
  const char * file_t = file;

  TFile f(file, "recreate");
  typedef struct{
    double _tempmin;
    double _tempmax;
    uint _N;
    uint _dim;
    uint _dim_vector;
    uint _I0;
    uint _iter;
    uint _tempstep;
    uint _datime_t;
  }Info;
  static Info info;
  
  info._N = N;
  info._dim = dim;
  info._dim_vector =  dim_vector;
  info._I0 = I0;
  info._iter = iter;
  info._tempstep = tempstep;
  info._datime_t = datime_t;
  info._tempmin = tempmin;
  info._tempmax = tempmax;
  
  cout << "info.tempmax : " << info._tempmax << endl;
  cout << "tempmax : " << tempmax << endl;
  
  TTree * tree = new TTree("info", "Info");
  tree -> Branch("Info", &info._tempmin, "_tempmin/D:_tempmax:_N/i:_dim:_dim_vector:_I0:_iter:_tempstep:_datime_t");
  tree -> Fill();


  TTree * tree = new TTree("info", "Info");
  tree -> Branch("N", N);
  tree -> Branch("dim", dim);
  tree -> Branch("dim_vector", dim_vector);
  tree -> Branch("I0", I0);
  tree -> Branch("iter", iter);
  tree -> Branch("tempmin", tempmin);
  tree -> Branch("tempmax", tempmax);
  tree -> Branch("tempstep", tempstep);
  tree -> Branch("datime", datime_t);
  
  TClonesArray * array = new TClonesArray("Block", iter + 1);

  int E_tmp;
  double M_tmp;
  double S_tmp;
  double T_tmp;
  double* data = new double[4];
  double tempN = (tempmax - tempmin) / (tempstep - 1);
 

  for(uint i = 0; i < dim_vector; i++){
    cout << "i: " << i << endl << flush;
    TString treeName(TString("Lattice_") + TString(std::to_string( i ).c_str() ));
    TString treeTitle(TString("LatticeTree_") + TString(std::to_string( i ).c_str() ));
    TTree * tree = new TTree(treeName, treeTitle);
    
    for(double t = tempmin; t <= tempmax; t += tempN){

      cout << "t: " << t << endl << flush;
      setT(t);
      tree -> Branch(TString("T_") + TString( TString(std::to_string( (int) ( (t-tempmin) / tempN ) ).c_str() ) ), &array);
  for(double t = tempmin; t <= tempmax; t += tempN ){
    setT(t);

    TString treeName(TString("T_") + TString(std::to_string( (int) ( (t-tempmin) / tempN ) ).c_str() ));
    TString treeTitle(TString("Tree_") + TString(std::to_string( (int) ( (t-tempmin) / tempN ) ).c_str() ));
    TTree * tree = new TTree(treeName, treeTitle);

    for(uint i = 0; i < dim_vector; i++){
      tree -> Branch(TString("Lattice_") + TString(std::to_string(i).c_str()), &array);

      lattice_vector[i].cooling(I0);
      E_tmp = lattice_vector[i].energy(false);
      M_tmp = lattice_vector[i].magnetization();
      S_tmp = ( (double) E_tmp )/ lattice_vector[i].getNumSpin();
      T_tmp = Lattice::getT();

      Block * b1 = (Block*)array -> ConstructedAt(0);
      b1 -> setBlock(i, T_tmp, E_tmp, M_tmp, S_tmp, I0);
      
      for(uint j = 0; j < iter; j++){
        data = lattice_vector[i].coolingPar();
        E_tmp += (int) data[1];
        M_tmp += data[2];
        S_tmp += data[3];
        T_tmp = data[0];

        Block * b2 = (Block*)array -> ConstructedAt(j + 1);
        b2 -> setBlock(i, T_tmp, E_tmp, M_tmp, S_tmp, j + 1);

      }
      tree -> GetBranch(TString("T_") + TString( TString(std::to_string( (int) ( (t-tempmin) / tempN ) ).c_str() ) ) ) -> Fill();
      }

      /*
      //DEBUG
      for (Int_t j=0; j< array->GetEntries(); j++){
	    Block *block=(Block*)array->At(j);
	    cout<< "E : " << block->E << " M : " << block->M << endl << flush;  
      }
      //DEBUG 
      */

    } 
    tree->Fill();
    std::cout << "Ho scritto il reticolo " << i << std::endl;
      std::cout << "Ho scritto il reticolo " << i << std::endl;
      tree->Fill();
      array->Clear();
    } 
    tree->Fill();
  }
  f.Write();
  f.Close();
}
