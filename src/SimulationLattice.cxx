#include "SimulationLattice.h"
#include "INFO.h"

#include "TTree.h"
#include "TBranch.h"
#include "TBranchClones.h"
#include "TFile.h"
#include "TString.h"
#include "Lattice.h"
#include "TDatime.h"
#include "TH1D.h"

#include "TDirectory.h"

//#include "omp.h"
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <utility>

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

uint SimulationLattice::getI0() const { return I0; }

double SimulationLattice::getTempMin() const { return tempmin; }

double SimulationLattice::getTempMax() const{ return tempmax; }

double SimulationLattice::getTempStep() const {return tempstep; }

uint SimulationLattice::getT() { return Lattice::getT(); }

void SimulationLattice::setFile(const TString& f) { file =f; }

void SimulationLattice::setIter(const uint& i) { iter = i; }

void SimulationLattice::setI0(const uint& _I0) { I0 = _I0; }

void SimulationLattice::setTempMin(const double& _tempmin) { tempmin = _tempmin; }

void SimulationLattice::setTempMax(const double& _tempmax) { tempmax = _tempmax; }

void SimulationLattice::setTempStep(const uint& _tempstep) { tempstep = _tempstep; }

void SimulationLattice::setT(const double& _T) { Lattice::setT(_T); }

void SimulationLattice::run(){

  std::cout << R"(
?=======================================?
!                                       !
!          Simulation Started           !
!                                       !
?=======================================?
!                                       !)";
  std::cout << endl
            << "! " << dim_vector << " lattices " << dim << "x" << dim << " with size " << N << " created\t!\n"
            << "! " << I0 << " iter. not collecting data\t!\n"
            << "! " << iter << " iter. collecting data  \t!\n"
            << "! from " << tempmin << " to " << tempmax << " in " << tempstep << "steps\t\t!\n"
            << R"(!                                       !
?=======================================?)" 
           << endl;
  
  TDatime datime;
  unsigned int datime_t = datime.Get();
  const char * file_t = file;

  TFile f(file, "recreate");
  static INFO info;

  info._N = N;
  info._dim = dim;
  info._dim_vector =  dim_vector;
  info._I0 = I0;
  info._iter = iter;
  info._tempstep = tempstep;
  info._datime_t = datime_t;
  info._tempmin = tempmin;
  info._tempmax = tempmax;

  TTree * tree = new TTree("info", "Info");
  tree -> Branch("Info", &info._tempmin, "_tempmin/D:_tempmax:_N/i:_dim:_dim_vector:_I0:_iter:_tempstep:_datime_t");
  tree -> Fill();
  tree -> Write();
  f.Write();
  tree -> Delete();
  Block * block = new Block[dim_vector];
  
  int    * E_tmp = new int[dim_vector];
  double * M_tmp = new double[dim_vector];
  double * S_tmp = new double[dim_vector];
  double * T_tmp = new double[dim_vector];
  double* data = new double[4];
  double tempN = (tempmax - tempmin) / (tempstep - 1);

  std::vector<double> temp_array;
  for(uint i = 0; i < tempstep; i++){
    temp_array.push_back( 0.);
  }
  double tc = 2.27 + gRandom -> Gaus(0, 0.1);
  double vc = 1;
  for(uint t = 0; t < tempstep / 4; t++){
    temp_array[ t * 4 + 0 ] =  tempmin + ( tc - vc - tempmin ) * ( (double) t / ((double)tempstep / 4) );
    temp_array[ t * 4 + 1 ] =  tc - vc + vc  * ( (double) t / ((double)tempstep / 4) );
    temp_array[ t * 4 + 2 ] =  tc + vc * ( (double) t / ((double)tempstep / 4) );
    temp_array[ t * 4 + 3 ] =  tc + vc + (tempmax - tc - vc) * ( (double) t / ((double)tempstep / 4) );
  }
  std::sort( temp_array.begin(), temp_array.end() );
  //for(uint i = 0; i < tempstep; i++){
  //  cout << "temp_array is : " << temp_array[i] << endl << flush;
  //}
  
  for(uint t = 0; t < tempstep; t++) {
    TString treeName(TString("T_") + TString(std::to_string( t ).c_str() ));
    TString treeTitle(TString("TemperatureTree_") + TString(std::to_string( t ).c_str() ));
    TTree * tree = new TTree(treeName, treeTitle);
    
    setT(temp_array[t]);

    for(uint i = 0; i < dim_vector; i++) {
      tree -> Branch(TString("Lattice_") + TString( TString(std::to_string( i ).c_str() ) ), "Block", &block[i]);
      lattice_vector[i].cooling(I0);
      E_tmp[i] = lattice_vector[i].energy();
      M_tmp[i] = lattice_vector[i].magnetization();
      S_tmp[i] = ( (double) E_tmp[i])/ lattice_vector[i].getNumSpin();
      T_tmp[i] = Lattice::getT();
   
      block[i].setBlock(i, T_tmp[i], E_tmp[i], M_tmp[i], S_tmp[i], 0, I0);
    }
    tree -> Fill();
    for(uint j = 0; j < iter; j++) {
      for(uint i = 0; i < dim_vector; i++) {
        data = lattice_vector[i].coolingPar();
        E_tmp[i] += (int) data[1];
        M_tmp[i] += data[2];
        S_tmp[i] += data[3];
        T_tmp[i] = data[0];
        if(data[2] < -1) cout << "ERROR" << endl << std::flush;
        block[i].setBlock(i, T_tmp[i], E_tmp[i], M_tmp[i], S_tmp[i], 0, j + 1);
      }
      tree -> Fill();
    }
    f.Write();
    tree->Delete();
    //gDirectory->ls();
    std::cout << "[ done "<< (int) ( ( (double) t / tempstep ) * 100 ) << "% ] T = " << temp_array[t] << std::endl << std::flush;
  }
  delete[] E_tmp;
  delete[] M_tmp;
  delete[] S_tmp;
  delete[] T_tmp;
  delete[] data;
  
  f.Write();
  f.Close();
}
