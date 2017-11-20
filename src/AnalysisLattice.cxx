#include "TString.h"
#include "TFile.h"
#include "TList.h"

#include "Lattice.h"
#include "AnalysisLattice.h"

ClassImp(AnalysisLattice)

/* Constructors */
AnalysisLattice::AnalysisLattice() : file_name(""), dim_c(1), dim_t(1) {
  energy             = new double*[1];
  magnetization      = new double*[1];
  site_energy        = new double*[1];
  site_magnetization = new double*[1];

  temperature        = new double[1];

  energy[1]             = new double;
  magnetization[1]      = new double;
  site_energy[1]        = new double;
  site_magnetization[1] = new double;
}

AnalysisLattice::AnalysisLattice(const TString& fname) : file_name(fname) {
  TFile f(fname);
  TList * li = f.GetListOfKeys();
  uint entries = li->GetEntries();
  for(uint i = 0; i < entries; i++){
    l.push_back(li->At(i)->GetName());
  }

  count(entries);

  creation();

  for(uint i = 0; i < dim_t; i++){
    temperature[i] = std::stod(l[i*dim_c].substr(l[i*dim_c].find_last_of("_") + 1));
    for(uint j = 0; j < dim_c; j++){
      energy[j][i] = (* (Lattice *)f.Get(l[i * dim_c +j].c_str())).energy2(false);
      magnetization[j][i] = (* (Lattice *)f.Get(l[i * dim_c +j].c_str())).magnetization();
      site_energy[j][i] = energy[j][i] / (* (Lattice *)f.Get(l[i * dim_c +j].c_str())).getNumSpin();
      site_magnetization[j][i] = magnetization[j][i] / (* (Lattice *)f.Get(l[i * dim_c +j].c_str())).getNumSpin();
    }
  }
  f.Close();
}


AnalysisLattice::~AnalysisLattice(){
  for(uint i = 0; i < dim_c; i++){
    delete []energy[i];
    delete []magnetization[i];
    delete []site_energy[i];
    delete []site_magnetization[i];
  }
  delete []temperature;
  delete []energy;
  delete []magnetization;
  delete []site_energy;
  delete []site_magnetization;
}

/* Private Methods */
void AnalysisLattice::count(const uint& li){
  uint temp = 0;
  dim_c = 0;
  for(uint i = 0; i < li; i++){
    temp = std::stoi(l[i].substr(14, l[i].find_last_of("_") - 14));
    temp > dim_c ? dim_c = temp : temp = dim_c;
  }
  dim_c++;
  dim_t = li / dim_c;
}

void AnalysisLattice::creation(){
  temperature = new double[dim_t];
  energy             = new double*[dim_c];
  magnetization      = new double*[dim_c];
  site_energy        = new double*[dim_c];
  site_magnetization = new double*[dim_c];
  for(uint i = 0; i < dim_c; i++){
    energy[i]             = new double[dim_t];
    magnetization[i]      = new double[dim_t];
    site_energy[i]        = new double[dim_t];
    site_magnetization[i] = new double[dim_t];
  }
}

/* Public Methods */
TString AnalysisLattice::getFileName() const { return file_name; }

vector<string> AnalysisLattice::getList() const { return l; }

void AnalysisLattice::analysis(const uint& x, const uint& y){}


void AnalysisLattice::print(){
  for(uint i = 0; i < dim_c; i++){
    for(uint j = 0; j < dim_t; j++){
      cout << energy[i][j] << '\t';
    }
    cout << '\n';
  }
  cout << '\n';

  for(uint i = 0; i < dim_c; i++){
    for(uint j = 0; j < dim_t; j++){
      cout << magnetization[i][j] << '\t';
    }
    cout << '\n';
  }
  cout << '\n';

  for(uint i = 0; i < dim_c; i++){
    for(uint j = 0; j < dim_t; j++){
      cout << site_energy[i][j] << '\t';
    }
    cout << '\n';
  }
  cout << '\n';

  for(uint i = 0; i < dim_c; i++){
    for(uint j = 0; j < dim_t; j++){
      cout << site_magnetization[i][j] << '\t';
    }
    cout << '\n';
  }
  cout << '\n';
}
