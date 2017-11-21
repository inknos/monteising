#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"

#include "Lattice.h"
#include "AnalysisLattice.h"

ClassImp(AnalysisLattice)

/* Constructors */
AnalysisLattice::AnalysisLattice() : file_name(""), dim_c(1), dim_t(1) {
  energy             = new double*[1];
  magnetization      = new double*[1];
  site_energy        = new double*[1];

  temperature        = new double[1];

  energy[1]             = new double;
  magnetization[1]      = new double;
  site_energy[1]        = new double;

  e_mean = new double[1];
  m_mean = new double[1];
  s_mean = new double[1];

  e_err  = new double[1];
  m_err  = new double[1];
  s_err  = new double[1];
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
    e_mean[i] = 0;
    m_mean[i] = 0;
    s_mean[i] = 0;
    
    for(uint j = 0; j < dim_c; j++){
      energy[j][i] = (* (Lattice *)f.Get(l[i * dim_c +j].c_str())).energy2(false);
      magnetization[j][i] = (* (Lattice *)f.Get(l[i * dim_c +j].c_str())).magnetization();
      site_energy[j][i] = energy[j][i] / (* (Lattice *)f.Get(l[i * dim_c +j].c_str())).getNumSpin();

      e_mean[i] += energy[j][i];
      m_mean[i] += magnetization[j][i];
      s_mean[i] += site_energy[j][i];
    }
    e_mean[i] /= dim_c;
    m_mean[i] /= dim_c;
    s_mean[i] /= dim_c;
  }
  
  for(uint i = 0; i < dim_t; i++){
    e_err[i] = 0;
    m_err[i] = 0;
    s_err[i] = 0;
    
    for(uint j = 0; j < dim_c; j++){
      e_err[i] += pow( e_mean[i] - energy[j][i] , 2);
      m_err[i] += pow( m_mean[i] - magnetization[j][i] , 2);
      s_err[i] += pow( s_mean[i] - site_energy[j][i] , 2);
    }
    e_err[i] /= ( dim_c <= 1 ? 1 : dim_c - 1);
    m_err[i] /= ( dim_c <= 1 ? 1 : dim_c - 1);
    s_err[i] /= ( dim_c <= 1 ? 1 : dim_c - 1);

    e_err[i] = sqrt(e_err[i]);
    m_err[i] = sqrt(m_err[i]);
    s_err[i] = sqrt(s_err[i]);
  }
  
}

f.Close();
}


AnalysisLattice::~AnalysisLattice(){
  for(uint i = 0; i < dim_c; i++){
    delete []energy[i];
    delete []magnetization[i];
    delete []site_energy[i];
  }
  delete []temperature;
  delete []energy;
  delete []magnetization;
  delete []site_energy;

  delete []e_mean;
  delete []m_mean;
  delete []s_mean;

  delete []e_err;
  delete []m_err;
  delete []s_err;
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

  e_mean = new double[dim_t];
  m_mean = new double[dim_t];
  s_mean = new double[dim_t];

  e_err  = new double[dim_t];
  m_err  = new double[dim_t];
  s_err  = new double[dim_t];

  for(uint i = 0; i < dim_c; i++){
    energy[i]             = new double[dim_t];
    magnetization[i]      = new double[dim_t];
    site_energy[i]        = new double[dim_t];
  }
}

/* Public Methods */
TString AnalysisLattice::getFileName() const { return file_name; }

vector<string> AnalysisLattice::getList() const { return l; }

void AnalysisLattice::analysisNoErr(const uint& x = 1, const uint& y = 0, const TString& name = "graph", const TString& title = "Graph"){
  TGraph ** gr = new TGraph * [dim_c];
  TMultiGraph * mg = new TMultiGraph(name, title);

  void * x_vec;
  void * y_vec;

  for(uint i = 0; i < dim_c; i++){
    switch(x){
    case 1:{ x_vec = energy[i]; break; }
    case 2:{ x_vec = temperature; break; }
    case 3:{ x_vec = magnetization[i]; break; }
    case 4:{ x_vec = site_energy[i]; break; }
    default:{
      std::cout << "Falling to default case for x\n" << std::flush;
      x_vec = &temperature;
      break;
    }
    }
    switch(y){
    case 1:{ y_vec = energy[i]; break; }
    case 2:{ y_vec = temperature; break; }
    case 3:{ y_vec = magnetization[i]; break; }
    case 4:{ y_vec = site_energy[i]; break; }
    default:{
      std::cout << "Falling to default case for y\n" << std::flush;
      y_vec = &energy[i];
      break;
    }
    }
    gr[i] = new TGraph(dim_t, (double*) x_vec, (double*) y_vec);
    gr[i] -> SetMarkerColor(kAzure - 10 + i);
    gr[i] -> SetMarkerStyle(20);
  }
  for(uint i = 0; i < dim_c; i++){
    mg -> Add(gr[i]);
  }

  mg -> Draw("AP");
}

void AnalysisLattice::analysisErr(const uint& x, const uint& y, const TString& name, const TString& title){
  return;
}

void AnalysisLattice::analysis(const uint& x, const uint& y, const uint& e, const TString& name, const TString& title){
  void (AnalysisLattice::*funct)(const uint& x, const uint& y, const TString& name, const TString& title);
  switch(e){
  case 1:{ func = &AnalysisLattice::analysisNoErr; break; }
  case 2:{ func = &AnalysisLattice::analysisErr;   break; }
  case 3:{ return; }
  default:{return; }
  }
}

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
}
