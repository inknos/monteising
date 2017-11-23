#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLine.h"

#include "Lattice.h"
#include "AnalysisLattice.h"

#include <vector>

ClassImp(AnalysisLattice)

/* Constructors */
AnalysisLattice::AnalysisLattice() : file_name(""), dim_c(1), dim_t(1) {
  energy             = new double*[1];
  magnetization      = new double*[1];
  site_energy        = new double*[1];
  e_mean             = new double*[DIM_ERR];
  m_mean             = new double*[DIM_ERR];
  s_mean             = new double*[DIM_ERR];
  e_err              = new double*[DIM_ERR];
  m_err              = new double*[DIM_ERR];
  s_err              = new double*[DIM_ERR];

  temperature        = new double[1];

  energy[1]             = new double;
  magnetization[1]      = new double;
  site_energy[1]        = new double;
  
  for(uint i = 0; i < DIM_ERR; i++){
    e_mean[i] = new double[1];
    m_mean[i] = new double[1];
    s_mean[i] = new double[1];

    e_err[i]  = new double[1];
    m_err[i]  = new double[1];
    s_err[i]  = new double[1];

    targetx_min = new double[DIM_ERR];
    targety_min = new double[DIM_ERR]; 
    targetx_max = new double[DIM_ERR];
    targety_max = new double[DIM_ERR];
 }
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
    }
  }
  f.Close();
}

AnalysisLattice::AnalysisLattice(const AnalysisLattice& obj) :
  file_name(obj.file_name), dim_c(obj.dim_c), dim_t(obj.dim_t) {
  creation();
  for(uint i = 0; i < dim_t; i++){
    temperature[i] = obj.temperature[i];
    for(uint j = 0; j < dim_c; i++){
      energy[j][i] = obj.energy[j][i];
      magnetization[j][i] = obj.magnetization[j][i];
      site_energy[j][i] = obj.site_energy[j][i];      
    }
    for(uint j = 0; j < DIM_ERR; j++){
      e_mean[j][i] = obj.e_mean[j][i];
      m_mean[j][i] = obj.e_mean[j][i];
      s_mean[j][i] = obj.e_mean[j][i];
      e_err[j][i]  = obj.e_err[j][i];
      m_err[j][i]  = obj.e_err[j][i];
      s_err[j][i]  = obj.e_err[j][i];
    }
  }
  for(uint i = 0; i < DIM_ERR; i++){
    targetx_min[i] = obj.targetx_min[i];
    targety_min[i] = obj.targety_min[i];
    targetx_max[i] = obj.targetx_max[i];
    targety_max[i] = obj.targetx_max[i];
  }
}

AnalysisLattice& AnalysisLattice::operator=(const AnalysisLattice& obj){
  if(&obj == this) return *this;
  this -> ~AnalysisLattice();
  new(this) AnalysisLattice(obj);
  return *this;
}

AnalysisLattice::~AnalysisLattice(){
  for(uint i = 0; i < dim_c; i++){
    delete []energy[i];
    delete []magnetization[i];
    delete []site_energy[i];
  }

  for(uint i = 0; i < DIM_ERR; i++){
    delete []e_mean[i];
    delete []m_mean[i];
    delete []s_mean[i];
    delete []e_err[i];
    delete []m_err[i];
    delete []s_err[i];
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

  delete []targetx_min;
  delete []targety_min;
  delete []targetx_max;
  delete []targety_max;
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
  temperature   = new double[dim_t];
  energy        = new double*[dim_c];
  magnetization = new double*[dim_c];
  site_energy   = new double*[dim_c];

  e_mean        = new double*[DIM_ERR];
  m_mean        = new double*[DIM_ERR];
  s_mean        = new double*[DIM_ERR];
  e_err         = new double*[DIM_ERR];
  m_err         = new double*[DIM_ERR];
  s_err         = new double*[DIM_ERR];

  targetx_min   = new double[DIM_ERR];
  targety_min   = new double[DIM_ERR];
  targetx_max   = new double[DIM_ERR];
  targety_max   = new double[DIM_ERR];
  
  for(uint i = 0; i < dim_c; i++){
    energy[i]        = new double[dim_t];
    magnetization[i] = new double[dim_t];
    site_energy[i]   = new double[dim_t];
  }

  for(uint i = 0; i < DIM_ERR; i++){
    e_mean[i] = new double[dim_t];
    m_mean[i] = new double[dim_t];
    s_mean[i] = new double[dim_t];
    e_err[i]  = new double[dim_t];
    m_err[i]  = new double[dim_t];
    s_err[i]  = new double[dim_t];
  }
}

/* Public Methods */
TString AnalysisLattice::getFileName() const { return file_name; }

vector<string> AnalysisLattice::getList() const { return l; }

TMultiGraph * AnalysisLattice::analysisNoErr(const uint& x, const uint& y,
                                             const uint& err,
                                             const TString& name, const TString& title,
                                             const bool& target) {
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

  if(target){
    TLine *** line = new TLine ** [DIM_ERR];
    for(uint i = 0; i < DIM_ERR; i++){
      line[i] = new TLine * [4];
    }
    for(uint i = 0; i < DIM_ERR; i++){
      line[i][0] = new TLine(targetx_min[i], targety_min[i], targetx_max[i], targety_min[i]);
      line[i][1] = new TLine(targetx_min[i], targety_max[i], targetx_max[i], targety_max[i]);
      line[i][2] = new TLine(targetx_max[i], targety_min[i], targetx_max[i], targety_max[i]);
      line[i][3] = new TLine(targetx_min[i], targety_min[i], targetx_min[i], targety_max[i]);
    }

    for(uint i = 0; i < DIM_ERR; i++){
      for(uint j = 0; j < 4; j++){
        line[i][j] -> SetLineColorAlpha(2 + i*2, 0.8);
        line[i][j] -> Draw("SAME");
      }
    }
  }

  return mg;
}

TMultiGraph * AnalysisLattice::analysisErr(const uint& x, const uint& y,
                                           const uint& err,
                                           const TString& name, const TString& title,
                                           const bool& target){
  void * x_vec;
  void * x_err;
  void * y_vec;
  void * y_err;
  double * terr = new double[dim_t];
  for(uint i = 0; i < dim_t; i++) terr[i] = 0.;

  TGraphErrors ** gr = new TGraphErrors*[DIM_ERR];
  
  TMultiGraph * mg = new TMultiGraph(name, title);

  uint init = 0, end = 0;
  switch(err){
  case 0: { init = 0; end = 1; break; }
  case 1: { init = 1; end = 2; break; }
  case 2: { init = 0; end = DIM_ERR; break; }
  default: {init = 0; end = DIM_ERR; break; }
  }

  for(uint i = init; i < end; i++){
    switch(x){
    case 1:{ x_vec = e_mean[i]; x_err = e_err[i]; break; }
    case 2:{ x_vec = temperature; x_err = terr; break; }
    case 3:{ x_vec = m_mean[i]; x_err = m_err[i]; break; }
    case 4:{ x_vec = s_mean[i]; x_err = s_err[i]; break; }
    default:{
      std::cout << "Falling to default case for x\n" << std::flush;
      x_vec = temperature; x_err = terr;
      break;
    }
    }
    switch(y){
    case 1:{ y_vec = e_mean[i]; y_err = e_err[i]; break; }
    case 2:{ y_vec = temperature; y_err = terr; break; }
    case 3:{ y_vec = m_mean[i]; y_err = m_err[i]; break; }
    case 4:{ y_vec = s_mean[i]; y_err = s_err[i]; break; }
    default:{
      std::cout << "Falling to default case for y\n" << std::flush;
      y_vec = e_mean[i]; y_err = e_err[i];
      break;
    }
    }
    gr[i] = new TGraphErrors(dim_t, (double*) x_vec, (double*) y_vec,
                             (double*) x_err, (double*) y_err);
    gr[i]->SetMarkerStyle(23 - i);
    gr[i]->SetMarkerColor(2 + i * 2);
  
    mg->Add(gr[i]);
  }

  mg->Draw("ALP");
  
  if(target){
    TLine *** line = new TLine ** [DIM_ERR];
    for(uint i = 0; i < DIM_ERR; i++){
      line[i] = new TLine * [4];
    }
    for(uint i = 0; i < DIM_ERR; i++){
      line[i][0] = new TLine(targetx_min[i], targety_min[i], targetx_max[i], targety_min[i]);
      line[i][1] = new TLine(targetx_min[i], targety_max[i], targetx_max[i], targety_max[i]);
      line[i][2] = new TLine(targetx_max[i], targety_min[i], targetx_max[i], targety_max[i]);
      line[i][3] = new TLine(targetx_min[i], targety_min[i], targetx_min[i], targety_max[i]);
    }

    for(uint i = 0; i < DIM_ERR; i++){
      for(uint j = 0; j < 4; j++){
        line[i][j] -> SetLineColorAlpha(2 + i*2, 0.8);
        line[i][j] -> Draw("SAME");
      }
    }
  }

  delete []terr;

  return mg;
}

TMultiGraph * AnalysisLattice::analysis(const uint& x, const uint& y, const uint& err, const TString& name, const TString& title, const bool& target = false){
  TMultiGraph * (AnalysisLattice::*func)(const uint& _x, const uint& _y, const uint& _err, const TString& _name, const TString& _title, const bool& _target);
  switch(err){
  case 0:{ func = &AnalysisLattice::analysisErr; break; }
  case 1:{ func = &AnalysisLattice::analysisErr; break; }
  case 2:{ func = &AnalysisLattice::analysisErr;   break; }
  case 3:{ func = &AnalysisLattice::analysisNoErr; break;}
  default:{ return 0; }
  }
  return (*this.*func)(x, y, err, name, title, target);
}

double* AnalysisLattice::getTargetX(const uint& spin){
  double * ret = new double[2];
  ret[0] = targetx_min[spin];
  ret[1] = targetx_max[spin];
  return ret;
}

double* AnalysisLattice::getTargetY(const uint& spin){
  double * ret = new double[2];
  ret[0] = targety_min[spin];
  ret[1] = targety_max[spin];
  return ret;
}

double* AnalysisLattice::getTarget(const uint& spin){
  double * ret = new double[4];
  ret[0] = targetx_min[spin];
  ret[1] = targetx_max[spin];
  ret[2] = targety_min[spin];
  ret[3] = targety_max[spin];
  return ret;
}

void AnalysisLattice::setTarget(const double& x_min, const double& x_max,
                                const double& y_min, const double& y_max,
                                const uint& dim, const uint& x_axis, const uint& y_axis){
  void * x_vec;
  void * y_vec;
  vector<uint> index;

  targetx_min[dim] = x_min;
  targety_min[dim] = y_min;
  targetx_max[dim] = x_max;
  targety_max[dim] = y_max;
  
  for(uint j = 0; j < dim_t; j++){
    e_mean[dim][j] = 0; m_mean[dim][j] = 0; s_mean[dim][j] = 0;
    e_err[dim][j]  = 0; m_err[dim][j]  = 0; s_err[dim][j]  = 0;
  }

  for(uint i = 0; i < dim_c; i++) {
    switch(x_axis){
    case 1:{ x_vec = energy[i]; break; }
    case 2:{ x_vec = temperature; break; }
    case 3:{ x_vec = magnetization[i]; break; }
    case 4:{ x_vec = site_energy[i]; break; }
    default:{
      std::cout << "Falling to default case for x\n" << std::flush;
      x_vec = temperature;
      break;
    }
    }
    switch(y_axis){
    case 1:{ y_vec = energy[i]; break; }
    case 2:{ y_vec = temperature; break; }
    case 3:{ y_vec = magnetization[i]; break; }
    case 4:{ y_vec = site_energy[i]; break; }
    default:{
      std::cout << "Falling to default case for y\n" << std::flush;
      y_vec = energy[i];
      break;
    }
    }
    for(uint j = 0; j < dim_t; j++){
      if( ((double*) x_vec)[j] >= x_min &&
          ((double*) x_vec)[j] <= x_max &&
          ((double*) y_vec)[j] >= y_min &&
          ((double*) y_vec)[j] <= y_max ) {
        // non contare
        index.push_back(i);
        break;
      }
    }
  }

  for(uint i = 0; i < dim_c; i++){
    if(std::find(index.begin(), index.end(), i) != index.end()) {
      for(uint j = 0; j < dim_t; j++){
        // somma i valori per la media
        e_mean[dim][j] += energy[i][j];
        m_mean[dim][j] += magnetization[i][j];
        s_mean[dim][j] += site_energy[i][j];
      }
    }
    else { continue; }
  }
  for(uint i = 0; i < dim_t; i++) {
    e_mean[dim][i] /= index.size();
    m_mean[dim][i] /= index.size();
    s_mean[dim][i] /= index.size();
  }

  for(uint i = 0; i < dim_c; i++){
    if(std::find(index.begin(), index.end(), i) != index.end()) {
      for(uint j = 0; j < dim_t; j++){ 
        e_err[dim][j] += pow( e_mean[dim][j] - energy[i][j] , 2);
        m_err[dim][j] += pow( m_mean[dim][j] - magnetization[i][j] , 2);
        s_err[dim][j] += pow( s_mean[dim][j] - site_energy[i][j] , 2);
      }
    }else{ continue; }
  }
  for(uint i = 0; i < dim_t; i++){
    e_err[dim][i] /= index.size();
    m_err[dim][i] /= index.size();
    s_err[dim][i] /= index.size();
    
    e_err[dim][i] = sqrt(e_err[dim][i]);
    m_err[dim][i] = sqrt(m_err[dim][i]);
    s_err[dim][i] = sqrt(s_err[dim][i]);
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

  for(uint i = 0; i < DIM_ERR; i++){
    for(uint j = 0; j < dim_t; j++){
      cout << e_mean[i][j] << '\t';
    }
    cout << '\n';
  }
  cout << '\n';

  for(uint i = 0; i < DIM_ERR; i++){
    for(uint j = 0; j < dim_t; j++){
      cout << m_mean[i][j] << '\t';
    }
    cout << '\n';
  }
  cout << '\n';

  for(uint i = 0; i < DIM_ERR; i++){
    for(uint j = 0; j < dim_t; j++){
      cout << s_mean[i][j] << '\t';
    }
    cout << '\n';
  }
  cout << '\n';

  for(uint i = 0; i < DIM_ERR; i++){
    for(uint j = 0; j < dim_t; j++){
      cout << e_err[i][j] << '\t';
    }
    cout << '\n';
  }
  cout << '\n';

  for(uint i = 0; i < DIM_ERR; i++){
    for(uint j = 0; j < dim_t; j++){
      cout << m_err[i][j] << '\t';
    }
    cout << '\n';
  }
  cout << '\n';

  for(uint i = 0; i < DIM_ERR; i++){
    for(uint j = 0; j < dim_t; j++){
      cout << s_err[i][j] << '\t';
    }
    cout << '\n';
  }
  cout << '\n';
}
