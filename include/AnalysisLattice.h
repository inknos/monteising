#ifndef ANALYSIS_LATTICE_H
#define ANALYSIS_LATTICE_H

#include "Lattice.h"
#include "TString.h"
#include "TFile.h"
#include "TMultiGraph.h"

#include <vector>
#include <string>

#define ENERGY        1
#define TEMPERATURE   2
#define MAGNETIZATION 3
#define SITE_ENERGY   4

#define NO_ERR        1
#define ERR           2
#define BOTH          3

#define DIM_ERR       2

#define SPIN_UP       1
#define SPIN_DOWN     0

typedef unsigned int uint;

class AnalysisLattice{
 private:
  TString file_name;
  std::vector<std::string> l;

  uint dim_c;
  uint dim_t;
  
  double *  temperature;
  double ** energy;
  double ** magnetization;
  double ** site_energy;

  double ** e_mean;
  double ** m_mean;
  double ** s_mean;

  double ** e_err;
  double ** m_err;
  double ** s_err;
  
  void count(const uint&);

  void creation();

  TMultiGraph * analysisNoErr(const uint&, const uint&, const TString&, const TString&);

  TMultiGraph *analysisErr(const uint&, const uint&, const TString&, const TString&);
  
 public:
  AnalysisLattice();

  AnalysisLattice(const TString&);

  AnalysisLattice(const AnalysisLattice&);
  
  ~AnalysisLattice();

  AnalysisLattice& operator=(const AnalysisLattice&);

  Lattice getLattice(const uint& i);
 
  TString getFileName() const;

  vector<string> getList() const;

  void setTarget(const double&, const double&, const double&, const double&, const uint&, const uint&, const uint&);
  
  void analysis(const uint&, const uint&, const uint&, const TString&, const TString&);

  void print();

  ClassDef(AnalysisLattice, 1)
    };

#endif
