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

#define DIM_ERR       2 // dimension of the analysis vectors = num spins

// for analysis()
#define BOTH         -1
#define SPIN_DOWN     0 // graph only spin down with errors
#define SPIN_UP       1 // graph only spin up   with errors
#define ERR           2 // graph both spins     with errors
#define NO_ERR        3 // graph raw data

#define TARGET        1
#define NO_TARGET     0

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

  double * targetx_min;
  double * targety_min;
  double * targetx_max;
  double * targety_max;
  
  void count(const uint&);

  void creation();

  TMultiGraph * analysisNoErr(const uint&, const uint&, const uint&, const TString&, const TString&, const bool&);

  TMultiGraph *analysisErr(const uint&, const uint&, const uint&, const TString&, const TString&, const bool&);
  
 public:
  AnalysisLattice();

  AnalysisLattice(const TString&);

  AnalysisLattice(const AnalysisLattice&);
  
  ~AnalysisLattice();

  AnalysisLattice& operator=(const AnalysisLattice&);

  Lattice getLattice(const uint& i);
 
  TString getFileName() const;

  vector<string> getList() const;

  double * getTargetX(const uint&);

  double * getTargetY(const uint&);

  double * getTarget(const uint&);
  
  void setTarget(const double&, const double&, const double&, const double&, const uint&, const uint&, const uint&);
  
  TMultiGraph * analysis(const uint&, const uint&, const uint&, const TString&, const TString&, const bool&);

  void print();

  ClassDef(AnalysisLattice, 1)
    };

#endif
