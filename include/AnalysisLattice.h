#ifndef ANALYSIS_LATTICE_H
#define ANALYSIS_LATTICE_H

#include "Lattice.h"
#include "TString.h"
#include "TFile.h"

#include <vector>
#include <string>

#define aENERGY        1
#define aTEMPERATURE   2
#define aMAGNETIZATION 3
#define aSITE_ENERGY   4

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
  double ** site_magnetization;
  Lattice getLattice(const uint& i);

  double * e_mean;
  double * m_mean;
  double * s_mean;

  double * e_err;
  double * m_err;
  double * s_err;
  
  void count(const uint&);

  void creation();

  void analysisNoErr(const uint&, const uint&, const TString&, const TString&);
  
  void analysisErr(const uint&, const uint&, const TString&, const TString&);
  
 public:
  AnalysisLattice();

  AnalysisLattice(const TString&);

  ~AnalysisLattice();

  TString getFileName() const;

  vector<string> getList() const;

  void analysis(const uint&, const uint&);

  void print();
  
  ClassDef(AnalysisLattice, 1)
    };

#endif
