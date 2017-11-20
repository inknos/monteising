#ifndef ANALYSIS_LATTICE_H
#define ANALYSIS_LATTICE_H

#include "Lattice.h"
#include "TString.h"
#include "TFile.h"

#include <vector>
#include <string>

#define ENERGY        1
#define TEMPERATURE   2
#define MAGNETIZATION 3
#define SITE_ENERGY   4

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

  void count(const uint&);

  void creation();

 public:
  AnalysisLattice();

  AnalysisLattice(const TString&);

  AnalysisLattice(const AnalysisLattice&);
  
  ~AnalysisLattice();

  AnalysisLattice& operator=(const AnalysisLattice&);

  Lattice getLattice(const uint& i);
 
  TString getFileName() const;

  vector<string> getList() const;

  void analysis(const uint&, const uint&, const TString&, const TString&);

  void print();

  ClassDef(AnalysisLattice, 1)
    };

#endif
