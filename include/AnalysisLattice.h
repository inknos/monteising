#ifndef ANALYSIS_LATTICE_H
#define ANALYSIS_LATTICE_H

#include "Lattice.h"
#include "TString.h"
#include "TFile.h"

#include <vector>
#include <string>

class AnalysisLattice{
 private:
  TString file_name;
  std::vector<std::string> l;

 public:
  AnalysisLattice();

  AnalysisLattice(const TFile&);  

  AnalysisLattice(const TString&);

  TString getFileName() const;

  vector<string> getList() const;
  
  ClassDef(AnalysisLattice, 1)
    };

#endif
