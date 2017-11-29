#ifndef ANALYSIS_LATTICE_H
#define ANALYSIS_LATTICE_H

#include "TString.h"

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

class AnalysisLattice {
 private:
  TString file_in;
  TString file_out;

 public:
  AnalysisLattice(const TString& file_input, const TString& file_output);


  TString getFileIn() const;

  TString getFileOut() const;

  void setFileIn(const TString& file_input);

  void setFileOut(const TString& file_output);

  void read();

  void write();
  
  void run();
 
};

#endif
