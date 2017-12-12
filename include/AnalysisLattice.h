#ifndef ANALYSIS_LATTICE_H
#define ANALYSIS_LATTICE_H

#include <vector>

#include "TString.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TGraph.h"

#define ENERGY         1
#define TEMPERATURE    2
#define TEMP           2
#define MAGNETIZATION  3
#define MAG            3
#define SITE_ENERGY    4
#define SENERGY        4
#define SUSCEPTIBILITY 5
#define SUSC           5


#define DIM_ERR       2 // dimension of the analysis vectors = num spins

// for analysis()
#define BOTH         -1
#define SPIN_DOWN     0 // graph only spin down with errors
#define SPIN_UP       1 // graph only spin up   with errors
#define ERR           2 // graph both spins     with errors
#define NO_ERR        3 // graph raw data

#define ABS           1
#define NO_ABS        0

#define TARGET        1
#define NO_TARGET     0

typedef unsigned int       uint;
typedef const unsigned int cuint;
typedef const TString      ctstring;
typedef const bool         cbool;

class AnalysisLattice {
 private:
  TString file_in;
  TString file_out;
  static double TempCritic;

  std::vector<double> bin(const std::vector<double>& vec);
  std::vector<double> binN(cuint& num_bin, const std::vector<double>& vec);

 public:

  AnalysisLattice(const TString& file_input, const TString& file_output);

  void run();

  TGraphErrors * drawLattice(cuint& lattice_number,
                             cuint& x_axis,
                             cuint& y_axis);

  TGraphErrors * drawLatticeMean(cuint& x_axis,
                                 cuint& y_axis);

  TMultiGraph * draw(cuint& x_axis, cuint& y_axis);


  static double analiticX(double * x, double * par);
  static double analiticM(double * x , double * par);

  void findTcritic();

  void fitLattice         ( bool mean,
                            cuint& x_axis,
                            cuint& y_axis,
                            double fit_temp_min,
                            double fit_temp_max,
                            int lat_number
                          );


  static double getTempCritic();

  static void setTempCritic(const double& _TempCritic);

  TString getFileIn() const;

  TString getFileOut() const;

  void setFileIn(const TString& file_input);

  void setFileOut(const TString& file_output);

  TGraph * evalBinning(cuint& nb);

};

#endif
