#include "AnalysisLattice.h"

#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"

#define MYRED     kRed - 4
#define MYGREEN   kGreen - 2
#define MYYELLOW  kYellow - 6
#define MYBLUE    kBlue - 4
#define MYPURPLE  kMagenta -2

#define FITPURPLE kMagenta -4
#define FITBLUE   kCyan + 1

void outputGrafici3D() {

  gStyle->SetOptFit();
  AnalysisLattice a40("prova.root", "provaOut3d40.root");

  TCanvas * c1 = new TCanvas("c1", "Canvas");
  c1 -> SetGrid();
  auto ET = new TMultiGraph();
  TGraphErrors * grET40 = a40.drawLattice(0, TEMPERATURE, ENERGY);
  grET40 -> SetMarkerColor(MYBLUE);

  ET -> SetTitle("Energy on Temperature for Different Size");
  ET -> Add(grET40);
  ET -> Draw("ALP");
  //gPad->Modified();
  ET -> GetXaxis() -> SetTitle("Temeprature (#frac{J}{k_{b}})");
  ET -> GetYaxis() -> SetTitle("Energy (J)");

  auto l1 = new TLegend(0.1,0.7,0.48,0.9);
  l1 -> SetHeader("Legend","C"); // option "C" allows to center the header
  l1 -> AddEntry(grET40,"Lattice 40x40","lep");
  l1 -> Draw();
  
  TCanvas * c2 = new TCanvas("c2", "Canvas");
  c2 -> SetGrid();
  auto MT = new TMultiGraph();
  TGraphErrors * grMT40 = a40.drawLatticeMean(TEMPERATURE, MAG);

  auto fMT40 = new TF1("fMT40" , AnalysisLattice::analiticM , 3.75, 4.1, 3);
  fMT40 -> SetParameters(4.5, 0.32, 1.6);
  fMT40 -> SetParNames("Tc", "exp", "coeff");
  //fMT100 -> SetLineColor(FITBLUE);
  grMT40 -> Fit(fMT40, "R");

  grMT40 -> SetMarkerColor(MYBLUE);

  MT -> SetTitle("Magnetization per Spin on Temperature for Different Size");
  MT -> Add(grMT40);
  MT -> Draw("ALP");

  c2 -> Update();

  MT -> GetXaxis() -> SetTitle("Temeperature (#frac{J}{k_{b}})");
  MT -> GetYaxis() -> SetTitle("Magnetization per Spin (#mu)");

  c2 -> Update();
  
  auto sMT40 = (TPaveStats*) grMT40 -> GetListOfFunctions() -> FindObject("stats");
  sMT40 -> SetTextColor(MYBLUE);
  c2 -> Modified();

  auto l2 = new TLegend(0.1,0.7,0.48,0.9);
  l2 -> SetHeader("Legend","C"); // option "C" allows to center the header
  l2 -> AddEntry(grMT40,"Lattice 40x40","lep");
  l2 -> AddEntry(fMT40, "Fit Function 40x40", "l");
  l2 -> Draw();

  TCanvas * c3 = new TCanvas("c3", "Canvas");
  c3 -> SetGrid();
  auto XT = new TMultiGraph();
  TGraphErrors * grXT40 = a40.drawLattice(0, TEMPERATURE, SUSC);
  grXT40 -> SetMarkerColor(MYBLUE);

  auto fXT40 = new TF1("fXT40" , AnalysisLattice::analiticX , 4, 4.4, 3);
  fXT40 -> SetParameters(4.5, 1.2, 1);
  fXT40 -> SetParNames("Tc", "exp", "coeff");
  //fMT100 -> SetLineColor(FITBLUE);
  grXT40 -> Fit(fXT40, "R");

  XT -> SetMaximum(180.);
  XT -> SetTitle("Magnetic Susceptibility on Temperature for Different Size");
  XT -> Add(grXT40);
  XT -> Draw("ALP");
  XT -> GetXaxis() -> SetTitle("Temeprature (#frac{J}{k_{b}})");
  XT -> GetYaxis() -> SetTitle("Magnetic Susceptibility (#frac{#mu}{k_{B}})");

  c3 -> Update();
  
  auto sXT40 = (TPaveStats*) grXT40 -> GetListOfFunctions() -> FindObject("stats");
  sXT40 -> SetTextColor(MYBLUE);
  c3 -> Modified();
  
  auto l3 = new TLegend(0.1,0.7,0.48,0.9);
  l3 -> SetHeader("Legend","C"); // option "C" allows to center the header
  l3 -> AddEntry(grXT40,"Lattice 40x40","lep");
  l3 -> AddEntry(fXT40, "Fit Function 40x40", "l");
  l3 -> Draw();

  //TCanvas * c4 = new TCanvas("c4", "Canvas");
  //a100.evalBinning(30) -> Draw("ALP");
  
}
