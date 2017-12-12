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

void outputGrafici() {

  gStyle->SetOptFit();

  //AnalysisLattice a20("prova.root", "anal_20.root");
  AnalysisLattice a40("prova.root", "anal_40.root");
  AnalysisLattice a60("prova.root", "anal_60.root");
  AnalysisLattice a80("sim_80.root", "anal_80.root");
  AnalysisLattice a100("sim_100.root", "anal_100.root");

  TCanvas * c1 = new TCanvas("c1", "Canvas");
  c1 ->SetGrid();
  auto ET = new TMultiGraph();
  //TGraphErrors * grET20 = a20.drawLattice(0, TEMPERATURE, ENERGY);
  TGraphErrors * grET40 = a40.drawLattice(0, TEMPERATURE, ENERGY);
  TGraphErrors * grET60 = a60.drawLattice(0, TEMPERATURE, ENERGY);
  TGraphErrors * grET80 = a80.drawLattice(0, TEMPERATURE, ENERGY);
  TGraphErrors * grET100 = a100.drawLattice(0, TEMPERATURE, ENERGY);
  //grET20 -> SetMarkerColor(MYGREEN);
  grET40 -> SetMarkerColor(MYYELLOW);
  grET60 -> SetMarkerColor(MYRED);
  grET80 -> SetMarkerColor(MYPURPLE);
  grET100 -> SetMarkerColor(MYBLUE);

  ET -> SetTitle("Energy on Temperature for Different Size");
  //ET -> Add(grET20);
  ET -> Add(grET40);
  ET -> Add(grET60);
  ET -> Add(grET80);
  ET -> Add(grET100);
  ET -> Draw("ALP");
  //gPad->Modified();
  ET -> GetXaxis() -> SetTitle("Temeprature (#frac{J}{k_{b}})");
  ET -> GetYaxis() -> SetTitle("Energy (J)");

  auto l1 = new TLegend(0.1,0.7,0.48,0.9);
  l1 -> SetHeader("Legend","C"); // option "C" allows to center the header
  //l1 -> AddEntry(grET20,"Lattice 20x20","lep");
  l1 -> AddEntry(grET40,"Lattice 40x40","lep");
  l1 -> AddEntry(grET60,"Lattice 60x60","lep");
  l1 -> AddEntry(grET80,"Lattice 80x80","lep");
  l1 -> AddEntry(grET100,"Lattice 100x100","lep");
  l1 -> Draw();
  
  TCanvas * c2 = new TCanvas("c2", "Canvas");
  c2 -> SetGrid();
  auto MT = new TMultiGraph();
  //TGraphErrors * grMT20 = a20.drawLatticeMean(TEMPERATURE, MAG);
  TGraphErrors * grMT40 = a40.drawLatticeMean(TEMPERATURE, MAG);
  TGraphErrors * grMT60 = a60.drawLatticeMean(TEMPERATURE, MAG);
  TGraphErrors * grMT80 = a80.drawLatticeMean(TEMPERATURE, MAG);
  TGraphErrors * grMT100 = a100.drawLatticeMean(TEMPERATURE, MAG);

  auto fMT100 = new TF1("fMT100" , AnalysisLattice::analiticM , 1.98, 2.2, 3);
  fMT100 -> SetParameters(2.2, 0.12, 1.1);
  fMT100 -> SetParNames("Tc", "exp", "coeff");
  //fMT100 -> SetLineColor(FITBLUE);
  grMT100 -> Fit(fMT100, "R");

  //grMT20 -> SetMarkerColor(MYGREEN);
  grMT40 -> SetMarkerColor(MYYELLOW);
  grMT60 -> SetMarkerColor(MYRED);
  grMT80 -> SetMarkerColor(MYPURPLE);
  grMT100 -> SetMarkerColor(MYBLUE);

  MT -> SetTitle("Magnetization per Spin on Temperature for Different Size");
  //MT -> Add(grMT20);
  MT -> Add(grMT40);
  MT -> Add(grMT60);
  MT -> Add(grMT80);
  MT -> Add(grMT100);
  MT -> Draw("ALP");

  c2 -> Update();

  MT -> GetXaxis() -> SetTitle("Temeperature (#frac{J}{k_{b}})");
  MT -> GetYaxis() -> SetTitle("Magnetization per Spin (#mu)");

  c2 -> Update();
  
  auto sMT100 = (TPaveStats*) grMT100 -> GetListOfFunctions() -> FindObject("stats");
  //sMT80  -> SetTextColor(MYPURPLE);
  sMT100 -> SetTextColor(MYBLUE);
  c2 -> Modified();

  auto l2 = new TLegend(0.1,0.7,0.48,0.9);
  l2 -> SetHeader("Legend","C"); // option "C" allows to center the header
  //l2 -> AddEntry(grMT20,"Lattice 20x20","lep");
  l2 -> AddEntry(grMT40,"Lattice 40x40","lep");
  l2 -> AddEntry(grMT60,"Lattice 60x60","lep");
  l2 -> AddEntry(grMT80,"Lattice 80x80","lep");
  //l2 -> AddEntry(fMT80, "Fit Function 80x80", "l");
  l2 -> AddEntry(grMT100,"Lattice 100x100","lep");
  l2 -> AddEntry(fMT100, "Fit Function 100x100", "l");
  l2 -> Draw();

  TCanvas * c3 = new TCanvas("c3", "Canvas");
  c3 -> SetGrid();
  auto XT = new TMultiGraph();
  //TGraphErrors * grXT20 = a20.drawLattice(0, TEMPERATURE, SUSC);
  TGraphErrors * grXT40 = a40.drawLattice(0, TEMPERATURE, SUSC);
  TGraphErrors * grXT60 = a60.drawLattice(0, TEMPERATURE, SUSC);
  TGraphErrors * grXT80 = a80.drawLattice(0, TEMPERATURE, SUSC);
  TGraphErrors * grXT100 = a100.drawLattice(0, TEMPERATURE, SUSC);
  //grXT20 -> SetMarkerColor(MYGREEN);
  grXT40 -> SetMarkerColor(MYYELLOW);
  grXT60 -> SetMarkerColor(MYRED);
  grXT80 -> SetMarkerColor(MYPURPLE);
  grXT100 -> SetMarkerColor(MYBLUE);

  auto fXT100 = new TF1("fXT100" , AnalysisLattice::analiticX , 1.9, 2.2, 3);
  fXT100 -> SetParameters(2.2, 1.175, 1);
  fXT100 -> SetParNames("Tc", "exp", "coeff");
  //fMT100 -> SetLineColor(FITBLUE);
  grXT100 -> Fit(fXT100, "R");

  XT -> SetTitle("Magnetic Susceptibility on Temperature for Different Size");
  //XT -> Add(grXT20);
  //XT -> Add(grXT40);
  //XT -> Add(grXT60);
  //XT -> Add(grXT80);
  XT -> Add(grXT100);
  XT -> Draw("ALP");
  //gPad->Modified();
  XT -> GetXaxis() -> SetTitle("Temeprature (#frac{J}{k_{b}})");
  XT -> GetYaxis() -> SetTitle("Magnetic Susceptibility (#frac{#mu}{k_{B}})");

  c3 -> Update();
  
  auto sXT100 = (TPaveStats*) grXT100 -> GetListOfFunctions() -> FindObject("stats");
  //sMT80  -> SetTextColor(MYPURPLE);
  sXT100 -> SetTextColor(MYBLUE);
  c3 -> Modified();
  
  auto l3 = new TLegend(0.1,0.7,0.48,0.9);
  l3 -> SetHeader("Legend","C"); // option "C" allows to center the header
  //l3 -> AddEntry(grXT20,"Lattice 20x20","lep");
  //l3 -> AddEntry(grXT40,"Lattice 40x40","lep");
  //l3 -> AddEntry(grXT60,"Lattice 60x60","lep");
  //l3 -> AddEntry(grXT80,"Lattice 80x80","lep");
  l3 -> AddEntry(grXT100,"Lattice 100x100","lep");
  l3 -> AddEntry(fXT100, "Fit Function 100x100", "l");
  l3 -> Draw();

  //TCanvas * c4 = new TCanvas("c4", "Canvas");
  //a100.evalBinning(30) -> Draw("ALP");
  
}
