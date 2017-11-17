#include "Lattice.h"

#include "TGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"

#include <vector>
#include <math.h>

using namespace std;

void draw(int jjj = 0){
  Lattice a(100,2);

  if(jjj != 0){
    for(int i = 0; i < jjj; i++){
      a.cooling();
    }
  }

  int n  = a.getN();
  int ns = a.getNumSpin();
  int index = 0;

  int zeros = 0;
  int ones  = 0;

  for(int i = 0; i < ns; i++){
    a.getSpin(i) ? ones++ : zeros++;
  }
  cout << zeros << " " << ones << " " << zeros + ones << endl;

  TCanvas* c1 = new TCanvas("c1","Ising Canvas",200,10,700,500);

  TMultiGraph *mg = new TMultiGraph();
   mg->SetTitle("Ising");
  
  TGraph * g0 = new TGraph(zeros);
  for(int i = 0; i < ns; i++){
    if(!a.getSpin(i)){
      g0->SetPoint(index,
                   (int) 1 + (i % n),
                   (int) (n - i / n)
                   );
      index++;
    }
  }
  g0->SetMarkerStyle(21);
  g0->SetMarkerColor(kBlue);
  //g0->GetXaxis()->SetRangeUser(0., 100.);
  //g0->GetYaxis()->SetRangeUser(0., 100.);

  index = 0;
  TGraph * g1 = new TGraph(ones);
  for(int i = 0; i < ns; i++){
    if(a.getSpin(i)){
      g1->SetPoint(index,
                   (int) 1 + (i % n),
                   (int) (n - i / n)
                   );
      index++;
    }
  }
  g1->SetMarkerStyle(21);
  g1->SetMarkerColor(kRed);
  //g1->GetXaxis()->SetRangeUser(0., 100.);
  //g1->GetYaxis()->SetRangeUser(0., 100.);

  mg->Add(g0);
  mg->Add(g1);
  mg->Draw("AP");
}

void draw2(int jjj = 0){
  Lattice a(30,3);

  if(jjj != 0){
    for(int i = 0; i < jjj; i++){
      a.cooling();
    }
  }
  int n  = a.getN();
  int ns = a.getNumSpin();
  int index = 0;

  //cout << count << endl;
  TCanvas* c1 = new TCanvas("c1","Ising Canvas",200,10,700,500);

  TGraph2D * g0 = new TGraph2D();
  for(int i = 0; i < ns; i++){
    if(!a.getSpin(i)){
      g0->SetPoint(index,
                   (int) 1 + (i % n),
                   (int) n - (i % (int) pow(n,2)) / n,
                   (int) n - (i / pow(n,2) )
                   );
      index++;
    }
  }
  g0->SetMarkerStyle(20);
  g0->SetMarkerSize(0.6);
  g0->SetMarkerColor(kBlue);
  g0->Draw("AP");


  index = 0;
  TGraph2D * g1 = new TGraph2D();
  for(int i = 0; i < ns; i++){
    if(a.getSpin(i)){
      g1->SetPoint(index,
                   (int) 1 + (i % n),
                   (int) n - (i % (int) pow(n,2)) / n,
                   (int) n - (i / pow(n,2) )
                   );
      index++;
    }
  }
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.6);
  g1->SetMarkerColor(kRed);
  g1->Draw("AP same");
}
