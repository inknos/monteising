#include "AnalysisLattice.h"
#include "Block.h"
#include "INFO.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TH1D.h"

#include <string>
#include <vector>

#include <iostream>
using namespace std;

AnalysisLattice::AnalysisLattice(const TString& file_input, const TString& file_output) :
  file_in(file_input), file_out(file_output) {}

TString AnalysisLattice::getFileIn() const { return file_in; }

TString AnalysisLattice::getFileOut() const { return file_out; }

void AnalysisLattice::setFileIn(const TString& file_input) { file_in = file_input; }

void AnalysisLattice::setFileOut(const TString& file_output) { file_out = file_output; }

void AnalysisLattice::run(){
  TH1D::SetDefaultSumw2(true);
  static INFO info;

  TFile ifile(file_in, "read");
  TTree * info_tree = (TTree*) ifile.Get("info");
  TBranch * info_branch = info_tree -> GetBranch("Info");
  info_branch -> SetAddress(&info._tempmin);
  info_tree -> GetEvent(0);
  double tempmin = info._tempmin;
  double tempmax = info._tempmax;
  uint N         = info._N;
  uint dim       = info._dim;
  uint dim_vector= info._dim_vector;
  uint I0        = info._I0;
  uint iter      = info._iter;
  uint tempstep  = info._tempstep;
  uint datime    = info._datime_t;

  double** Tv = new double*[dim_vector];
  double** Ev = new double*[dim_vector];
  double** Mv = new double*[dim_vector];
  double** Sv = new double*[dim_vector];

  for(uint i = 0; i < dim_vector; i++){
    Tv[i] = new double[iter + 1];
    Ev[i] = new double[iter + 1];
    Mv[i] = new double[iter + 1];
    Sv[i] = new double[iter + 1];
  }

  double** T = new double*[dim_vector];
  double** E = new double*[dim_vector];
  double** M = new double*[dim_vector];
  double** S = new double*[dim_vector];
  double** X = new double*[dim_vector];
  double** dT= new double*[dim_vector];
  double** dE= new double*[dim_vector];
  double** dM= new double*[dim_vector];
  double** dS= new double*[dim_vector];
  double** dX= new double*[dim_vector];
  for(uint i = 0; i < dim_vector; i++){
    T[i] = new double[tempstep];
    E[i] = new double[tempstep];
    M[i] = new double[tempstep];
    S[i] = new double[tempstep];
    X[i] = new double[tempstep];
    dT[i]= new double[tempstep];
    dE[i]= new double[tempstep];
    dM[i]= new double[tempstep];
    dS[i]= new double[tempstep];
    dX[i]= new double[tempstep];
  }

  for(uint i = 0; i < dim_vector; i++){
    for(uint t = 0; t < tempstep; t++){
      T[i][t] = 0; E[i][t] = 0; M[i][t] = 0; S[i][t] = 0;
      dT[i][t] = 0; dE[i][t] = 0; dM[i][t] = 0; dS[i][t] = 0;
    }
  }
  
  double * TCriticVector =  new double[dim_vector];
  
  
  Block * block = 0;
  TTree * tree;
  for(uint t = 0; t < tempstep; t++) {
    //cout << "t " << t << endl << flush;
    TString treeName(TString("T_") + TString(std::to_string(t).c_str() ) );
    tree = (TTree*) ifile.Get(treeName);

    TH1D ** histoT = new TH1D*[dim_vector];
    TH1D ** histoE = new TH1D*[dim_vector];
    TH1D ** histoM = new TH1D*[dim_vector];
    TH1D ** histoS = new TH1D*[dim_vector];
    for(uint i = 0; i < dim_vector; i++){
      histoT[i] = new TH1D;
      histoE[i] = new TH1D;
      histoM[i] = new TH1D;
      histoS[i] = new TH1D;
    }

    for(uint i = 0; i < dim_vector; i++) {
      //cout << "i " << i << endl << flush;
      TString branchName(TString("Lattice_") + std::to_string(i).c_str() );
      tree -> SetBranchAddress(branchName, &block);
      // X = <M2> - <M>2 / T
      for(uint j = 0; j < iter + 1; j++) {
        tree -> GetEntry(j);
        Tv[i][j] = block -> T;
        Ev[i][j] = block -> E;
        Mv[i][j] = block -> M;
        Sv[i][j] = block -> S;

        T[i][t] += Tv[i][j];
        E[i][t] += Ev[i][j];
        M[i][t] += Mv[i][j];
        S[i][t] += Sv[i][j];
      }
      T[i][t] /= (iter + 1);
      E[i][t] /= (iter + 1);
      M[i][t] /= (iter + 1);
      S[i][t] /= (iter + 1);
      for(uint j = 0; j < iter + 1; j++){
        dT[i][t] += ( T[i][t] - Tv[i][j] ) * ( T[i][t] - Tv[i][j] );
        dE[i][t] += ( E[i][t] - Ev[i][j] ) * ( E[i][t] - Ev[i][j] );
        dM[i][t] += ( M[i][t] - Mv[i][j] ) * ( M[i][t] - Mv[i][j] );
        dS[i][t] += ( S[i][t] - Sv[i][j] ) * ( S[i][t] - Sv[i][j] );
      }
      dT[i][t] = TMath::Sqrt( dT[i][t] / ( iter + 1 ) );
      dE[i][t] = TMath::Sqrt( dE[i][t] / ( iter + 1 ) );
      dS[i][t] = TMath::Sqrt( dM[i][t] / ( iter + 1 ) );
      X[i][t] = dM[i][t] / T[i][t];
      dM[i][t] = TMath::Sqrt( dS[i][t] / ( iter + 1 ) );
    }
  }
  
  double TCriticMean = 0.;
  for(uint i = 0; i < dim_vector; i++){
    double X_max_tmp = 0.;
    uint t_max_tmp = 0;
    
    for(uint t = 0; t < tempstep; t++){
      if(X[i][t] > X_max_tmp){
        X_max_tmp = X[i][t];
        t_max_tmp = t;
      }
    }
    TCriticVector[i] = T[i][t_max_tmp];
    cout << "T critic of Lattice " << i << " is :" << TCriticVector[i] << endl << flush;
    
    TCriticMean += TCriticVector[i];
  }
  cout << "T critic mean is : " << TCriticMean / dim_vector << endl << flush;

  for(uint i = 0; i < dim_vector; i++){
    delete[] Tv[i];
    delete[] Ev[i];
    delete[] Mv[i];
    delete[] Sv[i];
  }
  delete[] Tv;
  delete[] Ev;
  delete[] Mv;
  delete[] Sv;
  
  delete[] TCriticVector;

  TFile ofile(file_out, "RECREATE");
  TTree * new_info_tree = info_tree -> CloneTree();

  ifile.Close();

  Block * block_out = new Block[dim_vector * 2];

  for(uint t = 0; t < tempstep; t++){

    TString treeName(TString("T_") + TString(std::to_string(t).c_str() ) );
    TString treeTitle(TString("TemperatureTree") + TString(std::to_string(t).c_str() ) );
    TTree * tree_out = new TTree(treeName, treeTitle);
    //tree_out -> SetDirectory(&ofile);

    for(uint i = 0; i < dim_vector; i++){
      TString t_index_m(TString("Lattice_mean_") + TString(std::to_string(i).c_str() ) );
      TString t_index_s(TString("Lattice_std_")  + TString(std::to_string(i).c_str() ) );

      tree_out -> Branch(t_index_m, &block_out[i]);
      tree_out -> Branch(t_index_s, &block_out[dim_vector + i]);
    }
    for(uint i = 0; i < dim_vector; i++){
      block_out[i].setBlock(i, T[i][t], E[i][t], M[i][t], S[i][t], X[i][t], t);
      block_out[dim_vector + i].setBlock(i, dT[i][t], dE[i][t], dM[i][t], dS[i][t], X[i][t], t);
    }
    tree_out -> Fill();
  }

  for(uint i = 0; i < dim_vector; i++){
    delete[] T[i];
    delete[] E[i];
    delete[] M[i];
    delete[] S[i];
    delete[] dT[i];
    delete[] dE[i];
    delete[] dM[i];
    delete[] dS[i];
  }
  delete[] T;
  delete[] E;
  delete[] M;
  delete[] S;
  delete[] dT;
  delete[] dE;
  delete[] dM;
  delete[] dS;
  ofile.Write();
  ofile.Close();
}
/*
  void AnalysisLattice::run() {
  for(uint j = 0; j < iter + 1; j++){
  T_mean += T_vec[j];
  E_mean += ( (Block*) arr_in -> At(j) ) -> E;
  M_mean += ( (Block*) arr_in -> At(j) ) -> M;
  S_mean += ( (Block*) arr_in -> At(j) ) -> S;
  }
  T_mean /= (iter + 1);
  E_mean /= (iter + 1);
  M_mean /= (iter + 1);
  S_mean /= (iter + 1);
  for(uint j = 0; j < iter + 1; j++){
  T_std += pow(T_mean - ( (Block*) arr_in -> At(j) ) -> T, 2);
  E_std += pow(E_mean - ( (Block*) arr_in -> At(j) ) -> E, 2);
  S_std += pow(M_mean - ( (Block*) arr_in -> At(j) ) -> M, 2);
  M_std += pow(S_mean - ( (Block*) arr_in -> At(j) ) -> S, 2);
  }
  T_std = sqrt( T_std / iter );
  E_std = sqrt( E_std / iter );
  M_std = sqrt( M_std / iter );
  S_std = sqrt( S_std / iter );

  block_out = (Block *)arr_out -> ConstructedAt(t);
  block_out -> setBlock(i, T_mean, E_mean, M_mean, S_mean, t);

  block_out = (Block *)arr_out -> ConstructedAt(t + tempstep);
  block_out -> setBlock(i, T_std,  E_std,  M_std,  S_std,  t + tempstep);
  }
  //filla tree_out
  tree_out -> GetBranch(TString("Lattice_") + i_index) -> Fill();
  }
  tree_out -> Fill();
  ifile.Close();
  ofile.Write();
  ofile.Close();
  }
*/

TGraphErrors * AnalysisLattice::drawLattice(cuint& lattice_number,
                                            cuint& x_axis,
                                            cuint& y_axis,
                                            cbool& abs){
  static INFO info;

  TFile f(file_out, "read");
  TTree * info_tree = (TTree*) f.Get("info");
  TBranch * info_branch = info_tree -> GetBranch("Info");
  info_branch -> SetAddress(&info._tempmin);
  info_tree -> GetEvent(0);
  double tempmin = info._tempmin;
  double tempmax = info._tempmax;
  uint N         = info._N;
  uint dim       = info._dim;
  uint dim_vector= info._dim_vector;
  uint I0        = info._I0;
  uint iter      = info._iter;
  uint tempstep  = info._tempstep;
  uint datime    = info._datime_t;
  double Tcritic;

  double* x = new double[tempstep];
  double* y = new double[tempstep];
  double* dx = new double[tempstep];
  double* dy = new double[tempstep];

  Block * block_mean = 0;
  Block * block_std = 0;

  for(uint t = 0; t < tempstep; t++){
    TTree * tree = (TTree*) f.Get(TString("T_") + TString(std::to_string(t).c_str() ) );

    tree -> SetBranchAddress(TString("Lattice_mean_") + TString(std::to_string(lattice_number).c_str()), &block_mean );
    tree -> SetBranchAddress(TString("Lattice_std_") + TString(std::to_string(lattice_number).c_str()), &block_std);
    tree -> GetEntry(0);
    switch(x_axis) {
    case 1:{
      abs ? TMath::Abs( x[t] = block_mean -> E ) : x[t]  = block_mean -> E;
      dx[t] = block_std  -> E;
      break;
    }
    case 2:{
      x[t]  = block_mean -> T;
      dx[t] = block_std  -> T;
      break;
    }
    case 3:{
      x[t]  = block_mean -> M;
      dx[t] = block_std  -> M;
      break;
    }
    case 4:{
      x[t]  = block_mean -> S;
      dx[t] = block_std  -> S;
      break;
    }
    case 5:{
      x[t]  = block_mean -> X;
      dx[t] = block_std  -> X;
      break;
    }
    default:{
      x[t]  = block_mean -> T;
      dx[t] = block_std  -> T;
      std::cout << "Falling to default for x axis (Temperature)\n" << std::flush;
      break;
    }
    }
    switch(y_axis){
    case 1:{
      y[t]  = block_mean -> E;
      dy[t] = block_std  -> E;
      break;
    }
    case 2:{
      y[t]  = block_mean -> T;
      dy[t] = block_std  -> T;
      break;
    }
    case 3:{
      abs ? y[t]  = TMath::Abs( block_mean -> M ) : y[t] = block_mean -> M;
      dy[t] = block_std  -> M;
      break;
    }
    case 4:{
      y[t]  = block_mean -> S;
      dy[t] = block_std  -> S;
      break;
    }
    case 5:{
      y[t]  = block_mean -> X;
      dy[t] = block_std  -> X;
      break;
    }
    default:{
      y[t]  = block_mean -> M;
      dy[t] = block_std  -> M;
      std::cout << "Falling to default for y axis (Magnetization)\n" << std::flush;
      break;
    }
    }
  }
  
  TGraphErrors * graph = new TGraphErrors(tempstep, x, y, dx, dy);
  //graph->Draw();

  delete[] x;
  delete[] y;
  delete[] dx;
  delete[] dy;
  return graph;
}

TMultiGraph * AnalysisLattice::draw(cuint& x_axis, cuint& y_axis, cbool& abs){

  static INFO info;

  TFile ifile(file_out, "read");
  TTree * info_tree = (TTree*) ifile.Get("info");
  TBranch * info_branch = info_tree -> GetBranch("Info");
  info_branch -> SetAddress(&info._tempmin);
  info_tree -> GetEvent(0);
  double tempmin = info._tempmin;
  double tempmax = info._tempmax;
  uint N         = info._N;
  uint dim       = info._dim;
  uint dim_vector= info._dim_vector;
  uint I0        = info._I0;
  uint iter      = info._iter;
  uint tempstep  = info._tempstep;
  uint datime    = info._datime_t;

  TMultiGraph * multi = new TMultiGraph("multigraph", "Multigraph");
  TGraphErrors * gr;
  for(uint i = 0; i < dim_vector; i++){
    gr = drawLattice(i, x_axis, y_axis, abs);
    multi -> Add(gr);
  }
  multi -> Draw("AP");

  return multi;


}
