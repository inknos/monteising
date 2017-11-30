#include "AnalysisLattice.h"
#include "Block.h"
#include "INFO.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBranchClones.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TH1D.h"
#include "TVectorT.h"

#include <string>

#include <iostream>
using namespace std;

AnalysisLattice::AnalysisLattice(const TString& file_input, const TString& file_output) :
  file_in(file_input), file_out(file_output) {}

TString AnalysisLattice::getFileIn() const { return file_in; }

TString AnalysisLattice::getFileOut() const { return file_out; }

void AnalysisLattice::setFileIn(const TString& file_input) { file_in = file_input; }

void AnalysisLattice::setFileOut(const TString& file_output) { file_out = file_output; }

void AnalysisLattice::run(){
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

  double* Tv= new double[iter + 1];
  double* Ev= new double[iter + 1];
  double* Mv= new double[iter + 1];
  double* Sv= new double[iter + 1];
  
  double** T = new double*[dim_vector];
  double** E = new double*[dim_vector];
  double** M = new double*[dim_vector];
  double** S = new double*[dim_vector];
  double** dT= new double*[dim_vector];
  double** dE= new double*[dim_vector];
  double** dM= new double*[dim_vector];
  double** dS= new double*[dim_vector];
  for(uint i = 0; i < dim_vector; i++){
    T[i] = new double[tempstep];
    E[i] = new double[tempstep];
    M[i] = new double[tempstep];
    S[i] = new double[tempstep];
    dT[i]= new double[tempstep];
    dE[i]= new double[tempstep];
    dM[i]= new double[tempstep];
    dS[i]= new double[tempstep];
  }
 
  for(uint i = 0; i < dim_vector; i++){
    for(uint t = 0; t < tempstep; t++){
      T[i][t] = 0; E[i][t] = 0; M[i][t] = 0; S[i][t] = 0;
      dT[i][t] = 0; dE[i][t] = 0; dM[i][t] = 0; dS[i][t] = 0;
    }
  }
  
  TClonesArray * array = new TClonesArray("Block", iter + 1);
  
  TFile ofile(file_out, "RECREATE");
  TTree * new_info_tree = info_tree -> CloneTree();

  for(uint i = 0; i < dim_vector; i++){
    TString treeName(TString("Lattice_") + TString(std::to_string(i).c_str() ) );
    TTree * tree = (TTree*) ifile.Get(treeName);

    for(uint t = 0; t < tempstep; t++){
      TString branchName(TString("T_") + std::to_string(t).c_str() );
      TBranchClones * branch = (TBranchClones*) tree -> GetBranch(branchName);
      branch -> SetAddress(&array);
      tree -> GetEntry(0);

      for(uint j = 0; j < iter + 1; j++){
        Tv[j] = ((Block*) array -> At(j)) -> T;
        Ev[j] = ((Block*) array -> At(j)) -> E;
        Mv[j] = ((Block*) array -> At(j)) -> M;
        Sv[j] = ((Block*) array -> At(j)) -> S;

        T[i][t] += Tv[j];
        E[i][t] += Ev[j];
        M[i][t] += Mv[j];
        S[i][t] += Sv[j];
      }
      T[i][t] /= (iter + 1);
      E[i][t] /= (iter + 1);
      M[i][t] /= (iter + 1);
      S[i][t] /= (iter + 1);
      for(uint j = 0; j < iter + 1; j++){
        dT[i][t] += ( T[i][t] - Tv[j] ) * ( T[i][t] - Tv[j] );
        dE[i][t] += ( E[i][t] - Ev[j] ) * ( E[i][t] - Ev[j] );
        dM[i][t] += ( M[i][t] - Mv[j] ) * ( M[i][t] - Mv[j] );
        dS[i][t] += ( S[i][t] - Sv[j] ) * ( S[i][t] - Sv[j] );
      }
      dT[i][t] = TMath::Sqrt( dT[i][t] / ( iter + 1 ) );
      dE[i][t] = TMath::Sqrt( dE[i][t] / ( iter + 1 ) );
      dS[i][t] = TMath::Sqrt( dM[i][t] / ( iter + 1 ) );
      dM[i][t] = TMath::Sqrt( dS[i][t] / ( iter + 1 ) );
    }
    cout << "Ho analizzato il reticolo " << i << endl << flush;
  }
  delete[] Tv;
  delete[] Ev;
  delete[] Mv;
  delete[] Sv;
  
  ifile.Close();
  
  TClonesArray * arr_out = new TClonesArray("Block", 1);
  Block * block_out;
  for(uint i = 0; i < dim_vector; i++){

    TString treeName(TString("Lattice_") + TString(std::to_string(i).c_str() ) );
    TString treeTitle(TString("LatticeTree") + TString(std::to_string(i).c_str() ) );
    TTree * tree_out = new TTree(treeName, treeTitle);
    //tree_out -> SetDirectory(&ofile);
    
    for(uint t = 0; t < tempstep; t++){
      TString t_index_m(TString("Tmean_") + TString(std::to_string(t).c_str() ) );
      TString t_index_s(TString("Tstd_")  + TString(std::to_string(t).c_str() ) );

      tree_out -> Branch(t_index_m, &arr_out);
      block_out = (Block*) arr_out -> ConstructedAt(0);
      block_out -> setBlock(i, T[i][t], E[i][t], M[i][t], S[i][t], t);
      tree_out -> GetBranch(t_index_m) -> Fill();
      //arr_out -> Clear();
      
      tree_out -> Branch(t_index_s, &arr_out);
      block_out = (Block*) arr_out -> ConstructedAt(0);
      block_out -> setBlock(i, dT[i][t], dE[i][t], dM[i][t], dS[i][t], t);
      tree_out -> GetBranch(t_index_s) -> Fill();
      //arr_out -> Clear();
    }
    tree_out -> Fill();
    //arr_out -> Clear();
    //tree_out -> Write();
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
                                            cuint& y_axis){
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

  TTree * tree = (TTree*) f.Get(TString("Lattice_") + TString(std::to_string(lattice_number).c_str() ) );
  TClonesArray * array_mean = new TClonesArray("Block", 1);
  TClonesArray * array_std  = new TClonesArray("Block", 1);

  double* x = new double[tempstep];
  double* y = new double[tempstep];
  double* dx = new double[tempstep];
  double* dy = new double[tempstep];

  TBranchClones * branch_mean;
  TBranchClones * branch_std;
  for(uint i = 0; i < tempstep; i++){
    branch_mean = (TBranchClones*) ( tree -> GetBranch(TString("Tmean_") + TString(std::to_string(i).c_str()))  );
    branch_mean -> SetAddress(&array_mean);
    branch_std  = (TBranchClones*) ( tree -> GetBranch(TString("Tstd_") + TString(std::to_string(i).c_str()))  );
    branch_std  -> SetAddress(&array_std);
    tree -> GetEntry(0);
    switch(x_axis){
    case 1:{
      x[i]  = ((Block*) array_mean -> At(0)) -> E;
      dx[i] = ((Block*) array_std  -> At(0)) -> E;
      break;
    }
    case 2:{
      x[i]  = ((Block*) array_mean -> At(0)) -> T;
      dx[i] = ((Block*) array_std  -> At(0)) -> T;
      break;
    }
    case 3:{
      x[i]  = ((Block*) array_mean -> At(0)) -> M;
      dx[i] = ((Block*) array_std  -> At(0)) -> M;
      break;
    }
    case 4:{
      x[i]  = ((Block*) array_mean -> At(0)) -> S;
      dx[i] = ((Block*) array_std  -> At(0)) -> S;
      break;
    }
    default:{
      x[i]  = ((Block*) array_mean -> At(0)) -> M;
      dx[i] = ((Block*) array_std  -> At(0)) -> M;
      std::cout << "Falling to default for x axis (Magnetization)\n" << std::flush;
      break;
    }
    }
    switch(y_axis){
    case 1:{
      y[i]  = ((Block*) array_mean -> At(0)) -> E;
      dy[i] = ((Block*) array_std  -> At(0)) -> E;
      break;
    }
    case 2:{
      y[i]  = ((Block*) array_mean -> At(0)) -> T;
      dy[i] = ((Block*) array_std  -> At(0)) -> T;
      break;
    }
    case 3:{
      y[i]  = ((Block*) array_mean -> At(0)) -> M;
      dy[i] = ((Block*) array_std  -> At(0)) -> M;
      break;
    }
    case 4:{
      y[i]  = ((Block*) array_mean -> At(0)) -> S;
      dy[i] = ((Block*) array_std  -> At(0)) -> S;
      break;
    }
    default:{
      y[i]  = ((Block*) array_mean -> At(0)) -> M;
      dy[i] = ((Block*) array_std  -> At(0)) -> M;
      std::cout << "Falling to default for y axis (Temperature)\n" << std::flush;
      break;
    }
    }
  }

  TGraphErrors * graph = new TGraphErrors(tempstep, x, y, dx, dy);
  //graph->Draw();

  delete[] x;
  delete[] y;

  return graph;
}

TMultiGraph * AnalysisLattice::draw(cuint& x_axis, cuint& y_axis){

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
    gr = drawLattice(i, x_axis, y_axis);
    multi -> Add(gr);
  }
  multi -> Draw("ALP");

  return multi;


}
