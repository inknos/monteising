#include "AnalysisLattice.h"
#include "Block.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
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

  typedef struct{
    double _tempmin;
    double _tempmax;
    uint _N;
    uint _dim;
    uint _dim_vector;
    uint _I0;
    uint _iter;
    uint _tempstep;
    uint _datime_t;
  }INFO;
  
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

  double * E_vec = new double[iter+1];
  double * M_vec = new double[iter+1];
  double * S_vec = new double[iter+1];
  double * T_vec = new double[iter+1];
  
  double E_mean = 0;
  double M_mean = 0;
  double S_mean = 0;
  double T_mean = 0;

  double E_std = 0;
  double M_std = 0;
  double S_std = 0;
  double T_std = 0;

  TString i_index;
  TString t_index;

  TFile ofile(file_out, "recreate");
  TClonesArray * arr_out = new TClonesArray("Block", 1);
  Block * block_out; 
  
  int counter = 0;
  
  for(uint i = 0; i < dim_vector; i++) {
    i_index = TString( std::to_string(i).c_str() );
    TTreeReader info_reader(TString("Lattice_") + i_index, &ifile);
    TTree * tree_out = new TTree(TString("Lattice_") + i_index, TString("LaticeTree_") + i_index);
    
    for(uint t = 0; t < tempstep; t++) {
      t_index = TString(std::to_string(t).c_str());
      info_reader.Restart();

      TTreeReaderArray<Block> ivalue_reader(info_reader, TString("T_") + t_index );
      
      while(info_reader.Next()) {
        counter = 0;
        for(const Block& value: ivalue_reader){
          T_vec[counter] = value.T;
          E_vec[counter] = value.E;
          M_vec[counter] = value.M;
          S_vec[counter] = value.S;
          counter++;
        }
      }
      /////
      T_mean = 0;
      E_mean = 0;
      M_mean = 0;
      S_mean = 0;

      T_std = 0;
      E_std = 0;
      M_std = 0;
      S_std = 0;
      /*
      TH1D Thist(TVectorT<double>(iter + 1, T_vec));
      TH1D Ehist(TVectorT<double>(iter + 1, E_vec));
      TH1D Mhist(TVectorT<double>(iter + 1, M_vec));
      TH1D Shist(TVectorT<double>(iter + 1, S_vec));

      T_mean = Thist.GetMean();
      E_mean = Ehist.GetMean();
      M_mean = Mhist.GetMean();
      S_mean = Shist.GetMean();

      T_std  = Thist.GetStdDev();
      E_std  = Ehist.GetStdDev(); 
      M_std  = Mhist.GetStdDev(); 
      S_std  = Shist.GetStdDev();
      */
      for(uint j = 0; j < iter + 1; j++) {
        T_mean += T_vec[j];
        E_mean += E_vec[j];
        M_mean += M_vec[j];
        S_mean += S_vec[j];
      }
      T_mean /= (iter + 1);
      E_mean /= (iter + 1);
      M_mean /= (iter + 1);
      S_mean /= (iter + 1);
      for(uint j = 0; j < iter + 1; j++){
        T_std += (T_mean - T_vec[j]) * (T_mean - T_vec[j]);
        E_std += (E_mean - E_vec[j]) * (E_mean - E_vec[j]);
        S_std += (M_mean - M_vec[j]) * (M_mean - M_vec[j]);
        M_std += (S_mean - S_vec[j]) * (S_mean - S_vec[j]);
      }
      T_std = TMath::Sqrt( T_std / (iter + 1));
      E_std = TMath::Sqrt( E_std / (iter + 1));
      M_std = TMath::Sqrt( M_std / (iter + 1));
      S_std = TMath::Sqrt( S_std / (iter + 1));
      
      tree_out -> Branch(TString("Tmean_") + t_index, &arr_out);
      block_out = (Block *)arr_out -> ConstructedAt(0);
      block_out -> setBlock(i, T_mean, E_mean, M_mean, S_mean, t);
      tree_out -> GetBranch(TString("Tmean_") + t_index) -> Fill();
     
      tree_out -> Branch(TString("Tstd_") + t_index, &arr_out);
      block_out = (Block *)arr_out -> ConstructedAt(0);
      block_out -> setBlock(i, T_std,  E_std,  M_std,  S_std,  t); 
      tree_out -> GetBranch(TString("Tstd_") + t_index) -> Fill();
     }
    tree_out -> Fill();
  }
  delete []E_vec;
  delete []M_vec;
  delete []S_vec;
  delete []T_vec;
  ifile.Close();
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
