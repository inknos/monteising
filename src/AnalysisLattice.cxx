#include "AnalysisLattice.h"
#include "Block.h"
#include "INFO.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TH1D.h"
#include "TLegend.h"

#include <string>
#include <vector>

#include <iostream>
using namespace std;

//setting TempCritic
double AnalysisLattice::TempCritic = 2.27;

AnalysisLattice::AnalysisLattice(const TString& file_input, const TString& file_output) :
  file_in(file_input), file_out(file_output) {}

TString AnalysisLattice::getFileIn() const { return file_in; }

TString AnalysisLattice::getFileOut() const { return file_out; }

void AnalysisLattice::setFileIn(const TString& file_input) { file_in = file_input; }

void AnalysisLattice::setFileOut(const TString& file_output) { file_out = file_output; }

double AnalysisLattice::getTempCritic() { return AnalysisLattice::TempCritic; }

void AnalysisLattice::setTempCritic(const double& _TempCritic){ AnalysisLattice::TempCritic = _TempCritic; }

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

  std::cout << R"(
?=======================================?
!                                       !
!           Analysis Started            !
!                                       !
?=======================================?
!                                       !)";
  std::cout << endl
            << "! " << dim_vector << " lattices " << dim << "x" << dim << " with size " << N << " created\t!\n"
            << "! " << I0 << " iter. not collecting data\t!\n"
            << "! " << iter << " iter. collecting data  \t!\n"
            << "! from " << tempmin << " to " << tempmax << " in " << tempstep << "steps\t\t!\n"
            << R"(!                                       !
?=======================================?)"
            << endl;

  uint num_bin = 8;
  uint dim_before_bin = iter + 1;
  uint dim_after_bin = (uint) ( dim_before_bin / TMath::Power((int)2, (int)num_bin) ) ;

  std::cout << dim_before_bin << " " << dim_after_bin << std::endl << std::flush;

  std::vector<std::vector<double> > Tw ( dim_vector, std::vector<double>(dim_before_bin, 0.) );
  std::vector<std::vector<double> > Ew ( dim_vector, std::vector<double>(dim_before_bin, 0.) );
  std::vector<std::vector<double> > Mw ( dim_vector, std::vector<double>(dim_before_bin, 0.) );
  std::vector<std::vector<double> > Sw ( dim_vector, std::vector<double>(dim_before_bin, 0.) );

  std::vector<std::vector<double> > Tv ( dim_vector, std::vector<double>(dim_after_bin, 0.) );
  std::vector<std::vector<double> > Ev ( dim_vector, std::vector<double>(dim_after_bin, 0.) );
  std::vector<std::vector<double> > Mv ( dim_vector, std::vector<double>(dim_after_bin, 0.) );
  std::vector<std::vector<double> > Sv ( dim_vector, std::vector<double>(dim_after_bin, 0.) );

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
      T[i][t]  = 0 ; E[i][t] = 0 ; M[i][t] = 0 ; S[i][t] = 0 ; X[i][t] = 0;
      dT[i][t] = 0; dE[i][t] = 0; dM[i][t] = 0; dS[i][t] = 0; dX[i][t] = 0;
    }
  }

  double * TCriticVector =  new double[dim_vector];

  double * Tml  = new double[tempstep];
  double * Eml  = new double[tempstep];
  double * Mml  = new double[tempstep];
  double * Sml  = new double[tempstep];
  double * Xml  = new double[tempstep];
  double * dTml = new double[tempstep];
  double * dEml = new double[tempstep];
  double * dMml = new double[tempstep];
  double * dSml = new double[tempstep];
  double * dXml = new double[tempstep];
  for(uint t = 0; t < tempstep; t++) {
    Tml[t] = 0.;Eml[t] = 0.;Mml[t] = 0.;Sml[t] = 0.;Xml[t] = 0.;
    dTml[t] = 0.;dEml[t] = 0.;dMml[t] = 0.;dSml[t] = 0.;dXml[t] = 0.;
  }

  Block * block = 0;
  TTree * tree;
  for(uint t = 0; t < tempstep; t++) {
    //cout << "t " << t << endl << flush;
    TString treeName(TString("T_") + TString(std::to_string(t).c_str() ) );
    tree = (TTree*) ifile.Get(treeName);

    for(uint i = 0; i < dim_vector; i++) {
      //cout << "i " << i << endl << flush;
      TString branchName(TString("Lattice_") + std::to_string(i).c_str() );
      tree -> SetBranchAddress(branchName, &block);
      // X = <M2> - <M>2 / T
      for(uint j = 0; j < dim_before_bin; j++) {
        tree -> GetEntry(j);
        Tw[i][j] = block -> T;
        Ew[i][j] = block -> E;
        Mw[i][j] = block -> M;
        Sw[i][j] = block -> S;
        //std::cout << Tw[i][j] << std::endl << std::flush;
      }
      Tv[i] = binN(num_bin, Tw[i]);
      Ev[i] = binN(num_bin, Ew[i]);
      Mv[i] = binN(num_bin, Mw[i]);
      Sv[i] = binN(num_bin, Sw[i]);
    }

    for(uint i = 0; i < dim_vector; i++){
      //std::cout << "i " << i << std::endl << std::flush;
      for(uint j = 0; j < dim_after_bin; j++){
        //std::cout << "j " << j << std::endl << std::flush;
        T[i][t] += Tv[i][j];
        E[i][t] += Ev[i][j];
        M[i][t] += Mv[i][j];
        S[i][t] += Sv[i][j];
        //std::cout << Tv[i][j] << std::endl << std::flush;
      }
      T[i][t] /= (dim_after_bin);
      E[i][t] /= (dim_after_bin);
      M[i][t] /= (dim_after_bin);
      S[i][t] /= (dim_after_bin);
      //std::cout << T[i][t] << std::endl << std::flush;
      for(uint j = 0; j < dim_after_bin; j++){
        //std::cout << "j " << j << std::endl << std::flush;
        dT[i][t] += ( T[i][t] - Tv[i][j] ) * ( T[i][t] - Tv[i][j] );
        dE[i][t] += ( E[i][t] - Ev[i][j] ) * ( E[i][t] - Ev[i][j] );
        dM[i][t] += ( M[i][t] - Mv[i][j] ) * ( M[i][t] - Mv[i][j] );
        dS[i][t] += ( S[i][t] - Sv[i][j] ) * ( S[i][t] - Sv[i][j] );
      }

      dT[i][t] = TMath::Sqrt( dT[i][t] / ( dim_after_bin ) );
      dE[i][t] = TMath::Sqrt( dE[i][t] / ( dim_after_bin ) );
      dS[i][t] = TMath::Sqrt( dS[i][t] / ( dim_after_bin ) );
      //
      //calcolo suscettivita
      X[i][t] = dM[i][t];
      //calcolo errore suscettività
      for(uint j = 0; j < dim_after_bin; j++){
        dX[i][t] += X[i][t] - ( M[i][t] - Mv[i][j] ) * ( M[i][t] - Mv[i][j] );
        dX[i][t] *= dX[i][t];
      }
      dX[i][t] = TMath::Sqrt( dX[i][t] / ( dim_after_bin ) );
      X[i][t] /= T[i][t] ;
      dX[i][t] /= T[i][t] ;
      //
      dM[i][t] = TMath::Sqrt( dS[i][t] / ( dim_after_bin ) );

      Tml[t] += TMath::Abs( T[i][t] );
      Eml[t] += TMath::Abs( E[i][t] );
      Mml[t] += TMath::Abs( M[i][t] );
      Sml[t] += TMath::Abs( S[i][t] );
      Xml[t] += TMath::Abs( X[i][t] );
    }
    //std::cout << "step1" << std::endl << std::flush;
    //
    Tml[t] /= dim_vector;
    Eml[t] /= dim_vector;
    Mml[t] /= dim_vector;
    Sml[t] /= dim_vector;
    Xml[t] /= dim_vector;
    for(uint i = 0 ; i < dim_vector ; i++){
      dTml[t] += ( Tml[t] - TMath::Abs(T[i][t]) )*( Tml[t] - TMath::Abs( T[i][t]) );
      dEml[t] += ( Eml[t] - TMath::Abs(E[i][t]) )*( Eml[t] - TMath::Abs( E[i][t]) );
      dMml[t] += ( Mml[t] - TMath::Abs(M[i][t]) )*( Mml[t] - TMath::Abs( M[i][t]) );
      dSml[t] += ( Sml[t] - TMath::Abs(S[i][t]) )*( Sml[t] - TMath::Abs( S[i][t]) );
      dXml[t] += ( Xml[t] - TMath::Abs(X[i][t]) )*( Xml[t] - TMath::Abs( X[i][t]) );
    }
    dTml[t] = TMath::Sqrt( dTml[t]/dim_vector );
    dEml[t] = TMath::Sqrt( dEml[t]/dim_vector );
    dMml[t] = TMath::Sqrt( dMml[t]/dim_vector );
    dSml[t] = TMath::Sqrt( dSml[t]/dim_vector );
    dXml[t] = TMath::Sqrt( dXml[t]/dim_vector );
    //
    std::cout << "[ done "<< (int) (((double) t / tempstep) * 100)  <<"% ] stepT :  " << t << std::endl << std::flush;
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

  //std::cout << "pre analysis complete \n" << std::flush;

  TFile ofile(file_out, "RECREATE");
  TTree * new_info_tree = info_tree -> CloneTree();

  ifile.Close();

  Block * block_out = new Block[dim_vector * 2];
  Block * block_ml = new Block[tempstep * 2];

  for(uint t = 0; t < tempstep; t++){

    TString treeName(TString("T_") + TString(std::to_string(t).c_str() ) );
    TString treeTitle(TString("TemperatureTree") + TString(std::to_string(t).c_str() ) );
    TTree * tree_out = new TTree(treeName, treeTitle);

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

    TString ml_index_m(TString("Lattice_ml_") + TString(std::to_string(t).c_str() ) );
    TString ml_index_s(TString("Lattice_ml_std_")  + TString(std::to_string(t).c_str() ) );
    tree_out -> Branch(ml_index_m, &block_ml[t]);
    tree_out -> Branch(ml_index_s, &block_ml[tempstep + t]);
    block_ml[t].setBlock(t, Tml[t], Eml[t], Mml[t], Sml[t], Xml[t], t);
    block_ml[tempstep + t].setBlock(t, dTml[t], dEml[t], dMml[t], dSml[t], dXml[t], t);

    tree_out -> Fill();
  }
  std::cout << "end analysis\n" << std::flush;

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
  //
  delete[] TCriticVector;
  delete[] Tml;
  delete[] Eml;
  delete[] Mml;
  delete[] Sml;
  delete[] Xml;
  delete[] dTml;
  delete[] dEml;
  delete[] dMml;
  delete[] dSml;
  delete[] dXml;
  //
  ofile.Write();
  ofile.Close();
  std::cout << R"(?=======================================?
!                 end                   !
?=======================================?)" << std::endl;
}

TGraphErrors * AnalysisLattice::drawLattice(cuint& lattice_number,
                                            cuint& x_axis,
                                            cuint& y_axis) {
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

  double* x  = new double[tempstep];
  double* y  = new double[tempstep];
  double* dx = new double[tempstep];
  double* dy = new double[tempstep];

  Block * block_mean = 0;
  Block * block_std = 0;

  TString title;
  TString xtitle;
  TString ytitle;

  for(uint t = 0; t < tempstep; t++){

    TTree * tree = (TTree*) f.Get(TString("T_") + TString(std::to_string(t).c_str() ) );
    tree -> SetBranchAddress(TString("Lattice_mean_") + TString(std::to_string(lattice_number).c_str()), &block_mean );
    tree -> SetBranchAddress(TString("Lattice_std_") + TString(std::to_string(lattice_number).c_str()), &block_std);
    tree -> GetEntry(0);
    
    switch(y_axis){
    case 1:{
      y[t]  = block_mean -> E;
      dy[t] = block_std  -> E;
      ytitle = TString("Energy (J)");
      title  = TString("Energy on ");
      break;
    }
    case 2:{
      y[t]  = block_mean -> T;
      dy[t] = block_std  -> T;
      ytitle = TString("Temperature (#frac{J}{k_{B}})");
      title  = TString("Temperature on ");
      break;
    }
    case 3:{
      //abs ? y[t]  = TMath::Abs( block_mean -> M ) : y[t] = block_mean -> M;
      y[t]  = block_mean -> M;
      dy[t] = block_std  -> M;
      ytitle = TString("Magnetization per Site (#mu)");
      title  = TString("Magnetization per Site on ");
      break;
    }
    case 4:{
      y[t]  = block_mean -> S;
      dy[t] = block_std  -> S;
      ytitle = TString("Energy per Site (J)");
      title  = TString("Energy per Site on ");
      break;
    }
    case 5:{
      y[t]  = block_mean -> X;
      dy[t] = block_std  -> X;
      ytitle = TString("Magnetic Susceptibility (#frac{#mu}{k_{B}})");
      title  = TString("Magnetic Susceptibility on ");
      break;
    }
    default:{
      y[t]  = block_mean -> M;
      dy[t] = block_std  -> M;
      std::cout << "Falling to default for y axis (Magnetization)\n" << std::flush;
      ytitle = TString("Magnetization per Site (#mu)");
      title  = TString("Magnetization per Site on ");
      break;
    }
    }
    switch(x_axis) {
    case 1:{
      //abs ? TMath::Abs( x[t] = block_mean -> E ) : x[t]  = block_mean -> E;
      x[t]  = block_mean -> E;
      dx[t] = block_std  -> E;
      xtitle = TString("Energy (J)");
      title += TString("Energy");
      break;
    }
    case 2:{
      x[t]  = block_mean -> T;
      dx[t] = block_std  -> T;
      xtitle = TString("Temperature (#frac{J}{k_{B}})");
      title += TString("Temperature");
      break;
    }
    case 3:{
      x[t]  = block_mean -> M;
      dx[t] = block_std  -> M;
      xtitle = TString("Magnetization (#mu)");
      title += TString("Magnetization");
      break;
    }
    case 4:{
      x[t]  = block_mean -> S;
      dx[t] = block_std  -> S;
      xtitle = TString("Energy per Site (J)");
      title += TString("Energy per Site");
      break;
    }
    case 5:{
      x[t]  = block_mean -> X;
      dx[t] = block_std  -> X;
      xtitle = TString("Magnetic Susceptibility (#frac{#mu}{k_{B}})");
      title += TString("Magnetic Susceptibility");
      break;
    }
    default:{
      x[t]  = block_mean -> T;
      dx[t] = block_std  -> T;
      xtitle = TString("Temperature (#frac{J}{k_{B}})");
      title += TString("Temperature on ");
      std::cout << "Falling to default for x axis (Temperature)\n" << std::flush;
      break;
    }
    }
    //cout << "y[" << t << "] " << y[t] << endl << flush;
    //cout << "\t" << "\t" << "dy[" << t << "] " << dy[t] << endl << flush;
  }

  title += TString(" of Lattice ") + std::to_string(lattice_number).c_str();

  TGraphErrors * graph = new TGraphErrors(tempstep, x, y, dx, dy);
  graph -> SetTitle(title);
  graph -> GetXaxis() -> SetTitle(xtitle);
  graph -> GetYaxis() -> SetTitle(ytitle);
  graph -> SetMarkerStyle(22);
  graph -> SetMarkerColor(kBlue + 3);
  //graph->Draw();

  delete[] x;
  delete[] y;
  delete[] dx;
  delete[] dy;
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
  multi -> Draw("AP");

  return multi;
}

TGraphErrors * AnalysisLattice::drawLatticeMean(cuint& x_axis,
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
  double Tcritic;

  double* x  = new double[tempstep];
  double* y  = new double[tempstep];
  double* dx = new double[tempstep];
  double* dy = new double[tempstep];

  Block * block_mean = 0;
  Block * block_std = 0;

  TString title;
  TString xtitle;
  TString ytitle;

  for(uint t = 0; t < tempstep; t++){

    TTree * tree = (TTree*) f.Get(TString("T_") + TString(std::to_string(t).c_str() ) );
    tree -> SetBranchAddress(TString("Lattice_ml_") + TString(std::to_string(t).c_str()), &block_mean );
    tree -> SetBranchAddress(TString("Lattice_ml_std_") + TString(std::to_string(t).c_str()), &block_std);
    tree -> GetEntry(0);
    xtitle = "";
    ytitle = "";
    switch(y_axis){
    case 1:{
      y[t]  = block_mean -> E;
      dy[t] = block_std  -> E;
      ytitle = TString("Energy (J)");
      title = TString("Energy on ");
      break;
    }
    case 2:{
      y[t]  = block_mean -> T;
      dy[t] = block_std  -> T;
      ytitle = TString("Temperature (#frac{J}{k_{B}})");
      title = TString("Temperature on ");
      break;
    }
    case 3:{
      //abs ? y[t]  = TMath::Abs( block_mean -> M ) : y[t] = block_mean -> M;
      y[t]  = block_mean -> M;
      dy[t] = block_std  -> M;
      ytitle = TString("Magnetization per Site (#mu)");
      title = TString("Magnetization per Site on ");
      break;
    }
    case 4:{
      y[t]  = block_mean -> S;
      dy[t] = block_std  -> S;
      ytitle = TString("Energy per Site (J)");
      title = TString("Energy per Site on ");
      break;
    }
    case 5:{
      y[t]  = block_mean -> X;
      dy[t] = block_std  -> X;
      ytitle = TString("Magnetic Susceptibility (#frac{#mu}{k_{B}})");
      title = TString("Magnetic Susceptibility on ");
      break;
    }
    default:{
      y[t]  = block_mean -> M;
      dy[t] = block_std  -> M;
      std::cout << "Falling to default for y axis (Magnetization)\n" << std::flush;
      ytitle = TString("Magnetization per Site (#mu)");
      title = TString("Magnetization per Site on ");
      break;
    }
    }
    switch(x_axis) {
    case 1:{
      //abs ? TMath::Abs( x[t] = block_mean -> E ) : x[t]  = block_mean -> E;
      x[t]  = block_mean -> E;
      dx[t] = block_std  -> E;
      xtitle = TString("Energy (J)");
      title += TString("Energy");
      break;
    }
    case 2:{
      x[t]  = block_mean -> T;
      dx[t] = block_std  -> T;
      xtitle = TString("Temperature (#frac{J}{k_{B}})");
      title += TString("Temperature");
      break;
    }
    case 3:{
      x[t]  = block_mean -> M;
      dx[t] = block_std  -> M;
      xtitle = TString("Magnetization (#mu)");
      title += TString("Magnetization");
      break;
    }
    case 4:{
      x[t]  = block_mean -> S;
      dx[t] = block_std  -> S;
      xtitle = TString("Energy per Site (J)");
      title += TString("Energy per Site");
      break;
    }
    case 5:{
      x[t]  = block_mean -> X;
      dx[t] = block_std  -> X;
      xtitle = TString("Magnetic Susceptibility (#frac{#mu}{k_{B}})");
      title += TString("Magnetic Susceptibility");
      break;
    }
    default:{
      x[t]  = block_mean -> T;
      dx[t] = block_std  -> T;
      xtitle = TString("Temperature (#frac{J}{k_{B}})");
      title += TString("Temperature on ");
      std::cout << "Falling to default for x axis (Temperature)\n" << std::flush;
      break;
    }
    }

    //cout << "y[" << t << "] " << y[t] << endl << flush;
    //cout << "\t" << "\t" << "dy[" << t << "] " << dy[t] << endl << flush;
  }

  title += TString(": Mean On All Lattices");

  TGraphErrors * graph = new TGraphErrors(tempstep, x, y, dx, dy);
  graph -> SetTitle(title);
  graph -> GetXaxis() -> SetTitle(xtitle);
  graph -> GetYaxis() -> SetTitle(ytitle);
  graph -> SetMarkerStyle(22);
  graph -> SetMarkerColor(kBlue + 3);
  //graph->Draw();
  
  delete[] x;
  delete[] y;
  delete[] dx;
  delete[] dy;
  return graph;
}

////FITTING

////FITTING
double AnalysisLattice::analiticX(double * x, double * par){
  return par[2]*TMath::Power( TMath::Abs( (x[0] - par[0])/par[0]) , -par[1] );
}

double AnalysisLattice::analiticM(double * x , double * par ){
  return par[2]*TMath::Power( TMath::Abs( (x[0] - par[0])/par[0]) , par[1] );
}

void AnalysisLattice::fitLattice( bool mean,
                                  cuint& x_axis,
                                  cuint& y_axis,
                                  double fit_temp_min = 1.5 ,
                                  double fit_temp_max = 3.0 ,
                                  int lat_number = 0
                                  ) {
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
  ifile.Close();


  TGraphErrors * g;

  if(mean){
    g = drawLatticeMean(x_axis , y_axis);
  }
  else{
    g = drawLattice(lat_number, x_axis , y_axis );
  }

  findTcritic();

  TF1 * f;
  if(y_axis == 3){
    f = new TF1("f" , AnalysisLattice::analiticM , fit_temp_min , fit_temp_max , 3);
    f -> SetParameters(TempCritic, 0.12, 1.1);
    f -> SetParNames("T critic", "exp critic", "coeff");
    g -> SetTitle("Temperature/Magnetization Fit");
    g -> GetYaxis() -> SetTitle("Magnetization per Site (#mu)");
  }
  if(y_axis == 5){
    f = new TF1("f" , AnalysisLattice::analiticX , fit_temp_min , fit_temp_max , 3);
    f -> SetParameters(TempCritic, 1.175, 6e-6);
    f -> SetParNames("T critic", "exp critic", "coeff");
    g -> SetTitle("Temperature/Magnetic Susceptibility Fit");
    g -> GetYaxis() -> SetTitle("Magnetic Susceptibility (#frac{#mu}{k_{B}})");
  }
  g -> SetMarkerStyle(22);
  g -> SetMarkerColor(kBlue + 3);
  g -> GetXaxis() -> SetTitle("Temperature (#frac{J}{k_{B}})");
  g -> Fit(f , "R");
  g -> Draw("ALP");
}

/*
  void AnalysisLattice::plotAnalitic(){
  TF1 *fteo = new TF1("fteo", analiticM ,0.5,2.4,1);
  fteo->SetLineColor(kRed);
  double zeropam = 0.;
  double alpha = 0.125;
  fteo->SetParameters(alpha , zeropam);
  //fteo->SetParNames("normalizzazione","coefficiente");
  fteo->Draw();
  }
*/

std::vector<double> AnalysisLattice::bin(const std::vector<double> & vec){
  uint bin = 2;
  uint dimbin = (uint) (vec.size() / bin);
  std::vector<double> binvec;
  for(uint i = 0; i < dimbin; i++){
    binvec.push_back( ( vec[i * 2] + vec[i * 2 + 1] ) * 0.5 );
  }
  return binvec;
}

std::vector<double> AnalysisLattice::binN(cuint& num_bin, const std::vector<double>& vec){
  std::vector<double> binvec( vec );
  for(uint i = 0; i < num_bin; i++){
    binvec = bin(binvec);
    //std::cout << "bin" << std::endl << std::flush;
  }
  //std::cout << "dim " << binvec.size() << std::endl << std::flush;
  return binvec;
}


void AnalysisLattice::findTcritic(){
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

  double* Temp  = new double[tempstep];
  double* Susc  = new double[tempstep];

  Block * block_mean = 0;

  double SuscMax = 0.;
  int tMax = 0;

  for(uint t = 0; t < tempstep; t++){

    TTree * tree = (TTree*) f.Get(TString("T_") + TString(std::to_string(t).c_str() ) );
    tree -> SetBranchAddress(TString("Lattice_ml_") + TString(std::to_string(t).c_str()), &block_mean );
    tree -> GetEntry(0);

    Temp[t] = block_mean -> T;
    Susc[t] = block_mean -> X;

    if(Susc[t] > SuscMax){
      SuscMax = Susc[t];
      tMax = t;
    }

  }
  cout << "TempCritic is : " << Temp[tMax] << endl;
  setTempCritic(Temp[tMax]);

  delete[] Temp;
  delete[] Susc;
}

void AnalysisLattice::evalBinning(const uint& nb){
  static INFO info;
  TFile f(file_in, "read");
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
  uint num_bin = nb;
  uint dim_before_bin = iter + 1;
  uint dim_after_bin = (uint) ( dim_before_bin / TMath::Power((int)2, (int)0) ) ;
  std::vector<double> Eb ( dim_before_bin, 0. );
  std::vector<double> Ea ( dim_after_bin,  0. );
  std::vector<double> E  ( num_bin, 0. );
  std::vector<double> dE ( num_bin, 0. );
  std::vector<double> bin_vec (num_bin, 0.);
  std::vector<double> bin_v (num_bin, 0.);
  Block * block = 0;
  TString treeName(TString("T_") + TString(std::to_string(50).c_str() ) );
  TTree * tree = (TTree*) f.Get(treeName);
  TString branchName(TString("Lattice_") + std::to_string(0).c_str() );
  tree -> SetBranchAddress(branchName, &block);
  for(uint j = 0; j < dim_before_bin; j++) {
    tree -> GetEntry(j);
    Eb[j] = block -> E;
  } // get data
  for(uint b = 0; b < num_bin; b++) {
    bin_vec[b] = TMath::Power((int)2, (int)b);
    bin_v[b] = b;
    dim_after_bin = (uint) ( dim_before_bin / (int) bin_vec[b] ) ;
    Ea = std::vector<double>( dim_after_bin,  0. );
    Ea = binN(b, Eb);
    for(uint j = 0; j < dim_after_bin; j++){
      E[b] += Ea[j];
    }
    E[b] /= (dim_after_bin);
    for(uint j = 0; j < dim_after_bin; j++){
      dE[b] += ( E[b] - Ea[j] ) * ( E[b] - Ea[j] );
    }
    dE[b] = TMath::Sqrt( dE[b] / ( dim_after_bin ) );
  }
  std::reverse(bin_v.begin(), bin_v.end());
  TGraph * gr = new TGraph(num_bin, bin_v.data(), dE.data());
  gr -> SetMarkerStyle(33);
  gr -> SetMarkerSize(2);
  gr -> SetMarkerColor(kBlue + 3);
  gr -> SetTitle("Energy Std. Dev. of Lattice_0 at T_50 over Binning");
  gr->GetXaxis()->SetTitle("Bin Folding");
  gr->GetYaxis()->SetTitle("Energy Std. Dev.");
  gr -> Draw("ALP");
}
