//
//  autoCorrelation.C
//  
//
//  Created by Federico on 06/12/17.
//

#include <stdio.h>
#include "Lattice.h"

using namespace std;

void binning();
double autoCorrelation(double);
void tauVec();



void binning(int bin){
  Lattice l(10, 2);
  Lattice::setT(1.7);
  l.cooling(10000);
  
  int steps = 1000; //mettiamo multipli di 8
  double *  E  = new double[steps];
  double * Eb2 = new double[steps/bin];
  for(unsigned int t = 0; t < steps; t++){
    E[t] = 0.;
  }
  for(unsigned int t = 0; t < steps/bin; t++){
    Eb2[t] = 0.;
  }
  
  E[0] = l.energy();
  for(unsigned int t = 1 ; t < steps ; t++){
    E[t] = E[t-1] + l.coolingPar()[1];
  }
  
  double Emean  = 0.;
  double Eqmean = 0.;
  double var    = 0.;
  for(unsigned int t = 0; t < steps; t++){
    Emean  += E[t];
    Eqmean += E[t]*E[t];
  }
  Emean  /= steps;
  Eqmean /= steps;
  var = Eqmean - Emean*Emean;
  
  for(unsigned int t = 0 ; t < steps/bin ; t++){
    for(unsigned int k = 0; k < bin; k++){
      Eb2[t] += E[ bin*t + k ];
    }
    Eb2[t] /= bin;
  }
  
  double Eb2mean  = 0.;
  double Eb2qmean = 0.;
  double varb2    = 0.;
  for(unsigned int t = 0; t < steps/bin; t++){
    Eb2mean  += Eb2[t];
    Eb2qmean += Eb2[t]*Eb2[t];
  }
  Eb2mean  /= steps;
  Eb2qmean /= steps;
  varb2 = Eb2qmean - Eb2mean*Eb2mean;
  
  cout << "bin : " << bin << "   var : " << varb2 << endl << flush;
  
  
  
}







double autoCorrelation(double T){
  Lattice l(10, 2);
  Lattice::setT(T);
  //cout << "T : " << Lattice::getT() << endl << flush;
  l.cooling(10000);
  
  int steps = 1000;
  double *   E  = new double[steps];
  for(unsigned int t = 0; t < steps; t++){
    E[t] = 0.;
  }
  E[0] = l.energy();
  //cout << "E[0] : " << E[0] << endl;
  
  for(unsigned int t = 1 ; t < steps ; t++){
    E[t] = E[t-1] + l.coolingPar()[1];
  }

  double Emean  = 0.;
  double E2mean = 0.;
  for(unsigned int t = 0; t < steps; t++){
    Emean  += E[t];
    E2mean += E[t]*E[t];
  }
  Emean  /= steps;
  E2mean /= steps;
  
  //cout << "Emean^2 : " <<  Emean*Emean << endl;
  //cout << "E2mean : " << E2mean << endl;
  
  double * corr = new double[ steps ];
  for(unsigned int t = 0; t < steps; t++){
    corr[t] = 0.;
  }
  double tau = 0.;
  for(unsigned int t = 1; t < steps; t++){
    for(unsigned int i=0; i < steps-t; i++){
      corr[t] += E[i] * E[i+t] - Emean*Emean;
    }
    corr[t] /= (steps-t) * (E2mean - Emean*Emean);
    tau += corr[t];
  }
  tau += 0.5;
  
  //cout << "tau : " << tau << endl;
  
  delete[] E;
  delete[] corr;
  
  return tau;
}

void tauVec(){
  for(double t = 0.2; t < 30; t+=0.2){
    cout << "T : " << t << "   Autocorrelation : " << autoCorrelation(t) << endl << flush;
  }
}

