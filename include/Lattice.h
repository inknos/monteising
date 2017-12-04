#ifndef LATTICE_H
#define LATTICE_H

#include "TRandom3.h"
#include "TString.h"
#include "TObject.h"

#include <math.h>
#include <iostream>

class Lattice : public TObject{
 private:
  const unsigned int dim;       // dimension of the lattice
  const unsigned int q;         //
  const unsigned int N;         // number of spins in one dimension
  const unsigned int num_spin;  // pow(N, dim) total number of spins
  bool * lattice;               // [num_spin] lattice
  static double T;              // Temperature

 public:
  /* Public Constructors */
  Lattice();

  Lattice(const unsigned int& _N , const unsigned int& _dim);

  Lattice(const Lattice& obj);

  ~Lattice();

  /* Public Operators */
  Lattice& operator=(const Lattice& obj);
  
  friend std::ostream &operator<<(std::ostream& out, const Lattice& lat);

  bool operator==(const Lattice& obj);

  /* Getters */
  int getDim() const;

  int getN() const;

  bool getSpin(const unsigned int & i) const;

  int getQ() const;

  int getNumSpin() const;

  static double getT();
  
  static void setT(const double& _T);

  bool flipSpin(const unsigned int& n);

  void printLatticeCSV(const TString& name) const;

  void printLatticeROOT(const TString& name , const TString& ln = "lat") const;

  /* Physical functions */
  int energy(const bool& p = false ) const;
  
  int dE(const unsigned int& spin, const bool& p = false) const;

  float magnetization() const ;
  
  /* numerical function cooling */
  void cooling();

  void cooling(const uint& iter);

  double * coolingPar();

  ClassDef(Lattice,1)
    };

#endif
