#ifndef LATTICE_H
#define LATTICE_H

#include "TRandom3.h"
#include "TString.h"
#include "TObject.h"

#include <math.h>
#include <iostream>

typedef unsigned int uint;

class Lattice : public TObject{
 private:
  const uint dim;       // dimension of the lattice
  const uint N;         // number of spins in one dimension
  const uint num_spin;  // pow(N, dim) total number of spins
  bool * lattice;               // [num_spin] lattice
  static double T;              // Temperature

 public:
  /* Public Constructors */
  Lattice();

  Lattice(const uint& _N , const uint& _dim);

  Lattice(const Lattice& obj);

  ~Lattice();

  /* Public Operators */
  Lattice& operator=(const Lattice& obj);
  
  friend std::ostream &operator<<(std::ostream& out, const Lattice& lat);

  bool operator==(const Lattice& obj);

  /* Getters */
  uint getDim() const;

  uint getN() const;

  bool getSpin(const uint & i) const;

  uint getNumSpin() const;

  static double getT();
  
  static void setT(const double& _T);

  bool flipSpin(const uint& n);

  void printLatticeCSV(const TString& name) const;

  void printLatticeROOT(const TString& name , const TString& ln = "lat") const;

  /* Physical functions */
  int energy() const;
  
  int dE(const uint& spin) const;

  float magnetization() const ;
  
  /* numerical function cooling */
  void cooling();

  void cooling(const uint& iter);

  double * coolingPar();

  ClassDef(Lattice,1)
    };

#endif
