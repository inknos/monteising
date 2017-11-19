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
  static double T = 0;              // Temperature

 public:
  /* Public Constructors */
  Lattice();

  Lattice(const unsigned int& _N, const unsigned int& _dim, unsigned const int& _q);

  Lattice(const Lattice &obj);

  ~Lattice();

  /* Public Operators */
  Lattice& operator=(const Lattice& obj);

  friend std::ostream &operator<<(std::ostream &out, const Lattice &lat);

  bool operator==(const Lattice&);

  /* Getters */
  int getDim() const;

  int getN() const;

  bool getSpin(const unsigned int & i) const;

  int getQ() const;

  int getNumSpin() const;

  static double getT() const;
  void setT(double);



  bool flipSpin(const unsigned int&);

  void printLatticeCSV(const TString&) const;

  void printLatticeROOT(const TString&, const TString&) const;

  /* Physical functions */
  int energy(const bool&) const;

  int energy2(const bool&) const;

  int energyFede() const ;

  int energyParallel(int) const;
  
  int dE(const unsigned int&, const bool&) const;

  float magnetization() const ;
  
  /* numerical function cooling */
  void cooling(const float&);

  

  ClassDef(Lattice,1)
    };

#endif
