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
    const uint N;         // number of spins along one edge 
    const uint dim;       // dimension
    const uint num_spin;  // total number of spins : N^dim
    bool * lattice;       // lattice array
    static double T;      // temperature 

  public:
 
    /* CONSTRUCTORS */
    
    Lattice();

    Lattice(const uint& _N , const uint& _dim);

    Lattice(const Lattice& obj);
    
    /* DESTRUCTOR */

    ~Lattice();
    
    /* PHYSICAL AND NUMERICAL FUNCTIONS */
    
    bool flipSpin(const uint& n);
    
    int dE(const uint& spin) const;
    
    int energy() const;

    float magnetization() const ;
    
    void cooling();

    void cooling(const uint& iter);

    double * coolingPar();

    /* OVERLOADED OPERATORS */
    
    Lattice& operator=(const Lattice& obj);
    
    friend std::ostream &operator<<(std::ostream& out, const Lattice& lat);

    bool operator==(const Lattice& obj);

    /* GETTERS AND SETTERS */
    
    uint getN() const;
    
    uint getDim() const;
    
    uint getNumSpin() const;

    bool getSpin(const uint & i) const;

    static double getT();
    
    static void setT(const double& _T);
    
    /* OTHERS */

    void printLatticeCSV(const TString& name) const;

    void printLatticeROOT(const TString& name , const TString& ln = "lat") const;


  ClassDef(Lattice,1)
    };

#endif
