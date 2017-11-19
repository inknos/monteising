#ifndef SIMULATION_LATTICE_H
#define SIMULATION_LATTICE_H

#include "Lattice.h"
#include "TString.h"
#include "TFile.h"

typedef unsigned int uint;

class SimulationLattice{
  
  private : 
   
    Lattice * lattice_vector;

    const uint dim_vector;
 
  public :

    SimulationLattice(); 
     
    SimulationLattice(const Lattice& , const uint&);
                 
    SimulationLattice(const uint& , const uint& , const uint&);   
   
    SimulationLattice(const SimulationLattice&);
    
    SimulationLattice& operator=(const SimulationLattice& obj); 

    ~SimulationLattice();                            

    void simulation();

  ClassDef(SimulationLattice, 1)
};

#endif
