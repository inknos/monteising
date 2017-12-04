#ifndef SIMULATION_LATTICE_H
#define SIMULATION_LATTICE_H

#include "Lattice.h"
#include "TString.h"
#include "TFile.h"
#include "Block.h"

typedef unsigned int uint;

class SimulationLattice{
  
  private : 
   
    Lattice * lattice_vector;
    
    const uint dim_vector;
    const uint N;
    const uint dim;

    TString file;
    uint iter;
    uint I0;
    double tempmin;
    double tempmax;
    uint tempstep;
    
  public :

    SimulationLattice(); 

    SimulationLattice(const uint& _N, const uint& _dim, const uint& _dim_vector);
    
    SimulationLattice(const uint& _N, const uint& _dim, const uint& _dim_vector,
                      const TString& _file, const uint& _i0, const uint& _iter,
                      const double& _tempmin, const double& _tempmax, const uint& _tempstep);   
   
    SimulationLattice(const Lattice& _lat , const uint& _dim_vector,
                      const TString& _file, const uint& _i0, const uint& _iter,
                      const double& _tempmin, const double& _tempmax, const uint& _tempstep); 

    SimulationLattice(const SimulationLattice& obj);
    
    SimulationLattice& operator=(const SimulationLattice& obj); 

    ~SimulationLattice();                            

    Lattice getLattice(const uint& i) const;

    uint getDimVector() const;

    uint getN() const;
    
    uint getDim() const;

    TString getFile() const;

    uint getIter() const;

    double getTempMin() const;

    double getTempMax() const;

    double getTempStep() const;

    uint getI0() const;
    
    static uint getT();

    void setFile(const TString& _file);

    void setIter(const uint& _iter);

    void setTempMin(const double& _tempmin);

    void setTempMax(const double& _tempmax);

    void setTempStep(const uint& _tempstep);

    void setI0(const uint& _i0);

    static void setT(const double& _T);

    void run();
    
  
    //void simulation(const TString&, const uint&, const double&, const double&, const uint&);

  ClassDef(SimulationLattice, 1)
};

#endif
