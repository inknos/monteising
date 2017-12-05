#ifndef SIMULATION_LATTICE_H
#define SIMULATION_LATTICE_H

#include "Lattice.h"
#include "TString.h"
#include "TFile.h"
#include "Block.h"

class SimulationLattice{
  
  private : 
   
    Lattice * lattice_vector;
    
    const uint N;
    const uint dim;
    const uint dim_vector;

    TString file;
    
    uint iter;
    uint I0;
    
    double tempmin;
    double tempmax;
    uint tempstep;
    
  public :

    /* CONSTRUCTORS */
  
    SimulationLattice(); 

    SimulationLattice(const uint& _N, const uint& _dim, const uint& _dim_vector);
    
    SimulationLattice(const uint& _N, const uint& _dim, const uint& _dim_vector,
                      const TString& _file, const uint& _i0, const uint& _iter,
                      const double& _tempmin, const double& _tempmax, const uint& _tempstep);   
   
    SimulationLattice(const Lattice& _lat , const uint& _dim_vector,
                      const TString& _file, const uint& _i0, const uint& _iter,
                      const double& _tempmin, const double& _tempmax, const uint& _tempstep); 

    SimulationLattice(const SimulationLattice& obj);
    
    /* ASSIGNMENT OPERATOR */
    
    SimulationLattice& operator=(const SimulationLattice& obj); 

    /* DESTRUCTOR */
    
    ~SimulationLattice();
    
    /* GETTERS AND SETTERS */                            

    Lattice getLattice(const uint& i) const;

    uint getN() const;
    
    uint getDim() const;
    
    uint getDimVector() const;

    TString getFile() const;

    uint getIter() const;
    
    uint getI0() const;

    double getTempMin() const;

    double getTempMax() const;

    double getTempStep() const;
    
    static uint getT();

    void setFile(const TString& _file);

    void setIter(const uint& _iter);
    
    void setI0(const uint& _i0);

    void setTempMin(const double& _tempmin);

    void setTempMax(const double& _tempmax);

    void setTempStep(const uint& _tempstep);
    
    static void setT(const double& _T);
    
    /* RUN FUNCTION */
    
    void run();

  ClassDef(SimulationLattice, 1)
};

#endif
