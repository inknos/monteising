#include "Lattice.h"
#include "SimulationLattice.h"
#include "AnalysisLattice.h"
#include "TStopwatch.h"

#include <iostream>
using namespace std;

void analysis(){
  gRandom->SetSeed(10000);

  TSopwatch timer;

  timer.Start();
  Lattice lat(100,2);
  cout << "[ done ] Lattice\n" << flush;
  timer.Stop();
  timer.Print();

  timer.Start();
  SimulationLattice sim(lat, 10);
  cout << "[ done ] SimulationLattice\n" << flush;
  timer.Stop();
  timer.Print();

  timer.Start();
  sim.simulation("simulation.root", 10000, 0, 5, 20);
  cout << "[ done ] simulation\n" << flush;
  timer.Stop();
  timer.Print();

  timer.Start();
  AnalysisLattice ana("simulation.root");
  cout << "[ done ] AnalysisLattice\n" << flush;
  timer.Stop();
  timer.Print();

  timer.Start();
  ana.analysis(TEMPERATURE, ENERGY);
  ana.analysis(TEMPERATURE, MAGNETIZATION);
  ana.analysis(TEMPERATURE, SITE_ENERGY);
  timer.Stop();
  timer.Print();


  
}
