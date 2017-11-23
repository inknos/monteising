#include "Lattice.h"
#include "SimulationLattice.h"
#include "AnalysisLattice.h"
#include "TStopwatch.h"

#include <iostream>
using namespace std;

AnalysisLattice analysis(){
  gRandom->SetSeed(10000);

  TStopwatch timer;

  timer.Start();
  Lattice lat(50,2);
  cout << "[ done ] Lattice\n" << flush;
  timer.Stop();
  timer.Print();

  timer.Start();
  SimulationLattice sim(lat, 50);
  cout << "[ done ] SimulationLattice\n" << flush;
  timer.Stop();
  timer.Print();

  timer.Start();
  sim.simulation("simulation.root", 1000000, 0, 5, 10);
  cout << "[ done ] simulation\n" << flush;
  timer.Stop();
  timer.Print();

  timer.Start();
  AnalysisLattice ana("simulation.root");
  cout << "[ done ] AnalysisLattice\n" << flush;
  timer.Stop();
  timer.Print();

  timer.Start();
  //ana.analysis(TEMPERATURE, ENERGY, "te", "T/E");
  ana.analysis(TEMPERATURE, MAGNETIZATION, NO_ERR, "tm", "T/M");
  //ana.analysis(TEMPERATURE, SITE_ENERGY, NO_ERR, "ts", "T/S; temperature; energy");
  //ana.analysis(ENERGY, MAGNETIZATION, "em", "E/M");
  timer.Stop();
  timer.Print();

  return ana;
}
