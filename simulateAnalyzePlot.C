#include "SimulationLattice.h"
#include "AnalysisLattice.h"

#include "TStopwatch.h"


void simulateAnalyzePlot(int n = 10, int dim = 2, int dim_v = 10, const char* nn = "test_cooling.root",
			 int i0 = 10000, int iter = 100000, double t1 = 0., double t2 = 4., int step = 20 ){
  TStopwatch timer;

  SimulationLattice s(n, dim, dim_v, nn, i0, iter, t1, t2, step);

  timer.Start();
  s.run();
  timer.Stop();
  timer.Print();

  AnalysisLattice a(nn, "out.root");
  
  timer.Start();
  a.run();
  timer.Stop();
  timer.Print();

  timer.Start();
  a.draw(TEMPERATURE, MAGNETIZATION);
  timer.Stop();
  timer.Print();

}
