/*-------------------------------------------------------------------*
 *                                                                   *
 * This macro explains how to perform a complete simulation          *
 * using th libraries.                                               *
 *                                                                   *
 *-------------------------------------------------------------------*/
#include "SimulationLattice.h" // Lattice.h included
#include "AnalysisLattice.h"   //
#include "TString.h"           // Other useful root inclusions
#include "TStopwatch.h"        //

{
  const unsigned int lattice_size(10);
  const unsigned int lattice_dimension(2);
  // this way a 10x10 lattice is created
  const unsigned int number_of_lattices(5);
  // simulation performed simultaneously on 5 lattices
  TString simulation_file("simulation_file.root");
  TString analysis_file("analysis_file.root");
  unsigned int iter_pre_simulation(1000000);
  unsigned int iter_for_simulation(500000);
  double min_temperature(0.5);
  double max_temperature(3.5);
  double steps_of_temperature(30);

  TStopwatch timer;

  SimulationLattice s(lattice_size, lattice_dimension,
                      number_of_lattices, simulation_file,
                      iter_pre_simulation, iter_for_simulation,
                      min_temperature, max_temperature,
                      steps_of_temperature);
  /*
   * Alternative way to call it:
   *
   * SimulationLattice s(lattice_size, lattice_dimension, number_of_lattices);
   * s.setFile(simulation_file);
   * s.setI0(iter_pre_simulation);
   * ...
   * s.run();
   */
  timer.Start();   // Simulation started with parameters set
  s.run();         // simulation_file always recreated
  timer.Stop();    //
  timer.Print();   //

  AnalysisLattice a(simulation_file, analysis_file);
  /*
   * input and output file must be provided even
   * if the output file already exists.
   * The output file is recreated only by run()
   * but is called in reading mode by the other
   * data members
   */
  timer.Start();   // Analysis started: parameters got
  a.run();         // from input file.
  timer.Stop();    // Output file Recreated.
  timer.Print();   //

  /*
   * Once the analysis is performed:
   * a.draw(TEMPERATURE, MAGNETIZATION);
   *
   * defined name   defined name   defined value
   *
   * ENERGY                        1
   * TEMPERATURE    TEMP           2
   * MAGNETIZATION  MAG            3
   * SITE_ENERGY    SENERGY        4
   * SUSCEPTIBILITY SUSC           5
   */

}
