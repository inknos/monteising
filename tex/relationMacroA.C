/*----------------------------------------------------------------*
 *                                                                *
 * This macro explains how to construct handle and draw a Lattice *
 *                                                                *                                                                
 *----------------------------------------------------------------*/

#include "Lattice.h" 
#include "DrawLattice.h"

{ 

  //LATTICE CLASS

  Lattice::setT(0.5); //Set static member Temperature to 0.5 
                      //Default is T = 0.
  
  int N = 10;
  int dim = 2;
  Lattice lat(N,dim); //Construct a bi-dimensional Lattice 
                      //with 10 spins along each edge
  
  cout << "lat visualization :" << endl;
  cout << lat << endl; //Print Lattice spins on the console
                       //The lattice spins will be disposed randomly
    
  cout << "energy E of lat : " << lat.energy() << endl; 
  cout << "magnetization M of lat : " << lat.magnetization() << endl; 
  
  lat.cooling(); //Perform a Metropolis step
  
  double * data = new double[4];
  data = lat.coolingPar(); //perform a Metropolis step collecting data
  
  cout << "T =  " << data[0] << " : temperature" << endl; 
  cout << "dE = " << data[1] << " : energy variation" << endl;
  cout << "dM = " << data[2] << " : magnetization variation" << endl;
  cout << "dS = " << data[3] << " : energy per site variation" << endl;
  delete[] data;  
  
  /*
   *
   * According to Metropolis algortihm if dE <0 
   * there's a probability [1 - exp(-dE/T)]  (K=1) 
   * that the step will not change lat state with a spinFlip
   * in that case dE=0 ; dM=0 ; dS=0;
   *
   */
  
  
 
 
 
 
  
  const unsigned iter = 10000;
  lat.cooling(iter); // Perform 10000 Metropolis steps
  
  cout << lat << endl; //Print Lattice spins on the console
                       //The lattice will now be ordered 
  
  cout << "energy E of lat : " << lat.energy() << endl; 
  // Energy will be around -200
  cout << "magnetization M of lat : " << lat.magnetization() << endl;
  // Magnetization per site will be about 1 or -1
  
  /*
   *
   * If T is set to a value greater than the critic temperature 
   * (for 2D Ising Model Tc = 2.27)
   * and the same Macro is executed, 
   * lat will not thermalize. 
   *
   */
  
  //DRAWLATTICE
 
 
  lat = Lattice(10,2); 
  
  DrawLattice drawLat1(lat); // Create a DrawLattice object 
                             // lat is passed in order to draw it
 

  cout << "lat with dim = " << lat.getDim() << endl;
  cout << "and edge = " << lat.getN() << "will be drawn" << endl;
  
  drawLat1.draw(); // Draws the three-dimensional Lattice lat
   
  
  Lattice lat3D(20,3); // Construct a three-dimensional Lattice 
                       // with 20 spins along each edge
 
  lat = lat3D;   // Use assignment operator to reassign lat
  
  
  cout << "lat with dim = " << lat.getDim() << endl;
  cout << "and edge = " << lat.getN() << "will be drawn" << endl;
  

  DrawLattice drawLat2(lat);
  drawLat2.draw();
   
}
