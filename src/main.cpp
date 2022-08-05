/** MPS-Moving Particle Semi-implicit
 * ***(c) 2008 USP/TPN, NYXknowledge
 *** Marcio Tsukamoto, Guilherme Rueda, Mariana Robortella, 
 *** Rubens A. Amaro Jr, Cheng Liang Yee, et al.
 **/

#include "struct.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <cstddef>
#include <omp.h>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	/// Declare structures (see file struct.h)
	particles part;			// Particles data
	simParameters sim;		// Simulation data
	matProperties *prop;	// Properties of materials (Materials represent fluid, fixed solid, free solid, ...)

	//////////////
	/// USER INPUT (Can be read from a input file)
	sim.nMat = 1;			// Number of materials, e.g., nMat = 7: fluid-0 -> m = 0, fluid-1 -> m = 1, ..., fixed-solid-1 -> m = 5, ..., free-solid-1 -> m = 7)
	sim.dim = 2;			// Dimension (2 or 3)
	sim.nPart = 1000;		// Number of particles (considering all materials)
	sim.nNeighMax = 64;		// Max number of neighbors
	sim.g[0] = 0.0;			// Gravity X
	sim.g[1] = -9.81;		// Gravity Y
	sim.g[2] = 0.0;			// Gravity Z
	sim.dPart = 0.1;		// Initial adjacent particle distance
	sim.reS = 2.1;			// Effective radius small
	sim.reL = 4.0;			// Effective radius large (2D: 4.0, 3D: 2.1)
	sim.dt = 0.0001;		// Actual time-step
	sim.tStart = 0.0; 		// Start time
	sim.tFinal = 1.0;		// End time
	sim.nIter = 10e10;		// Simulation maximum number of iterations (Only to avoid no-end simulations)

	initProperties(&sim, prop);

	/// USER INPUT
	//////////////
	sim.reS = sim.reS*sim->dPart;
	sim.reL = sim.reL*sim->dPart
	sim.reS2 = sim.reS*sim.reS;
	sim.reL2 = sim.reL*sim.reL;
	sim.t = sim.tStart;

	/// Alloc memory and set to zero
	part.Alloc(sim.nPart, sim.nNeighMax);

	fprintf(stderr, "\n*Initialization: Free Solid");
	initSolidMotion(&part, &sim, prop);
	initSolidCollision(&part, &sim, prop);
	fprintf(stderr, " done!");
	
	/// Find Neighbors
	fprintf(stderr, "\n*Initialization: Neighbor list");
	CreateNeighborList(&part, &sim, prop);
	fprintf(stderr, " done!");
	
	/// Particle Number Density (PND). Small and Large radius
	fprintf(stderr, "\n*Initialization: PND");
	CalcPND(&part.pndS, &part.pndS_mat, &part.nNeighS, &part.neighS, &part.dNeighS, &part.wS, sim.reS, &sim, prop, &part);
	CalcPND(&part.pndL, &part.nNeighL, &part.neighL, &part.dNeighL, &part.wL, sim.reL, &sim, prop);
	fprintf(stderr, " done!");

	///	Time main loop
	fprintf(stderr, "\n*Start main loop");
	for(long iter=0; iter<sim.nIter; iter++)
	{
		/// Find Neighbors
		CreateNeighborList(&part, &sim, prop);
		/// Instantaneous Particle Number Density (PND) - pnd_star
		CalcPND(&part.pndS, &part.pndS_mat, &part.nNeighS, &part.neighS, &part.dNeighS, &part.wS, sim.reS, &sim, prop, &part);
		/// Boundary Conditions (Internal of free surface particle)
		CheckBoundCond(&part, &sim, prop);
		/// Compute solid motion
		calcSolidMotion(&part, &sim, prop);
		/// Printing-out
		OutputData(&part, &sim, prop, iter);

		sim.t = sim.t + sim.dt;
		/// Already done?
		if(sim.t>=sim.tFinal)
			break;
	}

	/// Clear memory
	part.Free();

	fprintf(stderr, "Simulation done!");

	return 0;
}
