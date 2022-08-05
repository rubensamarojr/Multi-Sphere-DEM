/** MPS-Moving Particle Semi-implicit
 * ***(c) 2008 USP/TPN, NYXknowledge
 *** Marcio Tsukamoto, Guilherme Rueda, Mariana Robortella, 
 *** Rubens A. Amaro Jr, Cheng Liang Yee, et al.
 **/

/// Include headers
#include "sctruct.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "src_solid_transform.h"

using std::cerr;
using std::endl;

/// Declare functions
int initSolidMotion(particles *part, simParameters *sim, matProperties *prop);
int initSolidCollision(particles *part, simParameters *sim, matProperties *prop);
int calcSolidMotion(particles *part, simParameters *sim, matProperties *prop);
int calcSolidBC(particles *part, simParameters *sim, matProperties *prop);
int IsGoodDivision(int a, int b);
int solidTimeStep(simParameters *sim, particles *part, matProperties *prop, const int m);
int collisionForceMoment(simParameters *sim, particles *part, matProperties *prop, const int m, double &next_dtSolidRigid, int &next_nStepsSolid);
int calcPressureIntegration(particles *part, matProperties *prop, const int dim, const int m, real *force, real *moment);
int addGravityForce(simParameters *sim, matProperties *prop, const int m, real *force);
int calc2dMotionEquation(particles *part, simParameters *sim, matProperties *prop, const int m, real *force, real *moment);
int calc3dMotionEquation(particles *part, simParameters *sim, matProperties *prop, const int m, real *force, real *moment);
int update2dPosVel(particles *part, simParameters *sim, matProperties *prop, double trans[3], double rotate[3], const int m);
int update3dPosVel(particles *part, simParameters *sim, matProperties *prop, double trans[3], double rotate[3], const int m);

/// Global variable
static int collision_ok=false;

/////////////////////////////////////////////////////////////////////////////////////
/// Initialization of rigid body mechanical properties and kinematic variables
int initSolidMotion(particles *part, simParameters *sim, matProperties *prop)
{
	/// Initial solid time step (may change if collision occurs)
	sim->dtSolidRigid = sim->dt; 	// sim->dtSolidRigid: Solid time step. sim->dt: Fluid time step
	sim->nStepsSolid = 1;			// Number of Subcycling steps >= 1
	
	/// Loop for all nMat materials, e.g.:
	/// fluid-0 -> m = 0, fluid-1 -> m = 1, ..., fixed-solid-1 -> m = 5, ..., free-solid-1 -> m = 7)
	for(int m=0; m<sim->nMat; m++)
	{
		/// Global and Local geometrical limits
		real minGlobal[3], maxGlobal[3], minLocal[3], maxLocal[3];

		/// Verify if "m" is a free solid and have particles, i.e., if the solid exists
		if((prop[m].type==solid) && (prop[m].motion==freeMot) && (prop[m].part.n>0))
		{
			/// Allocate memory for ".n" variables
			// prop[m].part.n: Number of wall particles belonging to solid "m"
			// prop[m].partDummy.n: Number of dummy particles belonging to solid "m"
			makeVector(&prop[m].solid.normalA, prop[m].part.n); 		// Normal vector at each solid wall particle
			makeVector(&prop[m].solid.dCG, prop[m].part.n);				// Distance between wall particle and solid CG (center of geometry)
			makeVector(&prop[m].solid.dCGDummy, prop[m].partDummy.n);	// Distance between dummy wall particle and solid CG
			makeRealArray(&prop[m].solid.areaA, prop[m].part.n);		// Area of one wall particle. 2D -> area = lo. 3D -> area = lo*lo.
			
			/// Loop for all dimensions 2 or 3 and set values to zero
			for(int k=0; k<sim->dim; k++)
			{
				/// Variables with respect to CG of the solid
				prop[m].solid.vel[k]=prop[m].solid.vel1[k]=prop[m].solid.vel2[k]=0.0;			// Translation velocity (vel1 and vel2 are auxiliar!)
				prop[m].solid.velAng[k]=prop[m].solid.velAng1[k]=prop[m].solid.velAng2[k]=0.0;	// Angular velocity (velAng1 and velAng2 are auxiliar!)
				prop[m].solid.f1[k]=prop[m].solid.f2[k]=0.0;			// Auxiliar forces
				prop[m].solid.m1[k]=prop[m].solid.m2[k]=0.0;			// Auxiliar moments (torque)
				
				prop[m].solid.collisionForce[k]=0.0;					// Collision/contact force (Local coordinate system)
				prop[m].solid.collisionMoment[k]=0.0;					// Collision/contact moment (torque) (Local coordinate system)
				prop[m].solid.hydroForce[k]=0.0;						// Fluid pressure on wall - force (Local coordinate system)
				prop[m].solid.hydroMoment[k]=0.0;						// Fluid pressure on wall - moment (torque) (Local coordinate system)

				prop[m].solid.globalCollisionForce[k] = 0.0;			// Collision/contact force (Global coordinate system)
				prop[m].solid.globalCollisionMoment[k] = 0.0;			// Collision/contact moment (torque) (Global coordinate system)
				prop[m].solid.globalHydroForce[k] = 0.0;				// Fluid pressure on wall - force (Local coordinate system)
				prop[m].solid.globalHydroMoment[k] = 0.0;				// Fluid pressure on wall - moment (torque) (Global coordinate system)
			}
			
			/// Initial quaternion (Only for 3D)
			if(sim->dim==3)
			{
				real rot[3], dTheta;

				/// Orientation of the solid at initial time t = 0
				dTheta=sqrt(prop[m].solid.orient[0]*prop[m].solid.orient[0]+
							prop[m].solid.orient[1]*prop[m].solid.orient[1]+
							prop[m].solid.orient[2]*prop[m].solid.orient[2]);

				/// Rotate limits if solid is rotated at initial time t = 0
				if(dTheta!=0)
				{
					rot[0]=prop[m].solid.orient[0]/dTheta;
					rot[1]=prop[m].solid.orient[1]/dTheta;
					rot[2]=prop[m].solid.orient[2]/dTheta;
				}
				else
				{
					rot[0]=0.0;
					rot[1]=0.0;
					rot[2]=0.0;
				}

				prop[m].solid.quat[0]=cos(dTheta/2);
				prop[m].solid.quat[1]=rot[0]*sin(dTheta/2);
				prop[m].solid.quat[2]=rot[1]*sin(dTheta/2);
				prop[m].solid.quat[3]=rot[2]*sin(dTheta/2);
			}

			/// min-max position limits in Global coordinate system
			/// Values of 1st particle
			minGlobal[0]=part->r.x[prop[m].part.v[1]];
			minGlobal[1]=part->r.y[prop[m].part.v[1]];
			minGlobal[2]=part->r.z[prop[m].part.v[1]];
			maxGlobal[0]=part->r.x[prop[m].part.v[1]];
			maxGlobal[1]=part->r.y[prop[m].part.v[1]];
			maxGlobal[2]=part->r.z[prop[m].part.v[1]];
			
			/// Transform min-max, Global->Local Coordinate System
			if (sim->dim == 2)
			{
				/********** USING ANGLE ZZ FOR 2D ************/
				///	cos()	 sin()
				///	-sin()	 cos()
				minLocal[0]= minGlobal[0]*cos(prop[m].solid.orient[2])+minGlobal[1]*sin(prop[m].solid.orient[2]);
				minLocal[1]=-minGlobal[0]*sin(prop[m].solid.orient[2])+minGlobal[1]*cos(prop[m].solid.orient[2]);
				maxLocal[0]= maxGlobal[0]*cos(prop[m].solid.orient[2])+maxGlobal[1]*sin(prop[m].solid.orient[2]);
				maxLocal[1]=-maxGlobal[0]*sin(prop[m].solid.orient[2])+maxGlobal[1]*cos(prop[m].solid.orient[2]);
			}
			else if (sim->dim == 3)
			{
				/********** USING QUATERNNION FOR 3D ************/
				minLocal[0]=minGlobal[0];
				minLocal[1]=minGlobal[1];
				minLocal[2]=minGlobal[2];
				tranformRotation3dGlobal2Local(minLocal, prop[m].solid.quat);
				maxLocal[0]=maxGlobal[0];
				maxLocal[1]=maxGlobal[1];
				maxLocal[2]=maxGlobal[2];
				tranformRotation3dGlobal2Local(maxLocal, prop[m].solid.quat);
			}

			/// Compute normals vectors and areas 2D and 3D
			int n=0;
			real neighDist=sim->dPart*1.74; // sqrt(3)
			/// Loop for wall particles of solid material "m"
			for(int i=0; i<prop[m].part.n; i++)
			{
				//	i:  Particle ID considering only the solid material m
				//	iPart: Particle ID considering all materials, i.e., fluid, fixed solid, free solid, ...
				int iPart=prop[m].part.v[i];
				real orient[3], rAux[3];
				
				/// MIN - MAX Global Coordinate System
				minGlobal[0]=std::min(minGlobal[0], part->r.x[iPart]);
				minGlobal[1]=std::min(minGlobal[1], part->r.y[iPart]);
				minGlobal[2]=std::min(minGlobal[2], part->r.z[iPart]);
				maxGlobal[0]=std::max(maxGlobal[0], part->r.x[iPart]);
				maxGlobal[1]=std::max(maxGlobal[1], part->r.y[iPart]);
				maxGlobal[2]=std::max(maxGlobal[2], part->r.z[iPart]);
				
				/// MIN - MAX Local Coordinate System
				if (sim->dim == 2)
				{
					/********** USING ANGLE ZZ FOR 2D ************/
					///	cos()	 sin()
					///	-sin()	 cos()
					rAux[0]= part->r.x[iPart]*cos(prop[m].solid.orient[2])+part->r.y[iPart]*sin(prop[m].solid.orient[2]);
					rAux[1]=-part->r.x[iPart]*sin(prop[m].solid.orient[2])+part->r.y[iPart]*cos(prop[m].solid.orient[2]);
				}
				else if (sim->dim == 3)
				{
					/********** USING QUATERNNION FOR 3D ************/
					rAux[0]= part->r.x[iPart];
					rAux[1]= part->r.y[iPart];
					rAux[2]= part->r.z[iPart];
					tranformRotation3dGlobal2Local(rAux, prop[m].solid.quat);
				}
				
				/// Local Coordinate System
				minLocal[0]=std::min(minLocal[0], rAux[0]);
				minLocal[1]=std::min(minLocal[1], rAux[1]);
				minLocal[2]=std::min(minLocal[2], rAux[2]);
				maxLocal[0]=std::max(maxLocal[0], rAux[0]);
				maxLocal[1]=std::max(maxLocal[1], rAux[1]);
				maxLocal[1]=std::max(maxLocal[2], rAux[1]);
				
				// Variable used to orient the normal vector at each wall particle
				orient[0]=0.0;
				orient[1]=0.0;
				orient[2]=0.0;
				
				/// Loop for wall particles of solid material "m"
				for(int j=0; j<prop[m].part.n; j++)
				{
					real vec[3];
					int jPart=prop[m].part.v[j];

					if(sim->dim==3)
						prop[m].solid.areaA.v[i]=sim->dPart*sim->dPart; // dPart is the inital distance between particles
					else if(sim->dim==2)
						prop[m].solid.areaA.v[i]=sim->dPart;
					if(iPart==jPart)
						continue;

					vec[0]=part->r.x[iPart]-part->r.x[jPart];
					vec[1]=part->r.y[iPart]-part->r.y[jPart];
					vec[2]=part->r.z[iPart]-part->r.z[jPart];
					if((fabs(vec[0])<neighDist) && (fabs(vec[1])<neighDist) && (fabs(vec[2])<neighDist))
					{
						orient[0]+=vec[0];
						orient[1]+=vec[1];
						orient[2]+=vec[2];
						n++;
					}
				}
				// Loop for dummy particles of solid material "m"
				for(int j=0; j<prop[m].partDummy.n; j++)
				{
					real vec[3];
					int jPart=prop[m].partDummy.v[j];

					if(iPart==jPart)
						continue;
					vec[0]=part->r.x[iPart]-part->r.x[jPart];
					vec[1]=part->r.y[iPart]-part->r.y[jPart];
					vec[2]=part->r.z[iPart]-part->r.z[jPart];
					if((fabs(vec[0])<neighDist) && (fabs(vec[1])<neighDist) && (fabs(vec[2])<neighDist))
					{
						orient[0]+=vec[0];
						orient[1]+=vec[1];
						orient[2]+=vec[2];
						n++;
					}
				}

				/// Normalizing normals in Global Coordinate System
				real normalMod;

				normalMod=sqrt(orient[0]*orient[0]+orient[1]*orient[1]+orient[2]*orient[2]);
				if(normalMod!=0)
				{
					prop[m].solid.normalA.x[i]=orient[0]/normalMod;
					prop[m].solid.normalA.y[i]=orient[1]/normalMod;
					prop[m].solid.normalA.z[i]=orient[2]/normalMod;
				}
				else
				{
					fprintf(stderr, "\nWarning: The norm of a normal vector is 0 ! \n\n");
					prop[m].solid.normalA.x[i]=0;
					prop[m].solid.normalA.y[i]=0;
					prop[m].solid.normalA.z[i]=0;
				}
				if((isnan(prop[m].solid.normalA.x[i])!=0) || (isnan(prop[m].solid.normalA.y[i])!=0) || (isnan(prop[m].solid.normalA.z[i])!=0))
				{
					fprintf(stderr, "\nWarning: normal is NaN 0 ! \n\n");
					prop[m].solid.normalA.x[i]=0;
					prop[m].solid.normalA.y[i]=0;
					prop[m].solid.normalA.z[i]=0;
				}
                 
				real nLocal[3];
				
				/// Transform Normal, Global->Local Coordinate System
				if (sim->dim == 2)
				{
					/********** USING ANGLE ZZ FOR 2D ************/
					///	cos()	 sin()
					///	-sin()	 cos()
					nLocal[0]= prop[m].solid.normalA.x[i]*cos(prop[m].solid.orient[2])+prop[m].solid.normalA.y[i]*sin(prop[m].solid.orient[2]);
					nLocal[1]=-prop[m].solid.normalA.x[i]*sin(prop[m].solid.orient[2])+prop[m].solid.normalA.y[i]*cos(prop[m].solid.orient[2]);
					prop[m].solid.normalA.x[i]=nLocal[0];
					prop[m].solid.normalA.y[i]=nLocal[1];
					part->normal.x[iPart]=prop[m].solid.normalA.x[i];
					part->normal.y[iPart]=prop[m].solid.normalA.y[i];
				}
				else if (sim->dim == 3)
				{
					/********** USING QUATERNNION FOR 3D ************/
					nLocal[0]=prop[m].solid.normalA.x[i];
					nLocal[1]=prop[m].solid.normalA.y[i];
					nLocal[2]=prop[m].solid.normalA.z[i];
					tranformRotation3dGlobal2Local(nLocal, prop[m].solid.quat);
					prop[m].solid.normalA.x[i]=nLocal[0];
					prop[m].solid.normalA.y[i]=nLocal[1];
					prop[m].solid.normalA.z[i]=nLocal[2];
					part->normal.x[iPart]=prop[m].solid.normalA.x[i];
					part->normal.y[iPart]=prop[m].solid.normalA.y[i];
					part->normal.z[iPart]=prop[m].solid.normalA.z[i];
				}
			}
			

			/// External cte force (only for forced solid)
			prop[m].solid.constantLoadPos[0]=0.0;
			prop[m].solid.constantLoadPos[1]=0.0;
			prop[m].solid.constantLoadPos[2]=0.0;
			prop[m].solid.constantLoad[0]=0.0;
			prop[m].solid.constantLoad[1]=0.0;
			prop[m].solid.constantLoad[2]=0.0;
			prop[m].solid.constantLoad[3]=0.0;
			prop[m].solid.constantLoad[4]=0.0;
			prop[m].solid.constantLoad[5]=0.0;
			prop[m].solid.constantLoadPos[0]+=minGlobal[0];
			prop[m].solid.constantLoadPos[1]+=minGlobal[1];
			prop[m].solid.constantLoadPos[2]+=minGlobal[2];
			
			/// Center of gravity initial setting 2/2
			for(int k=0; k<sim->dim; k++)
			{
				// minLocal[k]: the Minimum (X, Y, Z) solid particle position at Global Coordinate System
				// Computed at initial of simulation (t = 0)
				// prop[m].solid.cg[k]: Center of gravity in the Local Coordinate System (User input)
				prop[m].solid.cg[k]+=minLocal[k];
				prop[m].solid.cg0[k]=prop[m].solid.cg[k];
			}
			
			/// Transform CG if the rigid body is rotated at initial time t = 0s. Local->Global Coordinate System
			real cgGlobal[3], cgGlobal0[3];
			if (sim->dim == 2)
			{
				/********** USING ANGLE ZZ FOR 2D ************/
				///	cos()	-sin()
				///	sin()	 cos()
				cgGlobal[0]= prop[m].solid.cg[0]*cos(prop[m].solid.orient[2])-prop[m].solid.cg[1]*sin(prop[m].solid.orient[2]);
				cgGlobal[1]= prop[m].solid.cg[0]*sin(prop[m].solid.orient[2])+prop[m].solid.cg[1]*cos(prop[m].solid.orient[2]);
				prop[m].solid.cg[0]=cgGlobal[0];
				prop[m].solid.cg[1]=cgGlobal[1];
				cgGlobal0[0]= prop[m].solid.cg0[0]*cos(prop[m].solid.orient[2])-prop[m].solid.cg0[1]*sin(prop[m].solid.orient[2]);
				cgGlobal0[1]= prop[m].solid.cg0[0]*sin(prop[m].solid.orient[2])+prop[m].solid.cg0[1]*cos(prop[m].solid.orient[2]);
				prop[m].solid.cg0[0]=cgGlobal0[0];
				prop[m].solid.cg0[1]=cgGlobal0[1];
			}
			else if (sim->dim == 3)
			{
				/********** USING QUATERNNION FOR 3D ************/
				cgGlobal[0]=prop[m].solid.cg[0];
				cgGlobal[1]=prop[m].solid.cg[1];
				cgGlobal[2]=prop[m].solid.cg[2];
				tranformRotation3dLocal2Global(cgGlobal, prop[m].solid.quat);
				prop[m].solid.cg[0]=cgGlobal[0];
				prop[m].solid.cg[1]=cgGlobal[1];
				prop[m].solid.cg[2]=cgGlobal[2];
				cgGlobal0[0]=prop[m].solid.cg0[0];
				cgGlobal0[1]=prop[m].solid.cg0[1];
				cgGlobal0[2]=prop[m].solid.cg0[2];
				tranformRotation3dLocal2Global(cgGlobal0, prop[m].solid.quat);
				prop[m].solid.cg0[0]=cgGlobal0[0];
				prop[m].solid.cg0[1]=cgGlobal0[1];
				prop[m].solid.cg0[2]=cgGlobal0[2];
			}
			
			/// Print initial solid properties
			fprintf(stderr, "\n  Floating Body:\n");
			fprintf(stderr, "\tTotal Mass: %.6f", prop[m].solid.m);
			if(sim->dim==2)
			{
				fprintf(stderr, "\tMoment of Inertia Izz: %e\n", prop[m].solid.I[2][2]);
				fprintf(stderr, "\tReference Coordinates (Global coordinate system) - (MinX, MinY): (%.6f, %.6f) (MaxX, MaxY): (%.6f, %.6f)\n", minGlobal[0], minGlobal[1], maxGlobal[0], maxGlobal[1]);
				fprintf(stderr, "\tCoordinates of weight center (Global coordinate system) - (CGx, CGy): (%.6f, %.6f)\n", prop[m].solid.cg[0], prop[m].solid.cg[1]);
			}
			else
			{
				fprintf(stderr, "\n\t                                             [ %e %e %e ]\n", prop[m].solid.I[0][0], prop[m].solid.I[0][1], prop[m].solid.I[0][2]);
				fprintf(stderr, "\tMoment of Inertia (Local coordinate system): [ %e %e %e ]\n", prop[m].solid.I[1][0], prop[m].solid.I[1][1], prop[m].solid.I[1][2]);
				fprintf(stderr, "\t                                             [ %e %e %e ]\n", prop[m].solid.I[2][0], prop[m].solid.I[2][1], prop[m].solid.I[2][2]);
				fprintf(stderr, "\tReference Coordinates (Global coordinate system) - (MinX, MinY, MinZ): (%.6f, %.6f, %.6f) (MaxX, MaxY, MaxZ): (%.6f, %.6f, %.6f)\n", minGlobal[0], minGlobal[1], minGlobal[2], maxGlobal[0], maxGlobal[1], maxGlobal[2]);
				fprintf(stderr, "\tCoordinates of weight center (Global coordinate system) - (CGx, CGy, CGz): (%.6f, %.6f, %.6f)\n", prop[m].solid.cg[0], prop[m].solid.cg[1], prop[m].solid.cg[2]);
			}
			fprintf(stderr, "\tYoung: %e Poisson: %.2f Damping: %.6f Friction: %.6f \n", prop[m].solid.youngMod, prop[m].solid.poisson, prop[m].solid.damping, prop[m].solid.friction);


			/// Compute distance between wall particle and solid CG, Local Coordinate System
			/// Since the the body is rigid, these values are cte
			/// Loop for wall particles of solid material "m"
			for(int i=0; i<prop[m].part.n; i++)
			{
				int iPart=prop[m].part.v[i];
				
				/// dCG, Global Coordinate System
				prop[m].solid.dCG.x[i]=part->r.x[iPart]-prop[m].solid.cg[0];
				prop[m].solid.dCG.y[i]=part->r.y[iPart]-prop[m].solid.cg[1];
				if(sim->dim==3)
					prop[m].solid.dCG.z[i]=part->r.z[iPart]-prop[m].solid.cg[2];
                
				if (prop[m].solid.orient[0]!=0 || prop[m].solid.orient[1]!=0 || prop[m].solid.orient[2]!=0)
				{
					real dCGLocal[3];
					
					/// Transform dCG, Global->Local Coordinate System
					if (sim->dim == 2)
					{
						/********** USING ANGLE ZZ FOR 2D ************/
						///	cos()	 sin()
						///	-sin()	 cos()
						dCGLocal[0]= prop[m].solid.dCG.x[i]*cos(prop[m].solid.orient[2])+prop[m].solid.dCG.y[i]*sin(prop[m].solid.orient[2]);
						dCGLocal[1]=-prop[m].solid.dCG.x[i]*sin(prop[m].solid.orient[2])+prop[m].solid.dCG.y[i]*cos(prop[m].solid.orient[2]);
						prop[m].solid.dCG.x[i]=dCGLocal[0];
						prop[m].solid.dCG.y[i]=dCGLocal[1];
					}
					else if (sim->dim == 3)
					{
						/********** USING QUATERNNION FOR 3D ************/
						dCGLocal[0]=prop[m].solid.dCG.x[i];
						dCGLocal[1]=prop[m].solid.dCG.y[i];
						dCGLocal[2]=prop[m].solid.dCG.z[i];
						tranformRotation3dGlobal2Local(dCGLocal, prop[m].solid.quat);
						prop[m].solid.dCG.x[i]=dCGLocal[0];
						prop[m].solid.dCG.y[i]=dCGLocal[1];
						prop[m].solid.dCG.z[i]=dCGLocal[2];
					}
				}                     
			}

			/// Compute distance between dummy particle and solid CG, Local Coordinate System
			/// Loop for dummy particles of solid material "m"
			for(int i=0; i<prop[m].partDummy.n; i++)
			{
				int iPart=prop[m].partDummy.v[i];
				
				/// dCG, Global Coordinate System
				prop[m].solid.dCGDummy.x[i]=part->r.x[iPart]-prop[m].solid.cg[0];
				prop[m].solid.dCGDummy.y[i]=part->r.y[iPart]-prop[m].solid.cg[1];
				if(sim->dim==3)
					prop[m].solid.dCGDummy.z[i]=part->r.z[iPart]-prop[m].solid.cg[2];
                
				if (prop[m].solid.orient[0]!=0 || prop[m].solid.orient[1]!=0 || prop[m].solid.orient[2]!=0)
				{
					real dCGDummyLocal[3];
					
					/// Transform dCGDummy, Global->Local Coordinate System
					if (sim->dim == 2)
					{
						/********** USING ANGLE ZZ FOR 2D ************/
						///	cos()	 sin()
						///	-sin()	 cos()
						dCGDummyLocal[0]= prop[m].solid.dCGDummy.x[i]*cos(prop[m].solid.orient[2])+prop[m].solid.dCGDummy.y[i]*sin(prop[m].solid.orient[2]);
						dCGDummyLocal[1]=-prop[m].solid.dCGDummy.x[i]*sin(prop[m].solid.orient[2])+prop[m].solid.dCGDummy.y[i]*cos(prop[m].solid.orient[2]);
						prop[m].solid.dCGDummy.x[i]=dCGDummyLocal[0];
						prop[m].solid.dCGDummy.y[i]=dCGDummyLocal[1];
					}
					else if (sim->dim == 3)
					{
						/********** USING QUATERNNION FOR 3D ************/
						dCGDummyLocal[0]=prop[m].solid.dCGDummy.x[i];
						dCGDummyLocal[1]=prop[m].solid.dCGDummy.y[i];
						dCGDummyLocal[2]=prop[m].solid.dCGDummy.z[i];
						tranformRotation3dGlobal2Local(dCGDummyLocal, prop[m].solid.quat);
						prop[m].solid.dCGDummy.x[i]=dCGDummyLocal[0];
						prop[m].solid.dCGDummy.y[i]=dCGDummyLocal[1];
						prop[m].solid.dCGDummy.z[i]=dCGDummyLocal[2];
					}
				}
			}

			/// Compute Principal Axes and Moments of Inertia (Only for 3D)
			if(sim->dim==3)
			{
				if((prop[m].solid.I[0][1]!=prop[m].solid.I[1][0]) || 
                                   (prop[m].solid.I[0][2]!=prop[m].solid.I[2][0]) ||
                                   (prop[m].solid.I[1][2]!=prop[m].solid.I[2][1]))
				{
					fprintf(stderr, "\nError: The Inertia Matrix is non-symmetric.\n\n");
					exit(0);
				}
				if((prop[m].solid.I[0][1]!=0) || (prop[m].solid.I[0][2]!=0) || (prop[m].solid.I[1][2]!=0))
				{
					real B[3][3];       // Principal Axes
					real Ip[3];         // Principal Moments
					real base[3];

					/// Eigenvalue and Eigenvector
					eigen_decomposition(prop[m].solid.I, B, Ip);

					/// Principal Moments of Inertia
					prop[m].solid.I[0][0]=Ip[0];
					prop[m].solid.I[1][1]=Ip[1];
					prop[m].solid.I[2][2]=Ip[2];

					/// Products of Inertia
					prop[m].solid.I[0][1]=prop[m].solid.I[1][0]=0.0;
					prop[m].solid.I[0][2]=prop[m].solid.I[2][0]=0.0;
					prop[m].solid.I[1][2]=prop[m].solid.I[2][1]=0.0;

					/// If the principal axes are acceptable
					base[0]=B[1][0]*B[2][1]-B[1][1]*B[2][0];
					base[1]=B[0][1]*B[2][0]-B[2][1]*B[0][0];
					base[2]=B[0][0]*B[1][1]-B[1][0]*B[0][1];
					if((base[0]!=B[0][2]) || (base[1]!=B[1][2]) || (base[2]!=B[2][2]))
					{
						B[0][0]=-B[0][0];
						B[1][0]=-B[1][0];
						B[2][0]=-B[2][0];
					}

					real trace=B[0][0]+B[1][1]+B[2][2];

					/// Matrix B to Quaternion
					if(trace>0)         // I changed M_EPSILON to 0
					{
						real s=0.5/sqrt(trace+1.0);

						prop[m].solid.quat[0]=0.25/s;
						prop[m].solid.quat[1]=(B[2][1]-B[1][2])*s;
						prop[m].solid.quat[2]=(B[0][2]-B[2][0])*s;
						prop[m].solid.quat[3]=(B[1][0]-B[0][1])*s;
					}
					else
					{
						if((B[0][0]>B[1][1]) && (B[0][0]>B[2][2]))
						{
							real s=2.0*sqrtf(1.0+B[0][0]-B[1][1]-B[2][2]);

							prop[m].solid.quat[0]=(B[2][1]-B[1][2])/s;
							prop[m].solid.quat[1]=0.25*s;
							prop[m].solid.quat[2]=(B[0][1]+B[1][0])/s;
							prop[m].solid.quat[3]=(B[0][2]+B[2][0])/s;
						}
						else if(B[1][1]>B[2][2])
						{
							real s=2.0*sqrtf(1.0+B[1][1]-B[0][0]-B[2][2]);

							prop[m].solid.quat[0]=(B[0][2]-B[2][0])/s;
							prop[m].solid.quat[1]=(B[0][1]+B[1][0])/s;
							prop[m].solid.quat[2]=0.25*s;
							prop[m].solid.quat[3]=(B[1][2]+B[2][1])/s;
						}
						else
						{
							real s=2.0*sqrtf(1.0+B[2][2]-B[0][0]-B[1][1]);

							prop[m].solid.quat[0]=(B[1][0]-B[0][1])/s;
							prop[m].solid.quat[1]=(B[0][2]+B[2][0])/s;
							prop[m].solid.quat[2]=(B[1][2]+B[2][1])/s;
							prop[m].solid.quat[3]=0.25*s;
						}
					}
					fprintf(stderr, "\n\t   Principal       [ %.3f %.3f %.3f ]\n", B[0][0], B[0][1], B[0][2]);
					fprintf(stderr, "\tAxes of Inertia =[ %.3f %.3f %.3f ]\n", B[1][0], B[1][1], B[1][2]);
					fprintf(stderr, "\t                   [ %.3f %.3f %.3f ]\n", B[2][0], B[2][1], B[2][2]);
					fprintf(stderr, "\tPrincipal Moments of Inertia=[ %.3f %.3f %.3f ]\n", prop[m].solid.I[0][0], prop[m].solid.I[1][1], prop[m].solid.I[2][2]);
					// fprintf(stderr,"\tQuaternion=[ %.3f %.3f %.3f %.3f]\n", prop[m].solid.quat[0], prop[m].solid.quat[1], prop[m].solid.quat[2], prop[m].solid.quat[3]);
					
					/// Loop for wall particles of solid material "m"
					for(int i=0; i<prop[m].part.n; i++)
					{
						int iPart=prop[m].part.v[i];
						real nLocal[3];
						
						/// Transform Normal, Global->Local Coordinate System
						nLocal[0]=prop[m].solid.normalA.x[i];
						nLocal[1]=prop[m].solid.normalA.y[i];
						nLocal[2]=prop[m].solid.normalA.z[i];
						tranformRotation3dGlobal2Local(nLocal, prop[m].solid.quat);
						prop[m].solid.normalA.x[i]=nLocal[0];
						prop[m].solid.normalA.y[i]=nLocal[1];
						prop[m].solid.normalA.z[i]=nLocal[2];
						part->normal.x[iPart]=prop[m].solid.normalA.x[i];
						part->normal.y[iPart]=prop[m].solid.normalA.y[i];
						part->normal.z[iPart]=prop[m].solid.normalA.z[i];

						real dCGLocal[3];
						/// Transform dCG, Global->Local Coordinate System
						dCGLocal[0]=prop[m].solid.dCG.x[i];
						dCGLocal[1]=prop[m].solid.dCG.y[i];
						dCGLocal[2]=prop[m].solid.dCG.z[i];
						tranformRotation3dGlobal2Local(dCGLocal, prop[m].solid.quat);
						prop[m].solid.dCG.x[i]=dCGLocal[0];
						prop[m].solid.dCG.y[i]=dCGLocal[1];
						prop[m].solid.dCG.z[i]=dCGLocal[2];
					}

					/// Loop for dummy particles of solid material "m"
					for(int i=0; i<prop[m].partDummy.n; i++)
					{
						real dCGDummyLocal[3];
						/// Transform dCGDummy, Global->Local Coordinate System
						dCGDummyLocal[0]=prop[m].solid.dCGDummy.x[i];
						dCGDummyLocal[1]=prop[m].solid.dCGDummy.y[i];
						dCGDummyLocal[2]=prop[m].solid.dCGDummy.z[i];
						tranformRotation3dGlobal2Local(dCGDummyLocal, prop[m].solid.quat);
						prop[m].solid.dCGDummy.x[i]=dCGDummyLocal[0];
						prop[m].solid.dCGDummy.y[i]=dCGDummyLocal[1];
						prop[m].solid.dCGDummy.z[i]=dCGDummyLocal[2];
					}
				}
			}
		}
	}

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Verify if the simulation has collision/contact between rigid bodies based on user input
int initSolidCollision(particles *part, simParameters *sim, matProperties *prop)
{
	/// Set rigid body surface wall features - Edge, Vertex and Face
	calcSolidBC(part, sim, prop);
	
	/// Check that the collision between rigid bodies can occur
	/// Loop for all nMat materials, e.g.:
	/// fluid-0 -> m = 0, fluid-1 -> m = 1, ..., fixed-solid-1 -> m = 5, ..., free-solid-1 -> m = 7)
	for(int m=0; m<sim->nMat; m++)
	{
		/// Verify if "m" is a free solid and have particles, i.e., if the rigid body exists
		if((prop[m].type==solid) && (prop[m].motion==freeMot) && (prop[m].part.n>0))
		{
			/// Verify if the rigid body has at least one collision/contact property (User input)
			if((prop[m].solid.youngMod!=0) || (prop[m].solid.damping!=0) || (prop[m].solid.friction!=0))
			{
				/// Compute the minium rigid body time step (dt_solid) as a multiple of the fluid time step (dt_fluid)
				/// If dt_solid < dt_fluid, rigid body motion is updated in a subcyling temporal integration
				solidTimeStep(sim, part, prop, m);
				
				collision_ok=true;
				break;
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
/// Rigid Body Motion - Main loop
int calcSolidMotion(particles *part, simParameters *sim, matProperties *prop)
{
	/// Next solid time step (auxiliar variables)
	double next_dtSolidRigid = sim->dt;	// If no collision, this value doesn't change
	int next_nStepsSolid = 1;			// If no collision, this value doesn't change
	
	/// Subcycling temporal integration of rigid bodies motion
	/// Loop for steps in the subcycling
	for(int ii=0; ii<sim->nStepsSolid; ii++)
	{
		/// Loop for all nMat materials
		for(int m=0; m<sim->nMat; m++)
		{
			/// Verify if "m" is a solid (fixed or free) and have particles, i.e., if exists
			if((prop[m].type==solid) && (prop[m].part.n>0))
			{
				/// Set collision forces and moments to Zero
				for(int k=0; k<3; k++)
				{
					prop[m].solid.collisionForce[k]=0.0;
					prop[m].solid.collisionMoment[k]=0.0;
					prop[m].solid.globalCollisionForce[k] = 0.0;
					prop[m].solid.globalCollisionMoment[k] = 0.0;
					//prop[m].solid.globalHydroForce[k] = 0.0;
				}
			}
		}

		/// Colision Detection. If simulation has collision/contact
		if(collision_ok)
		{
			/// Loop for all nMat materials
			for(int m=0; m<sim->nMat; m++)
			{
				/// Verify if "m" is a free solid and have particles, i.e., if the solid exists
				if((prop[m].type==solid) && (prop[m].motion==freeMot) && (prop[m].part.n>0))
				{
					/// Compute collision/contact force and moment (torque) at CG of rigid body "m"
					collisionForceMoment(sim, part, prop, m, next_dtSolidRigid, next_nStepsSolid);
				}
			}
		}

		/// Loop for all nMat materials
		for(int m=0; m<sim->nMat; m++)
		{
			/// Verify if "m" is a free solid and have particles, i.e., if the solid exists
			if((prop[m].type==solid) && (prop[m].motion==freeMot) && (prop[m].part.n>0))
			{
				real force[3], moment[3];
				/// Set collision/contact force and moment at CG of rigid body "m". (If no collision, the value is zero)
				for(int k=0; k<3; k++)
				{
					force[k]=prop[m].solid.collisionForce[k];
					moment[k]=prop[m].solid.collisionMoment[k];
				}

				/// Add hydrodynamic force and moment (torque) at CG of the Rigid Body "m"
				calcPressureIntegration(part, prop, sim->dim, m, force, moment);
				/// Add gravitional force at CG of the Rigid Body "m"
				addGravityForce(sim, prop, m, force);
				
				/// Rigid Body equations of motion
				/// Update velocity (translational and angular), position and rotation at CG of the Rigid Body "m"
				/// Update positions and velocities of Wall and Dummy Particles belonging to the Rigid Body "m"
				if(sim->dim==2)
				{
					calc2dMotionEquation(part, sim, prop, m, force, moment);
				}
				else
				{
					calc3dMotionEquation(part, sim, prop, m, force, moment);
				}
			}
		}
	}

	/// Velocity/Position correction
	/// Velocity and position are updated here in the subcycling and will be updated
	/// again in updateExplicitMotion function of MPS fluid. Therefore, this "weird" correction is necessary.
	/// Loop for all nMat materials
	for(int m=0; m<sim->nMat; m++)
	{
		/// Verify if "m" is a free solid and have particles, i.e., if the solid exists
		if((prop[m].type==solid) && (prop[m].motion==freeMot) && (prop[m].part.n>0))
		{
			/// Loop for wall particles of solid material "m"
#pragma omp parallel for
			for(int i=0; i<prop[m].part.n; i++)
			{
				//	i:  Particle ID considering only the solid material m
				//	iPart: Particle ID considering all materials, i.e., fluid, fixed solid, free solid, ...
				int iPart=prop[m].part.v[i];

				// Position and velocity correction
				part->rn.x[iPart]-=part->un.x[iPart]*sim->dt;
				part->un.x[iPart]-=part->du.x[iPart];
				part->rn.y[iPart]-=part->un.y[iPart]*sim->dt;
				part->un.y[iPart]-=part->du.y[iPart];
				if(sim->dim==3)
				{
					part->rn.z[iPart]-=part->un.z[iPart]*sim->dt;
					part->un.z[iPart]-=part->du.z[iPart];
				}
			}
			/// Loop for dummy particles of solid material "m"
			for(int i=0; i<prop[m].partDummy.n; i++)
			{
				//	i:  Particle ID considering only the solid material m
				//	iPart: Particle ID considering all materials, i.e., fluid, fixed solid, free solid, ...
				int iPart=prop[m].partDummy.v[i];

				// Position and velocity correction
				part->rn.x[iPart]-=part->un.x[iPart]*sim->dt;
				part->un.x[iPart]-=part->du.x[iPart];
				part->rn.y[iPart]-=part->un.y[iPart]*sim->dt;
				part->un.y[iPart]-=part->du.y[iPart];
				if(sim->dim==3)
				{
					part->rn.z[iPart]-=part->un.z[iPart]*sim->dt;
					part->un.z[iPart]-=part->du.z[iPart];
				}
			}
		}
	}
	
	/// Values used in the next step (may change if collision occurs!)
	sim->dtSolidRigid = next_dtSolidRigid;
	sim->nStepsSolid = next_nStepsSolid;
	
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////
/// Compute hydrodynamic force and moment (torque) on the surface of the Rigid Body "m", Local Coordinate System
/// Loads due to fluid pressure on the solid wall
int calcPressureIntegration(particles *part, matProperties *prop, const int dim, const int m, real *force, real *moment)
{
	/// Set forces and moments to zero
	for(int k=0; k<3; k++)
	{
		prop[m].solid.hydroForce[k]=0.0;
		prop[m].solid.hydroMoment[k]=0.0;
		prop[m].solid.globalHydroForce[k] = 0.0;
        prop[m].solid.globalHydroMoment[k] = 0.0;
	}
	
	if(dim==2)
	{
		/// Loop for wall particles of solid material "m"
		for(int i=0; i<prop[m].part.n; i++)
		{
			int iPart=prop[m].part.v[i];
			real forceAux[3];
			
			/// Force = sum(Press_i * Area_i * Normal_i), Local Coordinate System
			forceAux[0]=-part->p.v[iPart]*prop[m].solid.areaA.v[i]*prop[m].solid.normalA.x[i];
			forceAux[1]=-part->p.v[iPart]*prop[m].solid.areaA.v[i]*prop[m].solid.normalA.y[i];
			prop[m].solid.hydroForce[0]+=forceAux[0];
			prop[m].solid.hydroForce[1]+=forceAux[1];

			/// Moment = sum[dCG X (Press_i * Area_i * Normal_i)], Local Coordinate System
			prop[m].solid.hydroMoment[2]+=prop[m].solid.dCG.x[i]*forceAux[1]-prop[m].solid.dCG.y[i]*forceAux[0];
		}
	}
	else
	{
		/// Loop for wall particles of solid material "m"
		for(int i=0; i<prop[m].part.n; i++)
		{
			int iPart=prop[m].part.v[i];
			real forceAux[3];
			
			/// Force = sum(Press_i * Area_i * Normal_i), Local Coordinate System
			forceAux[0]=-part->p.v[iPart]*prop[m].solid.areaA.v[i]*prop[m].solid.normalA.x[i];
			forceAux[1]=-part->p.v[iPart]*prop[m].solid.areaA.v[i]*prop[m].solid.normalA.y[i];
			forceAux[2]=-part->p.v[iPart]*prop[m].solid.areaA.v[i]*prop[m].solid.normalA.z[i];
			prop[m].solid.hydroForce[0]+=forceAux[0];
			prop[m].solid.hydroForce[1]+=forceAux[1];
			prop[m].solid.hydroForce[2]+=forceAux[2];
                        
			/// Moment = sum[dCG X (Press_i * Area_i * Normal_i)], Local Coordinate System
			prop[m].solid.hydroMoment[0]+=prop[m].solid.dCG.y[i]*forceAux[2]-prop[m].solid.dCG.z[i]*forceAux[1];
			prop[m].solid.hydroMoment[1]+=prop[m].solid.dCG.z[i]*forceAux[0]-prop[m].solid.dCG.x[i]*forceAux[2];
			prop[m].solid.hydroMoment[2]+=prop[m].solid.dCG.x[i]*forceAux[1]-prop[m].solid.dCG.y[i]*forceAux[0];
		}
	}
	
	/// Update force and moment at CG of the Rigid Body "m", Local Coordinate System
	for(int k=0; k<3; k++)
	{
		force[k]+=prop[m].solid.hydroForce[k];
		moment[k]+=prop[m].solid.hydroMoment[k];
	}
	

	/// Global forces and moments saved in file _motion.txt
	/// Local->Global Coordinate System
    if (dim==2)
    {
        /********** USING ANGLE ZZ FOR 2D ************/
        ///	cos()	-sin()
		///	sin()	 cos()
        prop[m].solid.globalHydroForce[0] = prop[m].solid.hydroForce[0] * cos(prop[m].solid.orient[2])
                        - prop[m].solid.hydroForce[1] * sin(prop[m].solid.orient[2]);
        prop[m].solid.globalHydroForce[1] = prop[m].solid.hydroForce[0] * sin(prop[m].solid.orient[2])
                        + prop[m].solid.hydroForce[1] * cos(prop[m].solid.orient[2]);
		prop[m].solid.globalHydroMoment[2] = prop[m].solid.hydroMoment[2];
    }
    else
    {
		/********** USING QUATERNNION FOR 3D ************/
		for(int k=0; k<3; k++)
		{
			prop[m].solid.globalHydroForce[k] = prop[m].solid.hydroForce[k];
			prop[m].solid.globalHydroMoment[k] = prop[m].solid.hydroMoment[k];
		}
        tranformRotation3dLocal2Global(prop[m].solid.globalHydroForce, prop[m].solid.quat);
        tranformRotation3dLocal2Global(prop[m].solid.globalHydroMoment, prop[m].solid.quat);
    }

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Add gravity to solid force of the Rigid Body "m", Local Coordinate System
int addGravityForce(simParameters *sim, matProperties *prop, const int m, real *force)
{
	if(sim->dim==2)
	{
		/********** USING ANGLE ZZ FOR 2D ************/
		///	cos()	 sin()
		///	-sin()	 cos()
		real solidWeight[3];

		/// Transform gravitational force, Global->Local Coordinate System
		solidWeight[0]=sim->g[0]*prop[m].solid.m*cos(prop[m].solid.orient[2])+sim->g[1]*prop[m].solid.m*sin(prop[m].solid.orient[2]);
		solidWeight[1]=-sim->g[0]*prop[m].solid.m*sin(prop[m].solid.orient[2])+sim->g[1]*prop[m].solid.m*cos(prop[m].solid.orient[2]);
		/// Add to force at CG of the Rigid Body "m", Local Coordinate System
		force[0]+=solidWeight[0];
		force[1]+=solidWeight[1];
	}
	else
	{
        /********** USING QUATERNNION FOR 3D ************/
		real solidWeight[3];

		/// Transform gravitational force, Global->Local Coordinate System
		solidWeight[0]=sim->g[0]*prop[m].solid.m;
		solidWeight[1]=sim->g[1]*prop[m].solid.m;
		solidWeight[2]=sim->g[2]*prop[m].solid.m;
		tranformRotation3dGlobal2Local(solidWeight, prop[m].solid.quat);
		/// Add to force at CG of the Rigid Body "m", Local Coordinate System
		force[0]+=solidWeight[0];		
		force[1]+=solidWeight[1];
		force[2]+=solidWeight[2];
	}

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Compute 2D motion of the Rigid Body "m", Global Coordinate System
int calc2dMotionEquation(particles *part, simParameters *sim, matProperties *prop, const int m, real *force, real *moment)
{
	real trans[3];
	real transAux[3];
	real rotate[3];
	real wcAux[3];

	/// Transform translational velocity at CG of the Rigid Body "m", Global->Local Coordinate System
	/********** USING ANGLE ZZ FOR 2D ************/
	///	cos()	 sin()
	///	-sin()	 cos()
	wcAux[0]=prop[m].solid.vel[0]*cos(prop[m].solid.orient[2])+prop[m].solid.vel[1]*sin(prop[m].solid.orient[2]);
	wcAux[1]=-prop[m].solid.vel[0]*sin(prop[m].solid.orient[2])+prop[m].solid.vel[1]*cos(prop[m].solid.orient[2]);
	prop[m].solid.vel[0]=wcAux[0];
	prop[m].solid.vel[1]=wcAux[1];

	/// Time integration method (Euler=0, Runge-Kutta=1)
	int integrationMethod=0;

	/// Update translational velocity and position at CG of the Rigid Body "m", Local Coordinate System
	/// Velocity[t+dt] = Force/Mass*dt + Velocity[t]
	/// Euler integration
	if(integrationMethod==0)
	{
		for(int k=0; k<sim->dim; k++)
		{
			prop[m].solid.vel[k]+=force[k]/prop[m].solid.m*sim->dtSolidRigid;
			trans[k]=prop[m].solid.vel[k]*sim->dtSolidRigid;
		}
	}
	/// Runge-kutta integration
	else if(integrationMethod==1)
	{
		for(int k=0; k<sim->dim; k++)
		{
			real rk[4];

			rk[0]=prop[m].solid.f1[k]/prop[m].solid.m;
			rk[1]=prop[m].solid.f2[k]/prop[m].solid.m;
			rk[2]=prop[m].solid.f2[k]/prop[m].solid.m;
			rk[3]=force[k]/prop[m].solid.m;
			
			prop[m].solid.vel[k]=prop[m].solid.vel1[k]+2*sim->dtSolidRigid*(rk[0]+2*rk[1]+2*rk[2]+rk[3])/6;
			trans[k]=prop[m].solid.vel[k]*sim->dtSolidRigid;

			prop[m].solid.f1[k]=prop[m].solid.f2[k];
			prop[m].solid.vel1[k]=prop[m].solid.vel2[k];
			prop[m].solid.f2[k]=force[k];
			prop[m].solid.vel2[k]=prop[m].solid.vel[k];
		}
	}

	/// Update angular velocity and rotation at CG of the Rigid Body "m", Local Coordinate System
	/// Angular_Velocity[t+dt] = Moment/Inertia*dt + Angular_Velocity[t]
	/// Euler integration
	if(integrationMethod==0)
	{
		prop[m].solid.velAng[2]+=moment[2]/prop[m].solid.I[2][2]*sim->dtSolidRigid;
		rotate[2]=prop[m].solid.velAng[2]*sim->dtSolidRigid;
	}
	/// Runge-kutta integration
	else if(integrationMethod==1)
	{
		double rk[4];

		rk[0]=prop[m].solid.m1[2]/prop[m].solid.I[2][2];
		rk[1]=prop[m].solid.m2[2]/prop[m].solid.I[2][2];
		rk[2]=prop[m].solid.m2[2]/prop[m].solid.I[2][2];
		rk[3]=moment[2]/prop[m].solid.I[2][2];

		prop[m].solid.velAng[2]=prop[m].solid.velAng1[2]+2*sim->dtSolidRigid*(rk[0]+2*rk[1]+2*rk[2]+rk[3])/6;
		rotate[2]=prop[m].solid.velAng[2]*sim->dtSolidRigid;

		prop[m].solid.m1[2]=prop[m].solid.m2[2];
		prop[m].solid.velAng1[2]=prop[m].solid.velAng2[2];
		prop[m].solid.m2[2]=moment[2];
		prop[m].solid.velAng2[2]=prop[m].solid.velAng[2];
	}

	/// Motion directions restrictions, Local Coordinate System (User input)
	trans[0]*=prop[m].solid.freeDir[0];
	trans[1]*=prop[m].solid.freeDir[1];
	rotate[2]*=prop[m].solid.freeDir[5];

	/// Transform CG motion, Local->Global Coordinate System
	/********** USING ANGLE ZZ FOR 2D ************/
	///	cos()	-sin()
	///	sin()	 cos()
	transAux[0]=trans[0]*cos(prop[m].solid.orient[2])-trans[1]*sin(prop[m].solid.orient[2]);
	transAux[1]=trans[0]*sin(prop[m].solid.orient[2])+trans[1]*cos(prop[m].solid.orient[2]);
	trans[0]=transAux[0];
	trans[1]=transAux[1];

	/// External constant load, Global Coordinate System
	real pos[2], aux[2];
	aux[0]=prop[m].solid.constantLoadPos[0]-prop[m].solid.cg[0];
	aux[1]=prop[m].solid.constantLoadPos[1]-prop[m].solid.cg[1];
	pos[0]=aux[0]*cos(rotate[2])-aux[1]*sin(rotate[2]);
	pos[1]=aux[0]*sin(rotate[2])+aux[1]*cos(rotate[2]);
	prop[m].solid.constantLoadPos[0]=prop[m].solid.cg[0]+trans[0]+pos[0];
	prop[m].solid.constantLoadPos[1]=prop[m].solid.cg[1]+trans[1]+pos[1];

	/// Update Positions and Velocities of Wall and Dummy Particles belonging to the Rigid Body "m", Global Coordinate System
	for(int k=0; k<3; k++)
	{
		if(isnan(trans[k])!=0)
			trans[k]=0;
		if(isnan(rotate[k])!=0)
			rotate[k]=0;
	}
	update2dPosVel(part, sim, prop, trans, rotate, m);

	/// Update CG Position and orientation of the Rigid Body "m", Global Coordinate System
	prop[m].solid.cg[0]+=trans[0];
	prop[m].solid.cg[1]+=trans[1];
	prop[m].solid.orient[2]+=rotate[2];

	/// Update CG Velocity of the Rigid Body "m", Global Coordinate System
	prop[m].solid.vel[0]=trans[0]/sim->dtSolidRigid;
	prop[m].solid.vel[1]=trans[1]/sim->dtSolidRigid;
	prop[m].solid.velAng[2]=rotate[2]/sim->dtSolidRigid;

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Compute 3D motion of the of the Rigid Body "m", Global Coordinate System
int calc3dMotionEquation(particles *part, simParameters *sim, matProperties *prop, const int m, real *force, real *moment)
{
	real trans[3];
	real rotate[3];

	/// Transform translational velocity at CG of the Rigid Body "m", Global->Local Coordinate System
    /********** USING QUATERNNION FOR 3D ************/
	tranformRotation3dGlobal2Local(prop[m].solid.vel, prop[m].solid.quat);

	/// Set rigid body mass (homogeneous mass assumed)
	prop[m].solid.M[0][0]=prop[m].solid.m;
	prop[m].solid.M[0][1]=0.0;
	prop[m].solid.M[0][2]=0.0;
	prop[m].solid.M[1][0]=0.0;
	prop[m].solid.M[1][1]=prop[m].solid.m;
	prop[m].solid.M[1][2]=0.0;
	prop[m].solid.M[2][0]=0.0;
	prop[m].solid.M[2][1]=0.0;
	prop[m].solid.M[2][2]=prop[m].solid.m;

	/// Update translational velocity and position at CG of the Rigid Body "m", Local Coordinate System
	/// Velocity[t+dt] = Force/Mass*dt + Velocity[t]
	/// Euler integration
	prop[m].solid.vel[0]+=sim->dtSolidRigid*force[0]/prop[m].solid.M[0][0];
	prop[m].solid.vel[1]+=sim->dtSolidRigid*force[1]/prop[m].solid.M[1][1];
	prop[m].solid.vel[2]+=sim->dtSolidRigid*force[2]/prop[m].solid.M[2][2];
	trans[0]=prop[m].solid.vel[0]*sim->dtSolidRigid;
	trans[1]=prop[m].solid.vel[1]*sim->dtSolidRigid;
	trans[2]=prop[m].solid.vel[2]*sim->dtSolidRigid;

	/// Transform angular velocity at CG of the Rigid Body "m", Global->Local Coordinate System
	/********** USING QUATERNNION FOR 3D ************/
	tranformRotation3dGlobal2Local(prop[m].solid.velAng, prop[m].solid.quat);
	
	/// Update angular velocity and rotation at CG of the Rigid Body "m", Local Coordinate System
	/// Angular_Velocity[t+dt] = Moment/Inertia*dt + Angular_Velocity[t]
	/// Euler integration
	real velAng[3];

	velAng[0]=prop[m].solid.velAng[0];
	velAng[1]=prop[m].solid.velAng[1];
	velAng[2]=prop[m].solid.velAng[2];

	prop[m].solid.velAng[0]+=sim->dtSolidRigid*(moment[0]-velAng[1]*velAng[2]*(prop[m].solid.I[2][2]-prop[m].solid.I[1][1]))/prop[m].solid.I[0][0];
	prop[m].solid.velAng[1]+=sim->dtSolidRigid*(moment[1]-velAng[2]*velAng[0]*(prop[m].solid.I[0][0]-prop[m].solid.I[2][2]))/prop[m].solid.I[1][1];
	prop[m].solid.velAng[2]+=sim->dtSolidRigid*(moment[2]-velAng[0]*velAng[1]*(prop[m].solid.I[1][1]-prop[m].solid.I[0][0]))/prop[m].solid.I[2][2];

	rotate[0]=prop[m].solid.velAng[0]*sim->dtSolidRigid;
	rotate[1]=prop[m].solid.velAng[1]*sim->dtSolidRigid;
	rotate[2]=prop[m].solid.velAng[2]*sim->dtSolidRigid;

	/// Motion directions restrictions, Local Coordinate System (User input)
	trans[0]*=prop[m].solid.freeDir[0];
	trans[1]*=prop[m].solid.freeDir[1];
	trans[2]*=prop[m].solid.freeDir[2];
	rotate[0]*=prop[m].solid.freeDir[3];
	rotate[1]*=prop[m].solid.freeDir[4];
	rotate[2]*=prop[m].solid.freeDir[5];

	/// Transform CG translation and rotation, Local->Global Coordinate System
    /********** USING QUATERNNION FOR 3D ************/
	tranformRotation3dLocal2Global(trans, prop[m].solid.quat);
	tranformRotation3dLocal2Global(rotate, prop[m].solid.quat);

	/// Update Positions and Velocities of Wall and Dummy Particles belonging to the Rigid Body "m", Global Coordinate System
	update3dPosVel(part, sim, prop, trans, rotate, m);

	/// Update CG Position and orientation of the Rigid Body "m", Global Coordinate System
	prop[m].solid.cg[0]+=trans[0];
	prop[m].solid.cg[1]+=trans[1];
	prop[m].solid.cg[2]+=trans[2];
	prop[m].solid.orient[0]+=rotate[0];
	prop[m].solid.orient[1]+=rotate[1];
	prop[m].solid.orient[2]+=rotate[2];

	/// Update CG Velocity of the Rigid Body "m", Global Coordinate System
	prop[m].solid.vel[0]=trans[0]/sim->dtSolidRigid;
	prop[m].solid.vel[1]=trans[1]/sim->dtSolidRigid;
	prop[m].solid.vel[2]=trans[2]/sim->dtSolidRigid;
	prop[m].solid.velAng[0]=rotate[0]/sim->dtSolidRigid;
	prop[m].solid.velAng[1]=rotate[1]/sim->dtSolidRigid;
	prop[m].solid.velAng[2]=rotate[2]/sim->dtSolidRigid;

	// External constant load, Global Coordinate System
	real pos[3];
	pos[0]=prop[m].solid.constantLoadPos[0]-prop[m].solid.cg[0];
	pos[1]=prop[m].solid.constantLoadPos[1]-prop[m].solid.cg[1];
	pos[2]=prop[m].solid.constantLoadPos[2]-prop[m].solid.cg[2];
	tranformRotation3dLocal2Global(pos, prop[m].solid.quat);
	prop[m].solid.constantLoadPos[0]=prop[m].solid.cg[0]+pos[0];
	prop[m].solid.constantLoadPos[1]=prop[m].solid.cg[1]+pos[1];
	prop[m].solid.constantLoadPos[2]=prop[m].solid.cg[2]+pos[2];
	
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Update 2D Positions and Velocities of Wall and Dummy Particles belonging to the Rigid Body "m", Global Coordinate System
int update2dPosVel(particles *part, simParameters *sim, matProperties *prop, double trans[3], double rotate[3], const int m)
{
	/// Loop for wall particles of solid material "m"
	for(int i=0; i<prop[m].part.n; i++)
	{
		int iPart=prop[m].part.v[i];
		double drGlobal[3], drAux[3];
		
		/// Transform distance between wall particle and solid CG, Local->Global Coordinate System
		/********** USING ANGLE ZZ FOR 2D ************/
		///	cos()	-sin()
		///	sin()	 cos()
		drAux[0]=prop[m].solid.dCG.x[i]*cos(prop[m].solid.orient[2])-prop[m].solid.dCG.y[i]*sin(prop[m].solid.orient[2]);
		drAux[1]=prop[m].solid.dCG.x[i]*sin(prop[m].solid.orient[2])+prop[m].solid.dCG.y[i]*cos(prop[m].solid.orient[2]);
		
		/// Wall particle displacement, Global Coordinate System
		drGlobal[0]=trans[0]+drAux[0]*cos(rotate[2])-drAux[1]*sin(rotate[2]);
		drGlobal[1]=trans[1]+drAux[0]*sin(rotate[2])+drAux[1]*cos(rotate[2]);

		/// Update Particles Velocities and Positions, Global Coordinate System
		part->un.x[iPart]=(prop[m].solid.cg[0]+drGlobal[0]-part->rn.x[iPart])/sim->dtSolidRigid;
		part->un.y[iPart]=(prop[m].solid.cg[1]+drGlobal[1]-part->rn.y[iPart])/sim->dtSolidRigid;
		part->rn.x[iPart]=prop[m].solid.cg[0]+drGlobal[0];
		part->rn.y[iPart]=prop[m].solid.cg[1]+drGlobal[1];
	}
	/// Loop for dummy particles of solid material "m"
	for(int i=0; i<prop[m].partDummy.n; i++)
	{
		int iPart=prop[m].partDummy.v[i];
		double drGlobal[3], drAux[3];
		
		/// Transform distance between dummy particle and solid CG, Local->Global Coordinate System
		/********** USING ANGLE ZZ FOR 2D ************/
		///	cos()	-sin()
		///	sin()	 cos()
		drAux[0]=prop[m].solid.dCGDummy.x[i]*cos(prop[m].solid.orient[2])-prop[m].solid.dCGDummy.y[i]*sin(prop[m].solid.orient[2]);
		drAux[1]=prop[m].solid.dCGDummy.x[i]*sin(prop[m].solid.orient[2])+prop[m].solid.dCGDummy.y[i]*cos(prop[m].solid.orient[2]);
		
		/// Dummy particle displacement, Global Coordinate System
		drGlobal[0]=trans[0]+drAux[0]*cos(rotate[2])-drAux[1]*sin(rotate[2]);
		drGlobal[1]=trans[1]+drAux[0]*sin(rotate[2])+drAux[1]*cos(rotate[2]);

		/// Update Particles Velocities and Positions, Global Coordinate System
		part->un.x[iPart]=(prop[m].solid.cg[0]+drGlobal[0]-part->rn.x[iPart])/sim->dtSolidRigid;
		part->un.y[iPart]=(prop[m].solid.cg[1]+drGlobal[1]-part->rn.y[iPart])/sim->dtSolidRigid;
		part->rn.x[iPart]=prop[m].solid.cg[0]+drGlobal[0];
		part->rn.y[iPart]=prop[m].solid.cg[1]+drGlobal[1];
	}

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Update 3D Positions and Velocities of Wall and Dummy Particles belonging to the Rigid Body "m", Global Coordinate System
int update3dPosVel(particles *part, simParameters *sim, matProperties *prop, double trans[3], double rotate[3], const int m)
{
	/// Update Quaternion
	real rot[3], dTheta;

	/// 3D delta rotation
	dTheta=sqrt(rotate[0]*rotate[0]+rotate[1]*rotate[1]+rotate[2]*rotate[2]);
	if(dTheta!=0)
	{
		rot[0]=rotate[0]/dTheta;
		rot[1]=rotate[1]/dTheta;
		rot[2]=rotate[2]/dTheta;
	}
	else
	{
		rot[0]=0.0;
		rot[1]=0.0;
		rot[2]=0.0;
	}
	updateQuaternion(prop[m].solid.quat, rot, dTheta);
	
	/// Loop for wall particles of solid material "m"
	for(int i=0; i<prop[m].part.n; i++)
	{
		int iPart=prop[m].part.v[i];
		double drGlobal[3], drAux[3];

		/// Transform distance between wall particle and solid CG, Local->Global Coordinate System
		/********** USING QUATERNNION FOR 3D ************/
		drAux[0]=prop[m].solid.dCG.x[i];
		drAux[1]=prop[m].solid.dCG.y[i];
		drAux[2]=prop[m].solid.dCG.z[i];
		tranformRotation3dLocal2Global(drAux, prop[m].solid.quat);

		/// Wall particle displacement, Global Coordinate System
		drGlobal[0]=trans[0]+drAux[0];
		drGlobal[1]=trans[1]+drAux[1];
		drGlobal[2]=trans[2]+drAux[2];

		/// Update Particles Velocities and Positions, Global Coordinate System
		part->un.x[iPart]=(prop[m].solid.cg[0]+drGlobal[0]-part->rn.x[iPart])/sim->dtSolidRigid;
		part->un.y[iPart]=(prop[m].solid.cg[1]+drGlobal[1]-part->rn.y[iPart])/sim->dtSolidRigid;
		part->un.z[iPart]=(prop[m].solid.cg[2]+drGlobal[2]-part->rn.z[iPart])/sim->dtSolidRigid;
		part->rn.x[iPart]=prop[m].solid.cg[0]+drGlobal[0];
		part->rn.y[iPart]=prop[m].solid.cg[1]+drGlobal[1];
		part->rn.z[iPart]=prop[m].solid.cg[2]+drGlobal[2];
	}
	/// Loop for dummy particles of solid material "m"
	for(int i=0; i<prop[m].partDummy.n; i++)
	{
		int iPart=prop[m].partDummy.v[i];
		double drGlobal[3], drAux[3];

		/// Transform distance between dummy particle and solid CG, Local->Global Coordinate System
		/********** USING QUATERNNION FOR 3D ************/
		drAux[0]=prop[m].solid.dCGDummy.x[i];
		drAux[1]=prop[m].solid.dCGDummy.y[i];
		drAux[2]=prop[m].solid.dCGDummy.z[i];
		tranformRotation3dLocal2Global(drAux, prop[m].solid.quat);

		/// Dummy particle displacement, Global Coordinate System
		drGlobal[0]=trans[0]+drAux[0];
		drGlobal[1]=trans[1]+drAux[1];
		drGlobal[2]=trans[2]+drAux[2];

		/// Update Particles Velocities and Positions, Global Coordinate System
		part->un.x[iPart]=(prop[m].solid.cg[0]+drGlobal[0]-part->rn.x[iPart])/sim->dtSolidRigid;
		part->un.y[iPart]=(prop[m].solid.cg[1]+drGlobal[1]-part->rn.y[iPart])/sim->dtSolidRigid;
		part->un.z[iPart]=(prop[m].solid.cg[2]+drGlobal[2]-part->rn.z[iPart])/sim->dtSolidRigid;
		part->rn.x[iPart]=prop[m].solid.cg[0]+drGlobal[0];
		part->rn.y[iPart]=prop[m].solid.cg[1]+drGlobal[1];
		part->rn.z[iPart]=prop[m].solid.cg[2]+drGlobal[2];
	}

	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////
/// BEGINNING - Functions used only if the simulation has collision/contact between rigid bodies

/////////////////////////////////////////////////////////////////////////////////////
/// Set rigid body surface wall features - Edge, Vertex and Face
int calcSolidBC(particles *part, simParameters *sim, matProperties *prop)
{
	/// Loop for all nMat materials
	for(int m=0; m<sim->nMat; m++)
	{
		/// Verify if "m" is a solid (fixed or free) and have particles, i.e., if exists
		if((prop[m].type==solid) && (prop[m].part.n>0))
		{
			/// Loop for wall particles of solid material "m"
#pragma omp parallel for
			for(int n=0; n<prop[m].part.n; n++)
			{
				real dx, dy, dz, rijMod;
				
				/// Stencil position of the nearest neighbors (almost two for 2D and four for 3D)
				real **stencil_pos = new real *[4];
				stencil_pos[0] = new real[3]; // v0
				stencil_pos[1] = new real[3]; // v1
				stencil_pos[2] = new real[3]; // v2
				stencil_pos[3] = new real[3]; // v3
                            
				int i=prop[m].part.v[n];
				int nPos = 0;
				
				for(int l=0; l<part->nNeighS.v[i]; l++)
				{
					int j=part->neighS.v[sim->nNeighMax*i+l];

					// Only wall particles of the solid m
					if(prop[part->index.v[j]].type!=solid||i==j||part->index.v[i]!=part->index.v[j]||part->idMat.v[j]==prop[m].idDummy)
						continue;
					
					dx = part->r.x[j]-part->r.x[i];
					dy = part->r.y[j]-part->r.y[i];
					dz = 0.0;
					if (sim->dim==3)
						dz = part->r.z[j]-part->r.z[i];
					rijMod = sqrt(dx*dx + dy*dy + dz*dz);

					// Stencil vectors (non-dimensional)
					if(rijMod < 1.1*sim->dPart)
					{
						if(nPos < 4)
						{
							stencil_pos[nPos][0] = dx/sim->dPart;
							stencil_pos[nPos][1] = dy/sim->dPart;
							stencil_pos[nPos][2] = dz/sim->dPart;
						}
						nPos++;
					}
				}
				
				if (sim->dim==2) // NEED BE IMPROVED !!!
				{
					// Collinear vector - Cross product = 0
					real cross;
					if(nPos > 2)
					{
						// INTERN EDGE - Wall neighboors > 2 
						part->solidBC.v[i]=edge; // Correspond to the edge 3D
					}
					else if(nPos < 2)
					{
						// EDGE OR VERTEX - Wall neighboors < 2
						part->solidBC.v[i]=edge; // Correspond to the edge 3D
					}
					else
					{
						// FACE OR EDGE
						// v0 x v1
						cross = stencil_pos[0][0]*stencil_pos[1][1] - stencil_pos[0][1]*stencil_pos[1][0];
						if (fabs(cross) <= 1e-3)
							part->solidBC.v[i]=face; // Correspond to the edge 2D
						else
							part->solidBC.v[i]=vertex; // Correspond to the edge 3D
					}
				}
				if (sim->dim==3)
				{
					// Collinear vector - Cross product = 0
					real cross[3];
					if(nPos > 4)
					{
						// INTERN EDGE - Wall neighboors > 4 
						part->solidBC.v[i]=edge;
					}
					else if(nPos < 4)
					{
						// EDGE OR VERTEX - Wall neighboors < 4 - only v0, v1 and v2 are computed
						// v0 x v1
						cross[0] = stencil_pos[0][1]*stencil_pos[1][2] - stencil_pos[0][2]*stencil_pos[1][1];
						cross[1] = stencil_pos[0][2]*stencil_pos[1][0] - stencil_pos[0][0]*stencil_pos[1][2];
						cross[2] = stencil_pos[0][0]*stencil_pos[1][1] - stencil_pos[0][1]*stencil_pos[1][0];

						// v0 x v1 = 0
						if (fabs(cross[0]) <= 1e-3 && fabs(cross[1]) <= 1e-3 && fabs(cross[2]) <= 1e-3)
						{
							part->solidBC.v[i]=edge;
						}
						// vo x v1 != 0   
						else
						{
							// v0 x v2
							cross[0] = stencil_pos[0][1]*stencil_pos[2][2] - stencil_pos[0][2]*stencil_pos[2][1];
							cross[1] = stencil_pos[0][2]*stencil_pos[2][0] - stencil_pos[0][0]*stencil_pos[2][2];
							cross[2] = stencil_pos[0][0]*stencil_pos[2][1] - stencil_pos[0][1]*stencil_pos[2][0];
							
							// v0 x v2 = 0
							if (fabs(cross[0]) <= 1e-3 && fabs(cross[1]) <= 1e-3 && fabs(cross[2]) <= 1e-3)
							{
								part->solidBC.v[i]=edge;
							}
							// v0 x v2 != 0
							else
							{
								// v1 x v2
								cross[0] = stencil_pos[1][1]*stencil_pos[2][2] - stencil_pos[1][2]*stencil_pos[2][1];
								cross[1] = stencil_pos[1][2]*stencil_pos[2][0] - stencil_pos[1][0]*stencil_pos[2][2];
								cross[2] = stencil_pos[1][0]*stencil_pos[2][1] - stencil_pos[1][1]*stencil_pos[2][0];
								
								if (fabs(cross[0]) <= 1e-3 && fabs(cross[1]) <= 1e-3 && fabs(cross[2]) <= 1e-3)
									part->solidBC.v[i]=edge; // v1 x v2 = 0
								else
									part->solidBC.v[i]=vertex; // v1 x v2 != 0
							}
						}
					}
					else
					{
						// FACE OR EDGE
						// v0 x v1
						cross[0] = stencil_pos[0][1]*stencil_pos[1][2] - stencil_pos[0][2]*stencil_pos[1][1];
						cross[1] = stencil_pos[0][2]*stencil_pos[1][0] - stencil_pos[0][0]*stencil_pos[1][2];
						cross[2] = stencil_pos[0][0]*stencil_pos[1][1] - stencil_pos[0][1]*stencil_pos[1][0];

						// v0 x v1 = 0
						if (fabs(cross[0]) <= 1e-3 && fabs(cross[1]) <= 1e-3 && fabs(cross[2]) <= 1e-3)
						{
							// v2 x v3
							cross[0] = stencil_pos[2][1]*stencil_pos[3][2] - stencil_pos[2][2]*stencil_pos[3][1];
							cross[1] = stencil_pos[2][2]*stencil_pos[3][0] - stencil_pos[2][0]*stencil_pos[3][2];
							cross[2] = stencil_pos[2][0]*stencil_pos[3][1] - stencil_pos[2][1]*stencil_pos[3][0];

							if (fabs(cross[0]) <= 1e-3 && fabs(cross[1]) <= 1e-3 && fabs(cross[2]) <= 1e-3)
								part->solidBC.v[i]=face; // v2 x v3 = 0
							else
								part->solidBC.v[i]=edge; // v2 x v3 != 0
						}
						// v0 x v1 != 0
						else
						{
							// v0 x v2
							cross[0] = stencil_pos[0][1]*stencil_pos[2][2] - stencil_pos[0][2]*stencil_pos[2][1];
							cross[1] = stencil_pos[0][2]*stencil_pos[2][0] - stencil_pos[0][0]*stencil_pos[2][2];
							cross[2] = stencil_pos[0][0]*stencil_pos[2][1] - stencil_pos[0][1]*stencil_pos[2][0];
							
							// v0 x v2 = 0
							if (fabs(cross[0]) <= 1e-3 && fabs(cross[1]) <= 1e-3 && fabs(cross[2]) <= 1e-3)
							{
								// v1 x v3
								cross[0] = stencil_pos[1][1]*stencil_pos[3][2] - stencil_pos[1][2]*stencil_pos[3][1];
								cross[1] = stencil_pos[1][2]*stencil_pos[3][0] - stencil_pos[1][0]*stencil_pos[3][2];
								cross[2] = stencil_pos[1][0]*stencil_pos[3][1] - stencil_pos[1][1]*stencil_pos[3][0];

								if (fabs(cross[0]) <= 1e-3 && fabs(cross[1]) <= 1e-3 && fabs(cross[2]) <= 1e-3)
									part->solidBC.v[i]=face; // v1 x v3 = 0
								else
									part->solidBC.v[i]=edge; // v1 x v3 != 0
							}
							// v0 x v2 != 0
							else
							{
								// v0 x v3
								cross[0] = stencil_pos[0][1]*stencil_pos[3][2] - stencil_pos[0][2]*stencil_pos[3][1];
								cross[1] = stencil_pos[0][2]*stencil_pos[3][0] - stencil_pos[0][0]*stencil_pos[3][2];
								cross[2] = stencil_pos[0][0]*stencil_pos[3][1] - stencil_pos[0][1]*stencil_pos[3][0];

								// v0 x v3 = 0
								if (fabs(cross[0]) <= 1e-3 && fabs(cross[1]) <= 1e-3 && fabs(cross[2]) <= 1e-3)
								{
									// v0 x v3 = 0
									// v1 x v2
									cross[0] = stencil_pos[1][1]*stencil_pos[2][2] - stencil_pos[1][2]*stencil_pos[2][1];
									cross[1] = stencil_pos[1][2]*stencil_pos[2][0] - stencil_pos[1][0]*stencil_pos[2][2];
									cross[2] = stencil_pos[1][0]*stencil_pos[2][1] - stencil_pos[1][1]*stencil_pos[2][0];

									if (fabs(cross[0]) <= 1e-3 && fabs(cross[1]) <= 1e-3 && fabs(cross[2]) <= 1e-3)
										part->solidBC.v[i]=face; // v1 x v2 = 0
									else
										part->solidBC.v[i]=edge; // v1 x v2 != 0
								}
								// v0 x v3 != 0
								else
								{
									part->solidBC.v[i]=edge;
								}
							}
						}
					}
				}
				
				for(int ii = 0; ii < 4; ii++) 
				{
					delete [] stencil_pos[ii];
				}
				delete [] stencil_pos;
			}
		}
	}

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Verify if number is divisible
/// https://stackoverflow.com/questions/12097805/how-to-check-if-number-is-divisible-in-c
int IsGoodDivision(int a, int b)
{
	while (b % 2 == 0) { b /= 2; }
	while (b % 5 == 0) { b /= 5; }
	return a % b == 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Compute minium time step used for rigid body. Solid Time step is a multiple of fluid time step
int solidTimeStep(simParameters *sim, particles *part, matProperties *prop, const int m)
{
	/// Loop for wall particles of solid material "m"
	for(int i=0; i<prop[m].part.n; i++)
	{
		int iPart = prop[m].part.v[i];		// Wall particles of the solid
		
		// Find the closest particle "j" to particle "i" with rij < 1.3 sim->dPart
		for(int j=0; j<part->nNeighS.v[iPart]; j++)
		{
			int jPart = part->neighS.v[iPart*sim->nNeighMax+j];
			int mm = part->index.v[jPart];

			// Only free solids that were not considered before in the collision forces and moments calculations (mm>m)
			// or fixed or forced solids
			if( ((prop[mm].type==solid) && (prop[mm].motion==freeMot) && (mm>m)) || ((prop[mm].type==solid) && (prop[mm].motion!=freeMot)) )
			{
				// Solid time step according to the duration of a typical contact/collision
				real dtMin = sim->dtSolidRigid;
				// Elastic modulus
				real Eij = prop[m].solid.youngMod*prop[mm].solid.youngMod /
						((1-pow(prop[m].solid.poisson,2))*prop[mm].solid.youngMod + 
						(1-pow(prop[mm].solid.poisson,2))*prop[m].solid.youngMod);
				// Mass
				real Mij = prop[m].solid.m; // Mass of j equals infinite to fixed or forced solid
				if (prop[mm].motion==freeMot)
					Mij = prop[m].solid.m*prop[mm].solid.m /(prop[m].solid.m + prop[mm].solid.m);
				
				// Relative particle velocity, Global Coordinate System
				real vij[3], vijMod;
				vij[0] = part->un.x[jPart] - part->un.x[iPart];
				vij[1] = part->un.y[jPart] - part->un.y[iPart];
				vij[2] = part->un.z[jPart] - part->un.z[iPart];
				
				vijMod = sqrt(vij[0]*vij[0] + vij[1]*vij[1] + vij[2]*vij[2]);
				
				// If vijMod(t=0) == 0 -> Approximation by velocity due to gravity
				if (isinf(1/vijMod))
					vijMod = sqrt(sim->g[0]*sim->g[0] + sim->g[1]*sim->g[1] + sim->g[2]*sim->g[2])*sim->dt;
				
				// if vij != 0
				if(!isinf(1/vijMod))
				{
					int nstp = 100;
					// dt < 2.87*[M^2/(lo/2*E^2*v)]^1/5 / nstp
					// Based on Hertz’s contact theory. Ensures that at least nstp time steps are used for a typical contact/collision
					dtMin = 2.87*pow(Mij*Mij/(sim->dPart/2*Eij*Eij*vijMod),0.2) / nstp;
				}
				
				// Rigid solid time step multiple of fluid time step
				// After reducing the fraction as much as possible, the denominator can be expressed
				// as a power of 2 multiplied by a power of 5, then the decimal part terminates
				real dtAux, dtsAux;

				dtAux = sim->dt;
				dtsAux = dtMin;

				// Time step to integer number
				while (fmod(dtsAux,1) > 10e-10)
				{
					dtAux = dtAux*10;
					dtsAux = dtsAux*10;
				}

				int int_dt;
				int_dt = int(dtAux);

				sim->nStepsSolid = ceil(sim->dt/dtMin);

				int ok = IsGoodDivision(int_dt,sim->nStepsSolid);

				while (ok==0)
				{
					sim->nStepsSolid++;
					ok = IsGoodDivision(int_dt,sim->nStepsSolid);
				}

				dtMin = sim->dt/sim->nStepsSolid;
				sim->nStepsSolid = (int)(sim->dt/dtMin+0.5);
				
				// Choose the minimum time step
				sim->dtSolidRigid = std::min(dtMin,sim->dtSolidRigid);
			}
		}
	}
	
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Compute collision/contact force and moment (torque) at rigid body "m"
int collisionForceMoment(simParameters *sim, particles *part, matProperties *prop, const int m, double &next_dtSolidRigid, int &next_nStepsSolid)
{
	real dCollision = sim->dPart;	// Particle overlap limit
	real projCol;					// Auxiliar variable
	
	/// Loop for wall particles of solid material "m"
	for(int i=0; i<prop[m].part.n; i++)
	{
		int iPart = prop[m].part.v[i];		// Wall particles of the solid
		int jCol = -1; 						// Particle "j" colliding with "i" then jCol != -1
		real minDist;
		if (sim->dim == 2)
			minDist = 1.119*sim->dPart;		// < sqrt(1.25)*sim->dPart (2D)
		else
			minDist = 1.225*sim->dPart;		// < sqrt(1.50)*sim->dPart (3D)
		int jMin = -1;						// Closest particle ("j") of particle "i"
		real dxMin, dyMin, dzMin, rijMin;
		
		// Find the closest particle "j" to particle "i" with rij < minDIst
		for(int j=0; j<part->nNeighS.v[iPart]; j++)
		{
			real rijMod;
			int jPart = part->neighS.v[iPart*sim->nNeighMax+j];
			int mm = part->index.v[jPart];

			// Only free solids that were not considered before in the collision forces and moments calculations (mm>m)
			// or fixed or forced solids
			if( ((prop[mm].type==solid) && (prop[mm].motion==freeMot) && (mm>m)) || ((prop[mm].type==solid) && (prop[mm].motion!=freeMot)) )
			{
				real dx, dy, dz;
				dx = part->rn.x[jPart]-part->rn.x[iPart];
				dy = part->rn.y[jPart]-part->rn.y[iPart];
				dz = 0.0;
				if (sim->dim==3)
					dz = part->rn.z[jPart]-part->rn.z[iPart];
				rijMod = sqrt(dx*dx + dy*dy + dz*dz);
				
				if (rijMod < minDist)
				{
					minDist = rijMod;
					jMin = jPart;
					dxMin = dx; dyMin = dy; dzMin = dz;
					rijMin = rijMod;
				}
			}
		}
		
		real norm[3]; // Collision normal vector
		norm[0]=norm[1]=norm[2]=0.0;
		
		// If rij < sqrt(1.5)*sim->dPart -> jMin is different of -1
		// Collision normal vector (norm[0],norm[1],norm[2]) and
		// Distance vector (dxMin,dyMin,dzMin) projection on collision normal vector
		if (jMin != -1)
		{
			int mm = part->index.v[jMin];
			
			// Collision normal vector
			// If jMin is face -> Global normal of the neighbor is used as collision normal vector
			if (part->solidBC.v[jMin]==face)
			{
				// Free solid
				if((prop[mm].type==solid) && (prop[mm].motion==freeMot))
				{
					/// Transform Normal, Local->Global Coordinate System
					if (sim->dim==2)
					{
						// USING ANGLE ZZ FOR 2D //
						///	cos()	-sin()
						///	sin()	 cos()
						norm[0] = (part->normal.x[jMin]*cos(prop[mm].solid.orient[2])-part->normal.y[jMin]*sin(prop[mm].solid.orient[2]));
						norm[1] = (part->normal.x[jMin]*sin(prop[mm].solid.orient[2])+part->normal.y[jMin]*cos(prop[mm].solid.orient[2]));
					}
					else
					{
						// USING QUATERNNION FOR 3D //
						norm[0] = part->normal.x[jMin];
						norm[1] = part->normal.y[jMin]; 
						norm[2] = part->normal.z[jMin];
						tranformRotation3dLocal2Global(norm, prop[mm].solid.quat);
					}
				}
				// Fixed wall
				else
				{
					norm[0] = part->normal.x[jMin];
					norm[1] = part->normal.y[jMin]; 
					norm[2] = part->normal.z[jMin];
				}
			}
			// If iPart is face -> The reverse of global normal of iPart is used as collision normal vector
			else if (part->solidBC.v[iPart]==face)
			{
				/// Transform Normal, Local->Global Coordinate System
				if (sim->dim==2)
				{
					// USING ANGLE ZZ FOR 2D //
					///	cos()	-sin()
					///	sin()	 cos()
					norm[0] = -(part->normal.x[iPart]*cos(prop[m].solid.orient[2])-part->normal.y[iPart]*sin(prop[m].solid.orient[2]));
					norm[1] = -(part->normal.x[iPart]*sin(prop[m].solid.orient[2])+part->normal.y[iPart]*cos(prop[m].solid.orient[2]));
				}
				else
				{
					// USING QUATERNNION FOR 3D //
					norm[0] = -part->normal.x[iPart];
					norm[1] = -part->normal.y[iPart]; 
					norm[2] = -part->normal.z[iPart];
					tranformRotation3dLocal2Global(norm, prop[m].solid.quat);
				}
			}
			// Neither jMin or iPart are face -> The reverse of the unit distance between particles is used as collision normal vector
			else
			{
				norm[0] = -dxMin/rijMin;
				norm[1] = -dyMin/rijMin; 
				norm[2] = -dzMin/rijMin;
			}
			
			// Distance vector (dx,dy,dz) projection on collision normal vector
			real proj;
			if(norm[0]==0&&norm[1]==0&&norm[2]==0)
				proj = 0;
			else
				proj = (dxMin*norm[0]+dyMin*norm[1]+dzMin*norm[2])/(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
			
			// If there is collision then |proj| <= dCollision
			if (fabs(proj) <= dCollision && proj != 0)
			{
				jCol = jMin;
				projCol = fabs(proj);
			}
			
			// Solid time step according to the duration of a typical contact/collision
			real dtMin = sim->dtSolidRigid;
			// Elastic modulus
			real Eij = prop[m].solid.youngMod*prop[mm].solid.youngMod /
					((1-pow(prop[m].solid.poisson,2))*prop[mm].solid.youngMod + 
					(1-pow(prop[mm].solid.poisson,2))*prop[m].solid.youngMod);
			// Mass
			real Mij = prop[m].solid.m; // Mass of j equals infinite to fixed or forced solid
			if (prop[mm].motion==freeMot)
				Mij = prop[m].solid.m*prop[mm].solid.m /(prop[m].solid.m + prop[mm].solid.m);
			
			// Relative particle velocity, Global Coordinate System
			real vij[3], vijMod;
			vij[0] = part->un.x[jMin] - part->un.x[iPart];
			vij[1] = part->un.y[jMin] - part->un.y[iPart];
			vij[2] = part->un.z[jMin] - part->un.z[iPart];
			
			vijMod = sqrt(vij[0]*vij[0] + vij[1]*vij[1] + vij[2]*vij[2]);
			
			// if vij != 0
			if(!isinf(1/vijMod))
			{
				int nstp = 100;
				// dt < 2.87*[M^2/(lo/2*E^2*v)]^1/5 / nstp
				// Based on Hertz’s contact theory. Ensures that at least nstp time steps are used for a typical contact/collision
				dtMin = 2.87*pow(Mij*Mij/(sim->dPart/2*Eij*Eij*vijMod),0.2) / nstp;
			}
			
			// Rigid solid time step multiple of fluid time step
			// After reducing the fraction as much as possible, the denominator can be expressed
			// as a power of 2 multiplied by a power of 5, then the decimal part terminates
			real dtAux, dtsAux;

			dtAux = sim->dt;
			dtsAux = dtMin;

			// Time step to integer number
			while (fmod(dtsAux,1) > 10e-10)
			{
				dtAux = dtAux*10;
				dtsAux = dtsAux*10;
			}

			int int_dt;
			int_dt = int(dtAux);

			next_nStepsSolid = ceil(sim->dt/dtMin);

			int ok = IsGoodDivision(int_dt,next_nStepsSolid);

			while (ok==0)
			{
				next_nStepsSolid++;
				ok = IsGoodDivision(int_dt,next_nStepsSolid);
			}

			dtMin = sim->dt/next_nStepsSolid;
			next_nStepsSolid = (int)(sim->dt/dtMin+0.5);
			
			// Choose the minimum time step
			next_dtSolidRigid = std::min(dtMin,sim->dtSolidRigid);
		}
		else
			continue;
				
		// If "jCol" is different of -1 -> there is collision 
		if (jCol != -1)
		{
			// Pressure of collision particles is set to zero
			part->p.v[iPart] = 0.0;
			part->p.v[jCol] = 0.0;
			
			int mm = part->index.v[jCol];
			
			real springForce[3], dampingForce[3], shearForce[3], collisionMoment[3];
			springForce[0] = springForce[1] = springForce[2] = 0.0;
			dampingForce[0] = dampingForce[1] = dampingForce[2] = 0.0;
			shearForce[0] = shearForce[1] = shearForce[2] = 0.0;
			collisionMoment[0] = collisionMoment[1] = collisionMoment[2] = 0.0;
			
			real tang[3]; // Collision tangential vector
			tang[0] = tang[1] = tang[2] = 0.0;
			
			real dCG[3], PCij[3], vi[3], vj[3], vijt[3], proj; // rij[3], drAuxJ[3], dCGi[3], dCGj[3];

			// Elastic modulus
			real Eij = prop[m].solid.youngMod*prop[mm].solid.youngMod /
					((1-pow(prop[m].solid.poisson,2))*prop[mm].solid.youngMod + 
					(1-pow(prop[mm].solid.poisson,2))*prop[m].solid.youngMod);
			// Normal Stiffness constant
			real Knij = 4.0*Eij*sqrt(sim->dPart/2.0)/3.0;
			// Tangential Stiffness constant
			real Ktij = 2.0*Knij/7.0;
			// Mass
			real Mij = prop[m].solid.m; // Mass of j equals infinite to fixed or forced solid
			if (prop[mm].motion==freeMot)
				Mij = prop[m].solid.m*prop[mm].solid.m/(prop[m].solid.m + prop[mm].solid.m);
			// Normal Damping coefficient
			real Cnij = prop[m].solid.damping*sqrt(6.0*Mij*Eij*sqrt(sim->dPart/2));
			// Tangential Damping coefficient
			real Ctij = 2.0*Cnij/7.0;
			// Friction
			real friction = prop[mm].solid.friction; // Friction of fixed or forced solid
			if (prop[mm].motion==freeMot)
				friction = (prop[m].solid.friction+prop[mm].solid.friction)/2.0; // Ensures action and reaction
			
			// Collision forces and moments 
			if(sim->dim==2)
			{
				/// Spring force (Non-linear Hertzian model), Global Coordinate System
				springForce[0]=Knij*pow(dCollision-projCol,1.5)*norm[0];
				springForce[1]=Knij*pow(dCollision-projCol,1.5)*norm[1];
				
				/// Force mm->m
				prop[m].solid.globalCollisionForce[0] += springForce[0];
				prop[m].solid.globalCollisionForce[1] += springForce[1];
				/// Force m->mm
				prop[mm].solid.globalCollisionForce[0] -= springForce[0];
				prop[mm].solid.globalCollisionForce[1] -= springForce[1];
				
				/// Particle velocity, Global Coordinate System
				vi[0] = part->un.x[iPart];
				vi[1] = part->un.y[iPart];
				vj[0] = part->un.x[jCol];
				vj[1] = part->un.y[jCol];

				/// Damping force (Non-linear Hertzian model), Global Coordinate System
				dampingForce[0] = Cnij*pow(dCollision-projCol,0.25)*(vj[0]-vi[0]);
				dampingForce[1] = Cnij*pow(dCollision-projCol,0.25)*(vj[1]-vi[1]);

				if(norm[0]==0&&norm[1]==0)
					proj = 0;
				else
					proj = (dampingForce[0]*norm[0]+dampingForce[1]*norm[1])/(norm[0]*norm[0]+norm[1]*norm[1]);
				
				dampingForce[0] = proj*norm[0];
				dampingForce[1] = proj*norm[1];
				
				/// Force mm->m
				prop[m].solid.globalCollisionForce[0] += dampingForce[0];
				prop[m].solid.globalCollisionForce[1] += dampingForce[1];
				/// Force m->mm
				prop[mm].solid.globalCollisionForce[0] -= dampingForce[0];
				prop[mm].solid.globalCollisionForce[1] -= dampingForce[1];
				
				/// Tangential unit vector
				real vijtMod = 0;
				if(norm[0]==0&&norm[1]==0)
				{
					tang[0]=tang[1]=0;
					vijt[0] = (vj[0]-vi[0]);
					vijt[1] = (vj[1]-vi[1]);
				}
				else
				{
					proj = ((vj[0]-vi[0])*norm[0]+(vj[1]-vi[1])*norm[1])/(norm[0]*norm[0]+norm[1]*norm[1]);
					vijt[0] = (vj[0]-vi[0])-proj*norm[0];
					vijt[1] = (vj[1]-vi[1])-proj*norm[1];
					vijtMod = sqrt(vijt[0]*vijt[0]+vijt[1]*vijt[1]);
					if (vijtMod==0)
						tang[0]=tang[1]=0;
					else
					{
						tang[0] = vijt[0]/vijtMod;
						tang[1] = vijt[1]/vijtMod;
					}
				}
				
				/// Shear force, Global Coordinate System
				/// Modified Coulomb model
				real NormalForce[2], normalForceMod, CoulombForce[2], CoulombForceMod;
				NormalForce[0] = springForce[0] + dampingForce[0];
				NormalForce[1] = springForce[1] + dampingForce[1];
				normalForceMod = sqrt(NormalForce[0]*NormalForce[0] + NormalForce[1]*NormalForce[1]);
				
				// If simulation is unstable, maybe the problem could be the friction. Then, uncoment tanh(8*vijtMod) below
				CoulombForce[0] = friction*normalForceMod*tanh(8.0*vijtMod)*tang[0];
				CoulombForce[1] = friction*normalForceMod*tanh(8.0*vijtMod)*tang[1];
				CoulombForceMod = sqrt(CoulombForce[0]*CoulombForce[0] + CoulombForce[1]*CoulombForce[1]);
				
				/// Linear Mindlin model
				real MindlinForce[2], MindlinMod, rijt[2];
				// Tangential particle overlap
				rijt[0] = vijt[0]*sim->dtSolidRigid;
				rijt[1] = vijt[1]*sim->dtSolidRigid;

				MindlinForce[0] = Ktij*rijt[0] + Ctij*vijt[0];
				MindlinForce[1] = Ktij*rijt[1] + Ctij*vijt[1];
				MindlinMod = sqrt(MindlinForce[0]*MindlinForce[0] + MindlinForce[1]*MindlinForce[1]);
				
				if (CoulombForceMod < MindlinMod)
				{
					shearForce[0] = CoulombForce[0];
					shearForce[1] = CoulombForce[1];
				}
				else
				{
					shearForce[0] = MindlinForce[0];
					shearForce[1] = MindlinForce[1];
				}
				
				/// Force mm->m
				prop[m].solid.globalCollisionForce[0] += shearForce[0];
				prop[m].solid.globalCollisionForce[1] += shearForce[1];
				/// Force m->mm
				prop[mm].solid.globalCollisionForce[0] -= shearForce[0];
				prop[mm].solid.globalCollisionForce[1] -= shearForce[1];
				
				/// Moment, Global Coordinate System
				// Intermediary point of contact
				PCij[0] = (part->rn.x[iPart]+part->rn.x[jCol])/2.0;
				PCij[1] = (part->rn.y[iPart]+part->rn.y[jCol])/2.0;
				
				/// Moment mm->m
				dCG[0] = PCij[0]-prop[m].solid.cg[0];
				dCG[1] = PCij[1]-prop[m].solid.cg[1];
				
				collisionMoment[2] = dCG[0]*(springForce[1]+dampingForce[1]+shearForce[1])-dCG[1]*(springForce[0]+dampingForce[0]+shearForce[0]);
				
				prop[m].solid.globalCollisionMoment[2] += collisionMoment[2];
				
				/// Moment m->mm
				dCG[0] = PCij[0]-prop[mm].solid.cg[0];
				dCG[1] = PCij[1]-prop[mm].solid.cg[1];
				
				collisionMoment[2] = dCG[0]*(springForce[1]+dampingForce[1]+shearForce[1])-dCG[1]*(springForce[0]+dampingForce[0]+shearForce[0]);
				
				prop[mm].solid.globalCollisionMoment[2] -= collisionMoment[2];

			}
			else
			{
				/// Spring force (Non-linear Hertzian model), Global Coordinate System
				springForce[0] = Knij*pow(dCollision-projCol,1.5)*norm[0];
				springForce[1] = Knij*pow(dCollision-projCol,1.5)*norm[1];
				springForce[2] = Knij*pow(dCollision-projCol,1.5)*norm[2];
				
				/// Force mm->m
				prop[m].solid.globalCollisionForce[0] += springForce[0];
				prop[m].solid.globalCollisionForce[1] += springForce[1];
				prop[m].solid.globalCollisionForce[2] += springForce[2];
				/// Force m->mm
				prop[mm].solid.globalCollisionForce[0] -= springForce[0];
				prop[mm].solid.globalCollisionForce[1] -= springForce[1];
				prop[mm].solid.globalCollisionForce[2] -= springForce[2];
				
				/// Particle velocity, Global Coordinate System
				vi[0] = part->un.x[iPart];
				vi[1] = part->un.y[iPart];
				vi[2] = part->un.z[iPart];
				vj[0] = part->un.x[jCol];
				vj[1] = part->un.y[jCol];
				vj[2] = part->un.z[jCol];

				/// Damping force (Non-linear Hertzian model), Global Coordinate System
				dampingForce[0] = Cnij*pow(dCollision-projCol,0.25)*(vj[0]-vi[0]);
				dampingForce[1] = Cnij*pow(dCollision-projCol,0.25)*(vj[1]-vi[1]);
				dampingForce[2] = Cnij*pow(dCollision-projCol,0.25)*(vj[2]-vi[2]);
				
				if(norm[0]==0&&norm[1]==0&&norm[2]==0)
					proj = 0;
				else
					proj = (dampingForce[0]*norm[0]+dampingForce[1]*norm[1]+dampingForce[2]*norm[2])/(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
				
				dampingForce[0] = proj*norm[0];
				dampingForce[1] = proj*norm[1];
				dampingForce[2] = proj*norm[2];
				
				/// Force mm->m
				prop[m].solid.globalCollisionForce[0] += dampingForce[0];
				prop[m].solid.globalCollisionForce[1] += dampingForce[1];
				prop[m].solid.globalCollisionForce[2] += dampingForce[2];
				/// Force m->mm
				prop[mm].solid.globalCollisionForce[0] -= dampingForce[0];
				prop[mm].solid.globalCollisionForce[1] -= dampingForce[1];
				prop[mm].solid.globalCollisionForce[2] -= dampingForce[2];
				
				/// Tangential unit vector
				real vijtMod = 0;
				if(norm[0]==0&&norm[1]==0&&norm[2]==0)
				{
					tang[0]=tang[1]=tang[2]=0;
					vijt[0] = (vj[0]-vi[0]);
					vijt[1] = (vj[1]-vi[1]);
					vijt[2] = (vj[2]-vi[2]);
				}
				else
				{
					proj = ((vj[0]-vi[0])*norm[0]+(vj[1]-vi[1])*norm[1]+(vj[2]-vi[2])*norm[2])/(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
					vijt[0] = (vj[0]-vi[0])-proj*norm[0];
					vijt[1] = (vj[1]-vi[1])-proj*norm[1];
					vijt[2] = (vj[2]-vi[2])-proj*norm[2];
					vijtMod = sqrt(vijt[0]*vijt[0]+vijt[1]*vijt[1]+vijt[2]*vijt[2]);
					if (vijtMod==0)
						tang[0]=tang[1]=tang[2]=0;
					else
					{
						tang[0] = vijt[0]/vijtMod;
						tang[1] = vijt[1]/vijtMod;
						tang[2] = vijt[2]/vijtMod;
					}
				}
				
				/// Shear force, Global Coordinate System
				/// Modified Coulomb model
				real NormalForce[3], normalForceMod, CoulombForce[3], CoulombForceMod;
				NormalForce[0] = springForce[0] + dampingForce[0];
				NormalForce[1] = springForce[1] + dampingForce[1];
				NormalForce[2] = springForce[2] + dampingForce[2];
				normalForceMod = sqrt(NormalForce[0]*NormalForce[0] + NormalForce[1]*NormalForce[1] + NormalForce[2]*NormalForce[2]);
				
				// Sliding - comment tanh
				// If simulation is unstable, maybe the problem could be the friction. Then, uncoment tanh(8*vijtMod) below
				CoulombForce[0] = friction*normalForceMod*tanh(8.0*vijtMod)*tang[0];
				CoulombForce[1] = friction*normalForceMod*tanh(8.0*vijtMod)*tang[1];
				CoulombForce[2] = friction*normalForceMod*tanh(8.0*vijtMod)*tang[2];
				CoulombForceMod = sqrt(CoulombForce[0]*CoulombForce[0] + CoulombForce[1]*CoulombForce[1] + CoulombForce[2]*CoulombForce[2]);
				
				/// Linear Mindlin model
				real MindlinForce[3], MindlinMod, rijt[3];
				// Tangential particle overlap
				rijt[0] = vijt[0]*sim->dtSolidRigid;
				rijt[1] = vijt[1]*sim->dtSolidRigid;
				rijt[2] = vijt[2]*sim->dtSolidRigid;

				MindlinForce[0] = Ktij*rijt[0] + Ctij*vijt[0];
				MindlinForce[1] = Ktij*rijt[1] + Ctij*vijt[1];
				MindlinForce[2] = Ktij*rijt[2] + Ctij*vijt[2];
				MindlinMod = sqrt(MindlinForce[0]*MindlinForce[0] + MindlinForce[1]*MindlinForce[1] + MindlinForce[2]*MindlinForce[2]);
				
				if (CoulombForceMod < MindlinMod)
				{
					shearForce[0] = CoulombForce[0];
					shearForce[1] = CoulombForce[1];
					shearForce[2] = CoulombForce[2];
				}
				else
				{
					shearForce[0] = MindlinForce[0];
					shearForce[1] = MindlinForce[1];
					shearForce[2] = MindlinForce[2];
				}
				
				/// Force mm->m
				prop[m].solid.globalCollisionForce[0] += shearForce[0];
				prop[m].solid.globalCollisionForce[1] += shearForce[1];
				prop[m].solid.globalCollisionForce[2] += shearForce[2];
				/// Force m->mm
				prop[mm].solid.globalCollisionForce[0] -= shearForce[0];
				prop[mm].solid.globalCollisionForce[1] -= shearForce[1];
				prop[mm].solid.globalCollisionForce[2] -= shearForce[2];
				
				/// Moment, Global Coordinate System
				// Intermediary point of contact
				PCij[0] = (part->rn.x[iPart]+part->rn.x[jCol])/2.0;
				PCij[1] = (part->rn.y[iPart]+part->rn.y[jCol])/2.0;
				PCij[2] = (part->rn.z[iPart]+part->rn.z[jCol])/2.0;
				
				/// Moment mm->m
				dCG[0] = PCij[0]-prop[m].solid.cg[0];
				dCG[1] = PCij[1]-prop[m].solid.cg[1];
				dCG[2] = PCij[2]-prop[m].solid.cg[2];
				
				collisionMoment[0] = dCG[1]*(springForce[2]+dampingForce[2]+shearForce[2])-dCG[2]*(springForce[1]+dampingForce[1]+shearForce[1]);
				collisionMoment[1] = dCG[2]*(springForce[0]+dampingForce[0]+shearForce[0])-dCG[0]*(springForce[2]+dampingForce[2]+shearForce[2]);
				collisionMoment[2] = dCG[0]*(springForce[1]+dampingForce[1]+shearForce[1])-dCG[1]*(springForce[0]+dampingForce[0]+shearForce[0]);
				
				prop[m].solid.globalCollisionMoment[0] += collisionMoment[0];
				prop[m].solid.globalCollisionMoment[1] += collisionMoment[1];
				prop[m].solid.globalCollisionMoment[2] += collisionMoment[2];
				
				/// Moment m->mm
				dCG[0] = PCij[0]-prop[mm].solid.cg[0];
				dCG[1] = PCij[1]-prop[mm].solid.cg[1];
				dCG[2] = PCij[2]-prop[mm].solid.cg[2];
				
				collisionMoment[0] = dCG[1]*(springForce[2]+dampingForce[2]+shearForce[2])-dCG[2]*(springForce[1]+dampingForce[1]+shearForce[1]);
				collisionMoment[1] = dCG[2]*(springForce[0]+dampingForce[0]+shearForce[0])-dCG[0]*(springForce[2]+dampingForce[2]+shearForce[2]);
				collisionMoment[2] = dCG[0]*(springForce[1]+dampingForce[1]+shearForce[1])-dCG[1]*(springForce[0]+dampingForce[0]+shearForce[0]);
				
				prop[mm].solid.globalCollisionMoment[0] -= collisionMoment[0];
				prop[mm].solid.globalCollisionMoment[1] -= collisionMoment[1];
				prop[mm].solid.globalCollisionMoment[2] -= collisionMoment[2];
			}
		}
	}
	
	/// Local forces and moments
	if(sim->dim==2)
	{
		/// Transform collision force, Global->Local Coordinate System
		/********** USING ANGLE ZZ FOR 2D ************/
		///	cos()	 sin()
		///	-sin()	 cos()
		prop[m].solid.collisionForce[0] = prop[m].solid.globalCollisionForce[0]*cos(prop[m].solid.orient[2])+prop[m].solid.globalCollisionForce[1]*sin(prop[m].solid.orient[2]);
		prop[m].solid.collisionForce[1] = -prop[m].solid.globalCollisionForce[0]*sin(prop[m].solid.orient[2])+prop[m].solid.globalCollisionForce[1]*cos(prop[m].solid.orient[2]);
		
		prop[m].solid.collisionMoment[2] = prop[m].solid.globalCollisionMoment[2];
	}
	else // if(dim==3)
	{
		for (int k = 0; k < 3; k++) 
		{
			prop[m].solid.collisionForce[k] = prop[m].solid.globalCollisionForce[k];
			prop[m].solid.collisionMoment[k] = prop[m].solid.globalCollisionMoment[k];
		}
		/// Transform collision force, Global->Local Coordinate System
		/********** USING QUATERNNION FOR 3D ************/
		tranformRotation3dGlobal2Local(prop[m].solid.collisionForce, prop[m].solid.quat);
		tranformRotation3dGlobal2Local(prop[m].solid.collisionMoment, prop[m].solid.quat);
	}
	
	return 0;
}

/// END - Functions used only if the simulation has collision/contact between rigid bodies
/////////////////////////////////////////////////////////////////////////////////////