/** MPS-Moving Particle Semi-implicit
 * ***(c) 2008 USP/TPN, NYXknowledge
 *** Marcio Tsukamoto, Guilherme Rueda, Mariana Robortella, 
 *** Rubens A. Amaro Jr, Cheng Liang Yee, et al.
 **/

#ifndef __STRUCT_H__
#define __STRUCT_H__

#include <string>
#include <cassert>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

////////////////////////////////////////////////////////////////////////////////////////
/// Weight function
inline real mps_weight(real r, real re)
{
	return (re>r)?(re/r-1.0):0.0;
}

/////////////////////////////////////////////////////////////////////////////////////////
/// Type definitions
int const ghost=123456;					/// material id (idMat)
enum materialEnumType { flow=50, fluid=51, solid=53, boundCond=62 };
enum motionEnumType
{
	freeMot=200, forcedMot=201, fixedMot=202
};

/// Rigid body
struct solidProperties
{
	real M[6][6];						/// Mass matrix
	real m, I[3][3];					/// Mass, Inertia matrix, Local coordinate system (freeSolid)
	real cg[3], cg0[3];					/// cg: center of gravity (CG), Local coordinate system (freeSolid)
	real orient[3];						/// orient: CG angular orientation, Local coordinate system (freeSolid)
	real poisson, youngMod, damping;	/// Poisson's ratio, Young's modulus and damping (Solid collision)
	real hydroForce[3], hydroMoment[3];	/// Force and moment (torque) due Hydrodynamic pressure, Local coordinate system (freeSolid)
	real globalHydroForce[3], globalHydroMoment[3];	/// Force and moment (torque) due Hydrodynamic pressure, Global coordinate system (freeSolid)
	real f1[3], f2[3];					/// f1: force(k-1), f2: force(k-2) (freeSolid)
	real m1[3], m2[3];					/// m1: moment(k-1), m2: moment(k-2) (freeSolid)
	real vel[3], velAng[4];				/// vel: translational velocity, velAng: angular velocity (freeSolid, forceSolid)
	real vel1[3], vel2[3];				/// Auxiliar velocity vel1: vel(k-1), vel2: vel(k-2)
	real velAng1[4], velAng2[4];		/// Auxiliar velocity velAng1: velAng(k-1), velAng2: velAng(k-2)
	real freeDir[6], followDir[6];		/// Motion directions restrictions
	realArray areaA;					/// particle projected area (freeSolid)
	vector dCG, dCGDummy;				/// distance between particles and the center of gravity (freeSolid)
	vector normalA;						/// Normal vector at each solid wall particle (freeSolid)
	real quat[4];						/// quaternion
	real collisionForce[3], collisionMoment[3];				/// Force and moment (torque) due Collision between rigid bodies, Local coordinate system
	real globalCollisionForce[3], globalCollisionMoment[3];	/// Force and moment (torque) due Collision between rigid bodies,  Global coordinate system
	real constantLoad[6], constantLoadPos[3];				///	External applied cte load
};

/// Particles
struct particles
{
	intArray index;						/// Specifies the material index
	intArray idMat;						/// Specifies the material id
	intArray neighS, neighL;			/// Neighbor particles list (small radius/large radius)
	intArray nNeighS, nNeighL;			/// Number of neighbor particles (small radius/large radius)
	intArray solidBC;					/// Solid boundary type
	realArray p;						/// Pressure
	realArray pndS, pndL, pndS_mat;		/// Particle number density (small radius/large radius)
	realArray dNeighS, dNeighL;			/// r[i]-r[j] (small radius/large radius)
	realArray wS, wL;					/// Weight function (small radius/large radius)
	vector r, rn, dr;					/// Position r(k) and r(*), and delta position
	vector u, un, du;					/// Velocities u(k), u(*) and u', and delta velocity
	vector normal;						/// Normal vector for solid wall particles
	
	void Alloc(integer nSize, integer nNeighMax);	/// Alloc memory to particles
};

// Simulation
struct simParameters
{
	integer dim;				/// Dimension (2 or 3)
	integer nIter;				/// Simulation maximum number of iterations
	integer nMat, nMaxMat;		/// Number of particles types (materials)
	integer nPart;				/// Number of particles, max number of particles and initial val
	integer nNeighMax;			/// Max number of neighbors
	integer nStepsSolid;			/// Number of Subcycling steps (>= 1) used for rigid body
	real g[3];					/// Gravity
	real reS, reS2, reL, reL2;	/// Effective radius
	real dPart;					/// Initial adjacent particle distance
	real t, dt;					/// Actual time, Actual time-step
	real tStart, tFinal;		/// Start and End time
	real dtSolidRigid;			/// Time-step rigid body
};

/// Materials
struct matProperties
{
	solidProperties solid;			/// Rigid body structure
	integer type;					/// Material type
	integer motion;					/// Motion type
	integer idDummy;				/// Dummy particles material ID
	intArray part, partDummy;		/// Particle and dummy indexes (freeSolid, forcedSolid, inflow, outflow)
};


/////////////////////////////////////////////////////////////////////////////////////////
/// MPS Function prototypes
int CreateNeighborList(particles *part, simParameters *sim, matProperties *prop);
int CalcPND(realArray *pnd, realArray *pnd1, intArray *nNeigh, intArray *neigh, realArray *dNeigh, realArray *w, real re, simParameters *sim, matProperties *prop, particles *part);
int CalcPND(realArray *pnd, intArray *nNeigh, intArray *neigh, realArray *dNeigh, realArray *w, real re, simParameters *sim, matProperties *prop);
int CheckBoundCond(particles *part, simParameters *sim, matProperties *prop);
int OutputData(particles *part, simParameters *sim, matProperties *prop, long iter);

////////////////////////////////////////////////////////////////////////////////////////
/// Array and vector
/// Type: real
typedef double real;
/// Type: real
typedef long integer;
/// Type: integer
typedef int integer32;
/// Type: array
struct intArray
{
	integer n;
	integer nMax;
	integer *v;
};
struct realArray
{
	long n;
	real *v;
};
// Type: vector
/struct vector
{
	long n;
	integer nMax;
	real *x, *y, *z;
};
/// Type: mask
typedef realArray mask;

/// Array Creation/Destruction
void makeIntArray(intArray *a, const integer n);
void freeIntArray(intArray *a);
void zeroIntArray(intArray *a, const integer n);

void makeRealArray(realArray *a, const integer n);
void freeRealArray(realArray *a);
void zeroRealArray(realArray *a, const integer n);

// Vector Creation/Destruction
void makeVector(vector *v, const integer n);
void freeVector(vector *v);
void zeroVector(vector *v, const integer n);

#endif /* __STRUCT_H__ */