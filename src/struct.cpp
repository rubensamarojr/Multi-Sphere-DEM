/** MPS-Moving Particle Semi-implicit
 * ***(c) 2008 USP/TPN, NYXknowledge
 *** Marcio Tsukamoto, Guilherme Rueda, Mariana Robortella, 
 *** Rubens A. Amaro Jr, Cheng Liang Yee, et al.
 **/

#include "struct.h"
#include <iostream>
#include <mkl.h>
#include <mkl_blas.h>

/////////////////////////////////////////////////////////////////////////////////////////
void makeIntArray(intArray *a, const integer n)
{
	a->v=(integer *) MKL_malloc(n*sizeof(integer), 128);
}

void freeIntArray(intArray *a)
{
	MKL_free(a->v);
}

void zeroIntArray(intArray *a, const integer n)
{
#pragma omp for
	for(int i=0; i<n; i++)
		a->v[i]=0;
}

/////////////////////////////////////////////////////////////////////////////////////////
void makeRealArray(realArray *a, const integer n)
{
	a->v=(real *) MKL_malloc(n*sizeof(real), 128);
}
void freeRealArray(realArray *a)
{
	MKL_free(a->v);
}
void zeroRealArray(realArray *a, const integer n)
{
#pragma omp for
	for(int i=0; i<n; i++)
		a->v[i]=0.0;
}

/////////////////////////////////////////////////////////////////////////////////////////
void makeVector(vector *v, const integer n)
{
	v->n=n;
	v->nMax=n;
	v->x=(real *) MKL_malloc(v->n*sizeof(real), 128);
	v->y=(real *) MKL_malloc(v->n*sizeof(real), 128);
	v->z=(real *) MKL_malloc(v->n*sizeof(real), 128);
	zeroVector(v, n);
}
void freeVector(vector *v)
{
	MKL_free(v->x);
	MKL_free(v->y);
	MKL_free(v->z);
	v->n=0;
}
void zeroVector(vector *v, const integer n)
{
#pragma omp parallel for
	for(int i=0; i<n; i++)
	{
		v->x[i]=0;
		v->y[i]=0;
		v->z[i]=0;
	}
}

void particles::Alloc(integer nSize, integer nNeighMax)
{
	makeIntArray(&index, nSize);
	makeIntArray(&idMat, nSize);
	makeIntArray(&neighS, nSize*nNeighMax);
	makeIntArray(&neighL, nSize*nNeighMax);
	makeIntArray(&nNeighS, nSize);
	makeIntArray(&nNeighL, nSize);
	makeIntArray(&solidBC, nSize);
	makeRealArray(&p, nSize);
	makeRealArray(&pndS, nSize);
	makeRealArray(&pndL, nSize);
	makeRealArray(&pndS_mat, nSize);
	makeRealArray(&dNeighL, nSize*nNeighMax);
	makeRealArray(&dNeighS, nSize*nNeighMax);
	makeRealArray(&wL, nSize*nNeighMax);
	makeRealArray(&wS, nSize*nNeighMax);
	makeVector(&r, nSize);
	makeVector(&rn, nSize);
	makeVector(&dr, nSize);
	makeVector(&u, nSize);
	makeVector(&un, nSize);
	makeVector(&du, nSize);
	makeVector(&normal, nSize);

	zeroIntArray(&index, nSize);
	zeroIntArray(&idMat, nSize);
	zeroIntArray(&neighS, nSize*nNeighMax);
	zeroIntArray(&neighL, nSize*nNeighMax);
	zeroIntArray(&nNeighS, nSize);
	zeroIntArray(&nNeighL, nSize);
	zeroIntArray(&solidBC, nSize);
	zeroRealArray(&p, nSize);
	zeroRealArray(&pndS, nSize);
	zeroRealArray(&pndL, nSize);
	zeroRealArray(&pndS_mat, nSize);
	zeroRealArray(&dNeighL, nSize*nNeighMax);
	zeroRealArray(&dNeighS, nSize*nNeighMax);
	zeroRealArray(&wL, nSize*nNeighMax);
	zeroRealArray(&wS, nSize*nNeighMax);
	zeroVector(&r, nSize);
	zeroVector(&rn, nSize);
	zeroVector(&dr, nSize);
	zeroVector(&u, nSize);
	zeroVector(&un, nSize);
	zeroVector(&du, nSize);
	zeroVector(&normal, nSize);
}

void particles::Free()
{
	freeIntArray(&index);
	freeIntArray(&idMat);
	freeIntArray(&neighS);
	freeIntArray(&neighL);
	freeIntArray(&nNeighS);
	freeIntArray(&nNeighL);
	freeIntArray(&solidBC);
	freeRealArray(&p);
	freeRealArray(&pndS);
	freeRealArray(&pndL);
	freeRealArray(&pndS_mat);
	freeRealArray(&dNeighL);
	freeRealArray(&dNeighS);
	freeRealArray(&wL);
	freeRealArray(&wS);
	freeVector(&r);
	freeVector(&rn);
	freeVector(&dr);
	freeVector(&u);
	freeVector(&un);
	freeVector(&du);
	freeVector(&normal);
}

/// Materials
void initProperties(simParameters *sim, matProperties * &prop)
{
	prop=new matProperties[sim->nMat];

	for(int m=0; m<sim->nMat; m++)
	{
		/// User input. Can be read from a input file
		prop[m].idDummy=-1;			// Dummy particles material ID
		prop[m].part.n=1000;		// Number of Fluid or Wall particles belonging to the material "m"
		prop[m].partDummy.n=200;	// Number of Dummy particles belonging to the material "m"
		prop[m].type=53;			// flow=50, fluid=51, solid=53, boundCond=62
		prop[m].motion=200;			// freeMot=200, forcedMot=201, fixedMot=202
	
		/// Verify if "m" is a solid (fixed or free) and have particles, i.e., if exists
		if((prop[m].type==solid) && (prop[m].part.n>0))
		{
			prop[m].solid.m=10.0;
			prop[m].solid.poisson=0.0;
			prop[m].solid.youngMod=0.0;
			prop[m].solid.damping=0.0;

			/// Verify if "m" is a free solid
			if(prop[m].motion==freeMot)
			{
				// Inertia matrix in Local coordinate system
				prop[m].solid.I[0][0]=1.0;		
				prop[m].solid.I[0][1]=0.0;
				prop[m].solid.I[0][2]=0.0;
				prop[m].solid.I[1][0]=0.0;
				prop[m].solid.I[1][1]=1.0;
				prop[m].solid.I[1][2]=0.0;
				prop[m].solid.I[2][0]=0.0;
				prop[m].solid.I[2][1]=0.0;
				prop[m].solid.I[2][2]=1.0;

				prop[m].solid.cg[0]=0.0;		// CGx in Local coordinate system
				prop[m].solid.cg[1]=0.0;		// CGy in Local coordinate system
				prop[m].solid.cg[2]=0.0;		// CGz in Local coordinate system
				prop[m].solid.orient[0]=0.0;	// ANGLEx orientation in Local coordinate system
				prop[m].solid.orient[1]=0.0;	// ANGLEy orientation in Local coordinate system
				prop[m].solid.orient[2]=0.0;	// ANGLEz orientation in Local coordinate system
				for(int k=0; k<6; k++)
				{
					prop[m].solid.freeDir[k]=1.0;	// Motion directions restrictions (Free:1, Fixed: 0) in Local Coordinate System
					prop[m].solid.followDir[k]=0.0;	// Rigid body follow prescribed motion (Free:1, Fixed: 0) in Local Coordinate System
				}
			}
		}
	}
}