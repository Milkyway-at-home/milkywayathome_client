/* This program tests the generation of a 2-component dwarf galaxy.
* It creates a plummer-plummer dwarf and an NFW dwarf.
* The structure and stability of these two dwarfs are checked by calculating their virial ratios and checking their center of masses.
*/

#include "nbody_mixeddwarf.h"
#include "nbody_dwarf_potential.h"
#include "dSFMT.h"
#include <stdio.h>


static inline real doubleCompDensity(const Dwarf* comp1, const Dwarf* comp2, real r)
{
	//Function calls from nbody_dwarf_potential
	return get_density(comp1, r) + get_density(comp2, r);
}

static inline real doubleCompPotential(const Dwarf* comp1, const Dwarf* comp2, real r)
{
	//Function calls from nbody_dwarf_potential
	return get_potential(comp1, r) + get_potential(comp2, r);
}

//Computes the kinetic and potential energies and ensures that it obeys the virial theorem
int checkVirialRatio(const Dwarf* comp1, const Dwarf* comp2, const mwvector* pos, const mwvector* vel, const real* mass, unsigned int numBodies)
{
	int failed = 0;
	const real TOL = 0.05;

	real T  = 0; //Kinetic energy
	real U1 = 0; //Potential to be calculated using the dwarf's potential
	real U2 = 0; //potential to be calculated using newtonian gravity

	//now doing the math
	for(unsigned int i = 0; i < numBodies; i++)
	{
		real x  = pos[i].x;
		real y  = pos[i].y;
		real z  = pos[i].z;
		real vx = vel[i].x;
		real vy = vel[i].y;
		real vz = vel[i].z;
		real m  = mass[i];
		
		real r = mw_sqrt_0(sqr_0(x) + sqr_0(y) + sqr_0(z));
		real v2 = sqr_0(vx) + sqr_0(vy) + sqr_0(vz);
		
		T += m * v2;

		U1 += m * doubleCompPotential(comp1, comp2, r);

		for(unsigned int j = i + 1; j < numBodies; j++)
		{
			real rij = mw_sqrt_0(sqr_0(x - pos[j].x) + sqr_0(y - pos[j].y) + sqr_0(z - pos[j].z));
			U2 += 2.0 * m * mass[j] / rij;
		}
	}
	T  *= 0.5;
	U1 *= 0.5;
	U2 *= 0.5;

	//Compute the ratios. These should be 1 if they are in perfect virial equilibrium
	real ratio_1 = U1 / T / 2.0;
	real ratio_2 = U2 / T / 2.0;

	if (mw_fabs_0(1.0 - ratio_1) > TOL)
	{
		failed += 1;
		printf("\t Virial test 1 failed, unstable structure- T: %1f U1: %1f ratio_1: %1f\n", T, U1, ratio_1);
	}
	if (mw_fabs_0(1.0 - ratio_2) > TOL)
	{
		failed += 1;
		printf("\tVirial test 2 failed, gravitationally unstable- T: %1f U2: %1f ratio_2: %1f\n", T, U2, ratio_2);
	}
	
	//printf("T: %1f, U1: %1f U2: %1f ratio_1: %1f ratio_2: %1f\n", T, U1, U2, ratio_1, ratio_2);

	return failed;

}

//This function checks to ensure that the 2 components are individually and the dwarf as a whole are centered at 0
int checkCM(const Dwarf* comp1, const Dwarf* comp2, const mwvector* pos, const mwvector* vel, real* mass, unsigned int numBodies)
{
	//Note, this function relies on the fact that half the bodies are baryonic and half dark matter
	int failed = 0;
	const real TOL = 0.000001;

	real totalMass_l = comp1->mass;
	real totalMass_d = comp2->mass;
	real totalMass = totalMass_l + totalMass_d;
	
	real cm_x_comp1  = 0.0;
    real cm_y_comp1  = 0.0;
    real cm_z_comp1  = 0.0;
    real cm_vx_comp1 = 0.0;
    real cm_vy_comp1 = 0.0;
    real cm_vz_comp1 = 0.0;

	real cm_x_comp2  = 0.0;
    real cm_y_comp2  = 0.0;
    real cm_z_comp2  = 0.0;
    real cm_vx_comp2 = 0.0;
    real cm_vy_comp2 = 0.0;
    real cm_vz_comp2 = 0.0;

	for(unsigned int i = 0; i < numBodies / 2; i++)
	{
		cm_x_comp1 += mass[i] * pos[i].x;
		cm_y_comp1 += mass[i] * pos[i].y;
		cm_z_comp1 += mass[i] * pos[i].z;

		cm_vx_comp1 += mass[i] * vel[i].x;
		cm_vy_comp1 += mass[i] * vel[i].y;
		cm_vy_comp1 += mass[i] * vel[i].y;
	}

	for(unsigned int i = numBodies / 2; i < numBodies; i++)
	{
		cm_x_comp2 += mass[i] * pos[i].x;
		cm_y_comp2 += mass[i] * pos[i].y;
		cm_z_comp2 += mass[i] * pos[i].z;

		cm_vx_comp2 += mass[i] * vel[i].x;
		cm_vy_comp2 += mass[i] * vel[i].y;
		cm_vy_comp2 += mass[i] * vel[i].y;
	}

	cm_x_comp1 /= totalMass_l;
	cm_y_comp1 /= totalMass_l;
	cm_z_comp1 /= totalMass_l;
	
	cm_vx_comp1 /= totalMass_l;
	cm_vy_comp1 /= totalMass_l;
	cm_vz_comp1 /= totalMass_l;
	
	if(mw_fabs_0(cm_x_comp1) > TOL || mw_fabs_0(cm_y_comp1) > TOL || mw_fabs_0(cm_z_comp1) > TOL
		|| mw_fabs_0(cm_vx_comp1) > TOL || mw_fabs_0(cm_vy_comp1) > TOL || mw_fabs_0(cm_vz_comp1) > TOL)
	{
		failed += 1;
		printf("\tFailed: Component 1 is off center, cm_x: %1f cm_y: %1f cm_z: %1f cm_vx: %1f cm_vy: %1f cm_z: %1f\n",
		cm_x_comp1, cm_y_comp1, cm_z_comp1, cm_vx_comp1, cm_vy_comp1, cm_vz_comp1);
	}

	cm_x_comp2 /= totalMass_d;
	cm_y_comp2 /= totalMass_d;
	cm_z_comp2 /= totalMass_d;
	
	cm_vx_comp2 /= totalMass_d;
	cm_vy_comp2 /= totalMass_d;
	cm_vz_comp2 /= totalMass_d;
	
	if(mw_fabs_0(cm_x_comp2) > TOL || mw_fabs_0(cm_y_comp2) > TOL || mw_fabs_0(cm_z_comp2) > TOL
		|| mw_fabs_0(cm_vx_comp2) > TOL || mw_fabs_0(cm_vy_comp2) > TOL || mw_fabs_0(cm_vz_comp2) > TOL)
	{
		failed += 1;
		printf("\tFailed: Component 2 is off center, cm_x: %1f cm_y: %1f cm_z: %1f cm_vx: %1f cm_vy: %1f cm_z: %1f\n",
		cm_x_comp2, cm_y_comp2, cm_z_comp2, cm_vx_comp2, cm_vy_comp2, cm_vz_comp2);
	}

	real cm_x_total = (cm_x_comp1 + cm_x_comp2) / totalMass;
	real cm_y_total = (cm_y_comp1 + cm_y_comp2) / totalMass;
	real cm_z_total = (cm_z_comp1 + cm_z_comp2) / totalMass;
	
	real cm_vx_total = (cm_vx_comp1 + cm_vx_comp2) / totalMass;
	real cm_vy_total = (cm_vy_comp1 + cm_vy_comp2) / totalMass;
	real cm_vz_total = (cm_vz_comp1 + cm_vz_comp2) / totalMass;

	if(mw_fabs_0(cm_x_total) > TOL || mw_fabs_0(cm_y_total) > TOL || mw_fabs_0(cm_z_total) > TOL
		|| mw_fabs_0(cm_vx_total) > TOL || mw_fabs_0(cm_vy_total) > TOL || mw_fabs_0(cm_vz_total) > TOL)
	{
		failed += 1;
		printf("\tFailed: Total Dwarf is off center, cm_x: %1f cm_y: %1f cm_z: %1f cm_vx: %1f cm_vy: %1f cm_z: %1f\n",
		cm_x_total, cm_y_total, cm_z_total, cm_vx_total, cm_vy_total, cm_vz_total);
	}

	return failed;
}

//make a two component plummer-plummer dwarf and check it
int testPlummerPlummer()
{
	int failed = 0;
	
	unsigned int numBodies = 10000;
	mwvector* positions    = mwCalloc(numBodies, sizeof(mwvector));
	mwvector* velocities   = mwCalloc(numBodies, sizeof(mwvector));
	real* masses           = mwCalloc(numBodies, sizeof(real));
	
	//we want the dwarf to be at the origin
	mwvector rshift, vshift;
	rshift.x = 0;
	rshift.y = 0;
	rshift.z = 0;
	vshift.x = 0;
	vshift.y = 0;
	vshift.z = 0;

        Dwarf* comp1       = mwMalloc(sizeof(Dwarf));
	comp1->type        = Plummer;
	comp1->mass        = 12.0;
	comp1->scaleLength = 0.2;

	Dwarf* comp2       = mwMalloc(sizeof(Dwarf));
	comp2->type        = comp1->type;
	comp2->mass        = 24.0;
	comp2->scaleLength = 0.4;

	dsfmt_t prng;
	dsfmt_init_gen_rand(&prng, 1234); //initialize the random variable

	//Actually generate the dwarf bodies by calling a special version of the actual generation function from nbody_mixeddwarf.c
	nbGenerateMixedDwarfCore_TESTVER(positions, velocities, masses, &prng, numBodies, comp1, comp2, &rshift, &vshift);
	//printf("x: %1f y: %1f z: %1f vx: %1f vy: %1f vz: %1f\n", positions[0].x, positions[0].y, positions[0].z, velocities[0].x, velocities[0].y, velocities[0].z);
	
	printf("Checking Virial stability of plummer-plummer\n");
	failed += checkVirialRatio(comp1, comp2, positions, velocities, masses, numBodies);

	printf("Checking center of mass and momentum of plummer-plummer\n");
	failed += checkCM(comp1, comp2, positions, velocities, masses, numBodies);

	free(positions);
	free(velocities);
	free(masses);
	free(comp1);
	free(comp2);

	return failed;
}

//make a plummer-NFW dwarf and check it
int testPlummerNFW()
{
	int failed = 0;
	
	unsigned int numBodies = 10000; //The more bodies, the better the virial ratio will be. This is a good number of bodies
	mwvector* positions    = mwCalloc(numBodies, sizeof(mwvector));
	mwvector* velocities   = mwCalloc(numBodies, sizeof(mwvector));
	real* masses           = mwCalloc(numBodies, sizeof(real));
	
	//we want the dwarf to be at the origin
	mwvector rshift, vshift;
	rshift.x = 0;
	rshift.y = 0;
	rshift.z = 0;
	vshift.x = 0;
	vshift.y = 0;
	vshift.z = 0;

	//The Plummer component
        Dwarf* comp1       = mwMalloc(sizeof(Dwarf));
        comp1->type        = Plummer;
        comp1->mass        = 12.0;
        comp1->scaleLength = 0.2;
	
	//The NFW component
	Dwarf* comp2       = mwMalloc(sizeof(Dwarf));
	comp2->type        = NFW;
	comp2->mass        = 48.0;
	comp2->scaleLength = 0.8;

	dsfmt_t prng;
	dsfmt_init_gen_rand(&prng, 1234); //initialize the random variable

	//Actually generate the dwarf bodies by calling a special version of the actual generation function from nbody_mixeddwarf.c
	nbGenerateMixedDwarfCore_TESTVER(positions, velocities, masses, &prng, numBodies, comp1, comp2, &rshift, &vshift);
	//printf("x: %1f y: %1f z: %1f vx: %1f vy: %1f vz: %1f\n", positions[0].x, positions[0].y, positions[0].z, velocities[0].x, velocities[0].y, velocities[0].z);
	
	printf("Checking Virial stability of plummer-NFW\n");
	failed += checkVirialRatio(comp1, comp2, positions, velocities, masses, numBodies);

	printf("Checking center of mass and momentum of plummer-NFW\n");
	failed += checkCM(comp1, comp2, positions, velocities, masses, numBodies);

	free(positions);
	free(velocities);
	free(masses);
	free(comp1);
	free(comp2);

	return failed;
}

//make a NFW-NFW dwarf and check it
int testNFWNFW()
{
	int failed = 0;
	
	unsigned int numBodies = 40000; //The more bodies, the better the virial ratio will be. This is test in particular requires a lot of bodies to fall within our virial threshold
	mwvector* positions    = mwCalloc(numBodies, sizeof(mwvector));
	mwvector* velocities   = mwCalloc(numBodies, sizeof(mwvector));
	real* masses           = mwCalloc(numBodies, sizeof(real));
	
	//we want the dwarf to be at the origin
	mwvector rshift, vshift;
	rshift.x = 0;
	rshift.y = 0;
	rshift.z = 0;
	vshift.x = 0;
	vshift.y = 0;
	vshift.z = 0;

	//The Plummer component
        Dwarf* comp1       = mwMalloc(sizeof(Dwarf));
        comp1->type        = NFW;
        comp1->mass        = 12.0;
        comp1->scaleLength = 0.2;
	
	//The NFW component
	Dwarf* comp2       = mwMalloc(sizeof(Dwarf));
	comp2->type        = NFW;
	comp2->mass        = 48.0;
	comp2->scaleLength = 0.8;

	dsfmt_t prng;
	dsfmt_init_gen_rand(&prng, 1234); //initialize the random variable

	//Actually generate the dwarf bodies by calling a special version of the actual generation function from nbody_mixeddwarf.c
	nbGenerateMixedDwarfCore_TESTVER(positions, velocities, masses, &prng, numBodies, comp1, comp2, &rshift, &vshift);
	//printf("x: %1f y: %1f z: %1f vx: %1f vy: %1f vz: %1f\n", positions[0].x, positions[0].y, positions[0].z, velocities[0].x, velocities[0].y, velocities[0].z);
	
	printf("Checking Virial stability of NFW-NFW\n");
	failed += checkVirialRatio(comp1, comp2, positions, velocities, masses, numBodies);

	printf("Checking center of mass and momentum of NFW-NFW\n");
	failed += checkCM(comp1, comp2, positions, velocities, masses, numBodies);

	free(positions);
	free(velocities);
	free(masses);
	free(comp1);
	free(comp2);

	return failed;
}


int main()
{

    int failed = 0;

	failed += testPlummerPlummer();
	failed += testPlummerNFW();
	failed += testNFWNFW();

	if(failed == 0)
	{
		printf("All tests successful!\n");
	}

    return failed;
}
