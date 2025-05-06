/* This program tests the generation of a 2-component dwarf galaxy.
* It creates a plummer-plummer dwarf and an NFW dwarf.
* The structure and stability of these two dwarfs are checked by calculating their virial ratios and checking their center of masses.
*/

#include "nbody_mixeddwarf.h"
#include "nbody_dwarf_potential.h"
#include "dSFMT.h"
#include <stdio.h>
#include <errno.h>
#include "test_env_util.h"
#include "nbody_lua.h"
#include "nbody_io.h"
#include "nbody_types.h"
#include "nbody_particle_data.h"
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

/* Dwarf galaxy parameters */
#define EVOLUTION_TIME "2.0"                  /* Evolution time in Gyr */
#define EVOLUTION_RATIO "0.0"                 /* Evolution time ratio */
#define BARYON_SCALE_RADIUS "0.181216"        /* Baryon Scale radius in kpc */
#define SCALE_RADIUS_RATIO "0.182799"         /* Scale radius ratio */
#define BARYON_MASS "1.22251"                 /* Baryon Mass in SMU */
#define MASS_RATIO "0.0126171"                /* Mass ratio */

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
		
		real r = mw_sqrt(sqr(x) + sqr(y) + sqr(z));
		real v = sqr(vx) + sqr(vy) + sqr(vz);
		
		T += m * v;

		U1 += m * doubleCompPotential(comp1, comp2, r);

		for(unsigned int j = i + 1; j < numBodies; j++)
		{
			real rij = mw_sqrt(sqr(x - pos[j].x) + sqr(y - pos[j].y) + sqr(z - pos[j].z));
			U2 += 2.0 * m * mass[j] / rij;
		}
	}
	T  *= 0.5;
	U1 *= 0.5;
	U2 *= 0.5;

	//Compute the ratios. These should be 1 if they are in perfect virial equilibrium
	real ratio_1 = U1 / T / 2.0;
	real ratio_2 = U2 / T / 2.0;

	if (mw_fabs(1.0 - ratio_1) > TOL)
	{
		failed += 1;
		printf("\t Virial test 1 failed, unstable structure- T: %1f U1: %1f ratio_1: %1f\n", T, U1, ratio_1);
	}
	if (mw_fabs(1.0 - ratio_2) > TOL)
	{
		failed += 1;
		printf("\tVirial test 2 failed, gravitationally unstable- T: %1f U2: %1f ratio_2: %1f\n", T, U2, ratio_2);
	}
	
	//printf("T: %1f, U1: %1f U2: %1f ratio_1: %1f ratio_2: %1f\n", T, U1, U2, ratio_1, ratio_2);

	return failed;

}

//This function checks to ensure that the 2 components are individually and the dwarf as a whole are centered at 0
int checkCM(const Dwarf* comp1, const Dwarf* comp2, const mwvector* pos, const mwvector* vel, real* mass, unsigned int numBodies, unsigned int numBodies_baryon)
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

	for(unsigned int i = 0; i < numBodies_baryon; i++)
	{
		cm_x_comp1 += mass[i] * pos[i].x;
		cm_y_comp1 += mass[i] * pos[i].y;
		cm_z_comp1 += mass[i] * pos[i].z;

		cm_vx_comp1 += mass[i] * vel[i].x;
		cm_vy_comp1 += mass[i] * vel[i].y;
		cm_vy_comp1 += mass[i] * vel[i].y;
	}

	for(unsigned int i = numBodies_baryon; i < numBodies; i++)
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
	
	if(mw_fabs(cm_x_comp1) > TOL || mw_fabs(cm_y_comp1) > TOL || mw_fabs(cm_z_comp1) > TOL
		|| mw_fabs(cm_vx_comp1) > TOL || mw_fabs(cm_vy_comp1) > TOL || mw_fabs(cm_vz_comp1) > TOL)
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
	
	if(mw_fabs(cm_x_comp2) > TOL || mw_fabs(cm_y_comp2) > TOL || mw_fabs(cm_z_comp2) > TOL
		|| mw_fabs(cm_vx_comp2) > TOL || mw_fabs(cm_vy_comp2) > TOL || mw_fabs(cm_vz_comp2) > TOL)
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

	if(mw_fabs(cm_x_total) > TOL || mw_fabs(cm_y_total) > TOL || mw_fabs(cm_z_total) > TOL
		|| mw_fabs(cm_vx_total) > TOL || mw_fabs(cm_vy_total) > TOL || mw_fabs(cm_vz_total) > TOL)
	{
		failed += 1;
		printf("\tFailed: Total Dwarf is off center, cm_x: %1f cm_y: %1f cm_z: %1f cm_vx: %1f cm_vy: %1f cm_z: %1f\n",
		cm_x_total, cm_y_total, cm_z_total, cm_vx_total, cm_vy_total, cm_vz_total);
	}

	return failed;
}

// make a mixed dwarf with certain potential and check it
int testMixedDwarf(const char* dwarf_potential_type)
{
	int failed = 0;
	char* input_lua_file = NULL;
	real nbody = 0.0;
    real nbody_baryon = 0.0;
    real nbody_dark = 0.0;
    Dwarf *comp1 = NULL;
    Dwarf *comp2 = NULL;
    real timestep = 0.0;

	const char* dwarf_params[] = {
        EVOLUTION_TIME,
        EVOLUTION_RATIO,
        BARYON_SCALE_RADIUS,
        SCALE_RADIUS_RATIO,
        BARYON_MASS,
        MASS_RATIO
    };

	// Clean up any output files from previous runs
    printf("Cleaning up output file from previous runs...\n");
    fflush(stdout);

    // Remove initial output file 
    const char* initial_output_files[] = {"initial.out", "initial.hist"};
    for (int i = 0; i < sizeof(initial_output_files) / sizeof(initial_output_files[0]); i++) {
        if (access(initial_output_files[i], F_OK) != -1) {
            if (remove(initial_output_files[i]) == 0) {
                printf("Removed old output file: %s\n", initial_output_files[i]);
            } else {
                fprintf(stderr, "Failed to remove old output file: %s - %s\n", 
                      initial_output_files[i], strerror(errno));
            }
        }
    }
    
    // Find the plummer_plummer.lua file
    input_lua_file = find_lua_file(dwarf_potential_type);
    if (!input_lua_file) {
        printf("Error: Could not find %s.lua in any expected location\n", dwarf_potential_type);
        failed = 1;
        return failed;
    }

	printf("Input file path: %s\n", input_lua_file);

	// First, read the parameters from the Lua file
    printf("Reading Lua parameters...\n");
	if (read_lua_parameters(input_lua_file, dwarf_params, &nbody, &nbody_baryon, &comp1, &comp2, &timestep) != 0) {
        printf("Error: Failed to read Lua parameters\n");
        failed = 1;
        free(input_lua_file);
        return failed;
    }
    
    unsigned int numBodies = (unsigned int)nbody;
    unsigned int numBodies_baryon = (unsigned int)nbody_baryon;
    
    // Initialize Lua state
	NBodyFlags nbf = {
        .inputFile = input_lua_file,
        .debugLuaLibs = TRUE,  // Enable debug output
        .forwardedArgs = dwarf_params,
        .numForwardedArgs = 6
    };
    
    lua_State* luaSt = nbOpenLuaStateWithScript(&nbf, NULL);
    if (!luaSt) {
        printf("Error: Failed to open Lua state\n");
        failed = 1;
        free(input_lua_file);
        return failed;
    }

    // Generate the mixed dwarf
    printf("Generating mixed dwarf galaxy...\n");
    
    // Instead of trying to create a PRNG directly, let Lua do the work
    // Execute Lua code to create a PRNG and store it in a global variable
    if (luaL_dostring(luaSt, "prng = mw.createRandom(1234)") != 0) {
        printf("Error creating PRNG in Lua: %s\n", lua_tostring(luaSt, -1));
        failed = 1;
        lua_close(luaSt);
        free(input_lua_file);
        return failed;
    }
    
    // Push a table onto the Lua stack with all the parameters needed for mixeddwarf generation
    lua_newtable(luaSt);
    
    // Add nbody parameter
    lua_pushnumber(luaSt, nbody);
    lua_setfield(luaSt, -2, "nbody");
    
    // Add nbody_baryon parameter
    lua_pushnumber(luaSt, nbody_baryon);
    lua_setfield(luaSt, -2, "nbody_baryon");
    
    // Add comp1 and comp2
    lua_pushlightuserdata(luaSt, comp1);
    lua_setfield(luaSt, -2, "comp1");
    
    lua_pushlightuserdata(luaSt, comp2);
    lua_setfield(luaSt, -2, "comp2");
    
    // Create and add position vector (at origin)
    mwvector* position = mwCalloc(1, sizeof(mwvector));
    position->x = 0.0;
    position->y = 0.0;
    position->z = 0.0;
    lua_pushlightuserdata(luaSt, position);
    lua_setfield(luaSt, -2, "position");
    
    // Create and add velocity vector (zero velocity)
    mwvector* velocity = mwCalloc(1, sizeof(mwvector));
    velocity->x = 0.0;
    velocity->y = 0.0;
    velocity->z = 0.0;
    lua_pushlightuserdata(luaSt, velocity);
    lua_setfield(luaSt, -2, "velocity");
    
    // Add ignore flag
    lua_pushboolean(luaSt, 0);  // false
    lua_setfield(luaSt, -2, "ignore");
    
    // Get the PRNG from the global variable we just created
    lua_getglobal(luaSt, "prng");
    lua_setfield(luaSt, -2, "prng");
    
    // Call the function with the table on the stack
    nbGenerateMixedDwarf(luaSt);
    
	printf("Mixed dwarf galaxy generated\n");
    
    // Read the generated particle data
    ParticleCollection* particles = read_particle_file("initial.out");
    if (!particles) {
        printf("Error: Failed to read initial.out file\n");
        failed = 1;
        free(position);
        free(velocity);
        lua_close(luaSt);
        free(input_lua_file);
        free(comp1);
        free(comp2);
        return failed;
    }
    
    // Verify correct number of particles
    if (particles->count != numBodies) {
        printf("Error: Expected %u particles, but got %zu\n", numBodies, particles->count);
        failed = 1;
        free_particle_collection(particles);
        free(position);
        free(velocity);
        lua_close(luaSt);
        free(input_lua_file);
        free(comp1);
        free(comp2);
        return failed;
    }
    
    printf("Successfully read %zu particles\n", particles->count);
    
    // Create arrays for position, velocity, and mass data
    mwvector* positions = mwCalloc(numBodies, sizeof(mwvector));
    mwvector* velocities = mwCalloc(numBodies, sizeof(mwvector));
    real* masses = mwCalloc(numBodies, sizeof(real));
    
    if (!positions || !velocities || !masses) {
        printf("Error: Failed to allocate memory for body data\n");
        failed = 1;
        if (positions) free(positions);
        if (velocities) free(velocities);
        if (masses) free(masses);
        free_particle_collection(particles);
        free(position);
        free(velocity);
        lua_close(luaSt);
        free(input_lua_file);
        free(comp1);
        free(comp2);
        return failed;
    }
    
    // Copy data from particles to our arrays
    for (unsigned int i = 0; i < numBodies; i++) {
        positions[i].x = particles->particles[i].x;
        positions[i].y = particles->particles[i].y;
        positions[i].z = particles->particles[i].z;
        
        velocities[i].x = particles->particles[i].vx;
        velocities[i].y = particles->particles[i].vy;
        velocities[i].z = particles->particles[i].vz;
        
        masses[i] = particles->particles[i].mass;
    }

	printf("Checking Virial stability of %s\n", dwarf_potential_type);
	failed += checkVirialRatio(comp1, comp2, positions, velocities, masses, numBodies);

	printf("Checking center of mass and momentum of %s\n", dwarf_potential_type);
	failed += checkCM(comp1, comp2, positions, velocities, masses, numBodies, numBodies_baryon);

    // Clean up
    free_particle_collection(particles);
    free(positions);
    free(velocities);
    free(masses);
    free(position);
    free(velocity);
	free(input_lua_file);
	free(comp1);
	free(comp2);
    lua_close(luaSt);

	return failed;
}

int main()
{

    int failed = 0;

	failed += testMixedDwarf("plummer_plummer.lua");
	// failed += testMixedDwarf("plummer_nfw.lua");
	// failed += testMixedDwarf("cored_cored.lua");

	if(failed == 0)
	{
		printf("All center of mass and virial stability tests successful!\n");
	}

    return failed;
}
