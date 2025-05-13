/* This program tests the generation of a 2-component dwarf galaxy.
* It creates a plummer-plummer dwarf and an NFW dwarf.
* The structure and stability of these two dwarfs are checked by calculating their virial ratios and checking their center of masses.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include "dSFMT.h"
#include "nbody_mixeddwarf.h"
#include "nbody_particle_data.h"
#include "test_env_util.h"

/* Dwarf galaxy parameters */
#define EVOLUTION_TIME "2.0"                  /* Evolution time in Gyr */
#define EVOLUTION_RATIO "0.0"                 /* Evolution time ratio */
#define BARYON_SCALE_RADIUS "0.181216"        /* Baryon Scale radius in kpc */
#define SCALE_RADIUS_RATIO "0.182799"         /* Scale radius ratio */
#define BARYON_MASS "1.22251"                 /* Baryon Mass in SMU */
#define MASS_RATIO "0.0126171"  

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
int checkCM(const Dwarf* comp1, const Dwarf* comp2, const mwvector* pos, const mwvector* vel, real* mass, unsigned int numBodies, unsigned int numBodies_light)
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

	for(unsigned int i = 0; i < numBodies_light; i++)
	{
		cm_x_comp1 += mass[i] * pos[i].x;
		cm_y_comp1 += mass[i] * pos[i].y;
		cm_z_comp1 += mass[i] * pos[i].z;

		cm_vx_comp1 += mass[i] * vel[i].x;
		cm_vy_comp1 += mass[i] * vel[i].y;
		cm_vy_comp1 += mass[i] * vel[i].y;
	}

	for(unsigned int i = numBodies_light; i < numBodies; i++)
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

// This function tests the nbGenerateMixedDwarf function (virial test and center of mass test) for a given dwarf potential type set in test lua files 
int testMixedDwarf(const char* dwarf_potential_type_lua)
{
    printf("Starting mixed dwarf test for %s\n", dwarf_potential_type_lua);
    fflush(stdout);

    int failed = 0;
    char* input_lua_file = NULL;
    
    // Variables to store the parameters from Lua
    real nbody = 0.0;              
    real nbody_baryon = 0.0;      
    Dwarf* comp1 = NULL;           
    Dwarf* comp2 = NULL;      
    real timestep = 0.0;      
    lua_State* lua_state = NULL;

    // Initialize particle arrays
    ParticleCollection* particle_data = NULL;
    mwvector* positions = NULL;
    mwvector* velocities = NULL;
    real* masses = NULL;

    const char* dwarf_params[] = {
        EVOLUTION_TIME,
        EVOLUTION_RATIO,
        BARYON_SCALE_RADIUS,
        SCALE_RADIUS_RATIO,
        BARYON_MASS,
        MASS_RATIO
    };

     // Clean up any output files from previous runs
    printf("Cleaning up any output files from previous runs...\n");
    fflush(stdout);

    // Remove initial output file 
    const char* initial_output_file = "initial.out";
    if (access(initial_output_file, F_OK) != -1) {
        if (remove(initial_output_file) == 0) {
            printf("Removed old output file: %s\n", initial_output_file);
        } else {
            fprintf(stderr, "Failed to remove old output file: %s - %s\n", 
                  initial_output_file, strerror(errno));
        }
    }

    // Find the lua file
    input_lua_file = find_lua_file(dwarf_potential_type_lua);
    if (!input_lua_file) {
        printf("Error: Could not find %s in any expected location\n", dwarf_potential_type_lua);
        failed = 1;
        return failed;
    }

    printf("Input file path: %s\n", input_lua_file);

    // Read the parameters from the Lua file
    if (read_lua_parameters(input_lua_file, dwarf_params, &nbody, &nbody_baryon, &comp1, &comp2, &timestep, &lua_state) != 0) {
        fprintf(stderr, "Error: Failed to read Lua parameters\n");
        fflush(stdout);
        failed = 1;
        return failed;
    }

    printf("Lua state: %p\n", lua_state);
    fflush(stdout);

    // Generate the mixed dwarf from the lua state
    registerGenerateMixedDwarf(lua_state);

    // Check if the initial output file exists
    const char* initial_output_filename = "initial.out";
    printf("Checking for initial output file: %s\n", initial_output_filename);
    fflush(stdout);

    // Read the initial output file
    particle_data = read_particle_file(initial_output_filename);
    if (!particle_data) {
        fprintf(stderr, "Error: Failed to read initial output file '%s'\n", initial_output_filename);
        fflush(stdout);
        failed = 1;
        return failed;
    }

    printf("Successfully read %zu particles from initial file\n", particle_data->count);
    fflush(stdout);

    // Create arrays for positions, velocities, and masses
    positions = (mwvector*)malloc(nbody * sizeof(mwvector));
    velocities = (mwvector*)malloc(nbody * sizeof(mwvector));
    masses = (real*)malloc(nbody * sizeof(real));

    for (unsigned int i = 0; i < nbody; i++) {
        positions[i] = (mwvector){particle_data->particles[i].x, particle_data->particles[i].y, particle_data->particles[i].z};
        velocities[i] = (mwvector){particle_data->particles[i].vx, particle_data->particles[i].vy, particle_data->particles[i].vz};
        masses[i] = particle_data->particles[i].mass;
    }

    fflush(stdout);

    // Center of mass test
    if (checkCM(comp1, comp2, positions, velocities, masses, nbody, nbody_baryon) != 0) {
        fprintf(stderr, "Error: Center of mass test failed\n");
        fflush(stdout);
        failed = 1;
        return failed;
    }

    printf("Center of mass test passed\n");

    // Virial test
    if (checkVirialRatio(comp1, comp2, positions, velocities, masses, nbody) != 0) {
        fprintf(stderr, "Error: Virial test failed\n");
        fflush(stdout);
        failed = 1;
        return failed;
    }
    
    printf("Virial test passed\n");

    // Free memory
    free_particle_collection(particle_data);
    free(input_lua_file);
    free(positions);
    free(velocities);
    free(masses);

    // Close the Lua state
    lua_close(lua_state);
    
    return failed;
}

int main()
{
    int failed = 0;
    int total_failed = 0;
    const char* dwarf_models[] = {
        "plummer_plummer.lua",
        "plummer_nfw.lua",
        //"cored_cored.lua" 
    };
    
    int num_models = sizeof(dwarf_models) / sizeof(dwarf_models[0]);
    
    for (int i = 0; i < num_models; i++) {
        printf("\n=== Testing %s ===\n", dwarf_models[i]);
        failed = testMixedDwarf(dwarf_models[i]);
        if (failed == 0) {
            printf("✓ %s passed all tests!\n", dwarf_models[i]);
        } else {
            printf("✗ %s failed %d tests!\n", dwarf_models[i], failed);
            total_failed += failed;
        }
    }

	if(total_failed == 0)
	{
		printf("\n=== SUMMARY ===\n");
		printf("All center of mass and virial stability tests successful!\n");
	}
	else
	{
		printf("\n=== SUMMARY ===\n");
		printf("Failed %d tests across all models\n", total_failed);
	}

    return total_failed > 0 ? 1 : 0;
}
