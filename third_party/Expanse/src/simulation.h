#ifndef SIMULATION_H
#define SIMULATION_H

#include "particle.h"
#include "bfe.h"

#define MAX_PARTICLES 4e6 // A simple, fixed-size array for now

typedef struct {
    Particle* particles;
    int particle_count;
    int particle_capacity; 
    BFEModel* model;
    double current_time;

    double* bfe_temp_cos_m_phi;
    double* bfe_temp_sin_m_phi;
    double* bfe_temp_Plm;
    double* bfe_temp_dPlm_dtheta;
} SimulationState;

// Creates a new SimulationState object
SimulationState* simulation_create(const char* filename);

// Destroys the SimulationState object and frees resources
void simulation_destroy(SimulationState* state);

// This is the main function that will advance the simulation by one step
void simulation_step(SimulationState* state, double dt);

// Ths function initializes the simulation state
int simulation_load_particles(SimulationState* state, const char* filename);

// This function applies forces to all particles based on their current positions
void apply_forces(SimulationState* state);

#endif // SIMULATION_H