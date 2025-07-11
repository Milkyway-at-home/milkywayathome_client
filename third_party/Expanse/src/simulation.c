// src/simulation.c
#include "simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bfe.h"

// Gravitational constant and a central mass at the origin
const double G = 6.67430e-11;
const double M_central = 1.989e30; // Mass of the Sun

// Creates a new SimulationState object
SimulationState* simulation_create(const char* filename) {
    SimulationState* state = (SimulationState*)malloc(sizeof(SimulationState));
    if (!state) {
        fprintf(stderr, "Failed to allocate memory for SimulationState\n");
        return NULL;
    }

    state->particle_count = 0;
    state->particle_capacity = MAX_PARTICLES;
    state->current_time = 0.0;
    state->model = bfe_create_from_file(filename);
    if (state->model == NULL) {
        fprintf(stderr, "CRITICAL: BFE model failed to load. Aborting simulation creation.\n");
        free(state); // Free the partially created state
        return NULL; // Return NULL to signal failure
    }

    int lmax = state->model->lmax;
    int legendre_size = (lmax + 1) * (lmax + 2) / 2;

    state->bfe_temp_cos_m_phi = malloc((lmax + 1) * sizeof(double));
    state->bfe_temp_sin_m_phi = malloc((lmax + 1) * sizeof(double));
    state->bfe_temp_Plm = malloc(legendre_size * sizeof(double));
    state->bfe_temp_dPlm_dtheta = malloc(legendre_size * sizeof(double));

    if (!state->bfe_temp_cos_m_phi || !state->bfe_temp_sin_m_phi || !state->bfe_temp_Plm || !state->bfe_temp_dPlm_dtheta) {
        // Handle memory error, destroy everything, and return NULL
        simulation_destroy(state); // A single cleanup function is useful here
        return NULL;
    }

    state->particles = (Particle*)malloc(state->particle_capacity * sizeof(Particle));
    if (!state->particles) {
        fprintf(stderr, "Failed to allocate memory for particles\n");
        bfe_destroy(state->model);
        free(state);
        return NULL;
    }

    return state;
}

void simulation_destroy(SimulationState* state) {
    if (state) {
        bfe_destroy(state->model);
        if (state->particles) {
            free(state->particles);
        }
        free(state->bfe_temp_cos_m_phi);
        free(state->bfe_temp_sin_m_phi);
        free(state->bfe_temp_Plm);
        free(state->bfe_temp_dPlm_dtheta);

        free(state);
    }
}

// Calculates forces and updates particle accelerations
void apply_forces(SimulationState* state) {
    for (int i = 0; i < state->particle_count; i++) {
        Particle* p = &state->particles[i];

        double f_vec[3];
        // Pass the pre-allocated arrays from the state
        bfe_calculate_force(p->position, state->model, f_vec,
                            state->bfe_temp_cos_m_phi,
                            state->bfe_temp_sin_m_phi,
                            state->bfe_temp_Plm,
                            state->bfe_temp_dPlm_dtheta);

        p->acceleration[0] = f_vec[0] / p->mass;
        p->acceleration[1] = f_vec[1] / p->mass;
        p->acceleration[2] = f_vec[2] / p->mass;
    }
}

// Updates positions and velocities using the Velocity Verlet algorithm
void update_particles(SimulationState* state, double dt) {
    for (int i = 0; i < state->particle_count; i++) {
        Particle* p = &state->particles[i];

        // 1. Update position: p_new = p_old + v * dt + 0.5 * a * dt^2
        p->position[0] += p->velocity[0] * dt + 0.5 * p->acceleration[0] * dt * dt;
        p->position[1] += p->velocity[1] * dt + 0.5 * p->acceleration[1] * dt * dt;
        p->position[2] += p->velocity[2] * dt + 0.5 * p->acceleration[2] * dt * dt;

        // Store old acceleration before recalculating
        double old_accel[3];
        memcpy(old_accel, p->acceleration, sizeof(old_accel));

        // We will recalculate forces/accelerations at the new position
        // For now, we assume it's done externally. We will do this in simulation_step.
        // (A placeholder for the force update)

        // 2. Update velocity: v_new = v_old + 0.5 * (a_old + a_new) * dt
        // We will complete this in the main step function.
    }
}

// The main simulation step using the Velocity Verlet algorithm
void simulation_step(SimulationState* state, double dt) {
    // First half-step velocity update
    for (int i = 0; i < state->particle_count; i++) {
        Particle* p = &state->particles[i];
        p->velocity[0] += 0.5 * p->acceleration[0] * dt;
        p->velocity[1] += 0.5 * p->acceleration[1] * dt;
        p->velocity[2] += 0.5 * p->acceleration[2] * dt;
    }

    // Update positions to full step
    for (int i = 0; i < state->particle_count; i++) {
        Particle* p = &state->particles[i];
        p->position[0] += p->velocity[0] * dt;
        p->position[1] += p->velocity[1] * dt;
        p->position[2] += p->velocity[2] * dt;
    }

    // Calculate new forces/accelerations at the new positions
    apply_forces(state);

    // Second half-step velocity update
    for (int i = 0; i < state->particle_count; i++) {
        Particle* p = &state->particles[i];
        p->velocity[0] += 0.5 * p->acceleration[0] * dt;
        p->velocity[1] += 0.5 * p->acceleration[1] * dt;
        p->velocity[2] += 0.5 * p->acceleration[2] * dt;
    }

    // Update BFE coefficients
    bfe_evolve_coeffs(state->model, state->particles, state->particle_count, dt);

    // Update current simulation time
    state->current_time += dt;
}

// Loads the simulation state from a file
int simulation_load_particles(SimulationState* state, const char* filename) {
    if (!state) return 0;

    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("[DEBUG] C: fopen failed");
        state->particle_count = 0;
        return 0;
    }

    state->particle_count = 0;
    double m, px, py, pz, vx, vy, vz;

    int items_read = fscanf(file, "# Mass (in 1e10 M_sun), Pos X, Y, Z (kpc), Vel X, Y, Z (km/s)\n%lf %lf %lf %lf %lf %lf %lf", &m, &px, &py, &pz, &vx, &vy, &vz);

    while (items_read == 7) {
        if (state->particle_count < MAX_PARTICLES) {
            int i = state->particle_count;
            state->particles[i] = (Particle){
                .id = i + 1, 
                .mass = m, 
                .position = {px, py, pz}, 
                .velocity = {vx, vy, vz},
                .acceleration = {0.0, 0.0, 0.0} // Also initialize acceleration
            };
            state->particle_count++;
        }
        items_read = fscanf(file, "%lf %lf %lf %lf %lf %lf %lf", &m, &px, &py, &pz, &vx, &vy, &vz);
    } 
    fclose(file);
    return state->particle_count;
}
