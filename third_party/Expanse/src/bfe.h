// src/bfe.h
#ifndef BFE_H
#define BFE_H

#include "particle.h"

typedef struct {
    int nmax; // Max expansion order
    int lmax;
    double scale_radius; // A scaling factor 'a'
    double* S_coeffs;    // Pointer to sine coefficients
    double* T_coeffs;    // Pointer to cosine coefficients
} BFEModel;

// Create a BFEModel from a file
BFEModel* bfe_create_from_file(const char* filename);

// Destroys given BFEModel
void bfe_destroy(BFEModel* model);

// Calculates force at a given position using the BFE model
void bfe_calculate_force(double pos[3], const BFEModel* model, double* force_out, 
                         double* cos_m_phi, double* sin_m_phi, double* Plm, double* dPlm_dtheta);

// Evolves the BFE coefficients based on particle positions and velocities
void bfe_evolve_coeffs(BFEModel* model, const Particle particles[], int num_particles, double dt);

#endif // BFE_H