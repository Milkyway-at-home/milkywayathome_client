#ifndef EXPANSE_H
#define EXPANSE_H

#include <stddef.h> 
#include "nbody_types.h" // For mwvector and Body structs

// Define the coefficient structure that your functions use internally.
// This must be defined here so that nbody_grav.c knows its size and layout.
struct bfe_coeffs {
    int nmax;
    int lmax;
    double* Snlm; 
    double* Tnlm;
    // Add any other members your struct requires
};

// --- Function Prototypes ---
// Declarations for functions implemented in your various .c files.

/**
 * @brief Computes the Basis Function Expansion coefficients for a set of bodies.
 * Implemented in: bfe_coefficients.c (or similar)
 */
struct bfe_coeffs* bfe_compute_coefficients(const Body* bodies, size_t nbody, int nmax, int lmax);

/**
 * @brief Calculates the acceleration at a specific position using pre-computed BFE coefficients.
 * Implemented in: bfe_acceleration.c (or similar)
 */
mwvector bfe_get_acceleration(mwvector pos, const struct bfe_coeffs* coeffs);

/**
 * @brief Frees the memory allocated for a bfe_coeffs struct.
 * Implemented in: bfe_memory.c (or similar)
 */
void bfe_free_coefficients(struct bfe_coeffs* coeffs);

#endif // EXPANSE_H