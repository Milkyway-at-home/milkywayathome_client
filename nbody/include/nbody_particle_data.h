#ifndef NBODY_PARTICLE_DATA_H
#define NBODY_PARTICLE_DATA_H

#include "nbody_types.h"  // For real type definition

// Structure to hold header information
typedef struct {
    int simple_output;  
    int has_milkyway;     
    real com_x, com_y, com_z; 
    real cmom_x, cmom_y, cmom_z;
} ParticleHeader;

// Structure to hold data for a single particle (one row)
typedef struct {
    int  type;       // Column 0
    int  id;         // Column 1 
    real x;          // Column 2 
    real y;          // Column 3 
    real z;          // Column 4 
    real l;          // Column 5 
    real b;          // Column 6 
    real r;          // Column 7 
    real vx;         // Column 8 
    real vy;         // Column 9 
    real vz;         // Column 10 
    real mass;       // Column 11 
    real v_los;      // Column 12
    real pm_ra;      // Column 13
    real pm_dec;     // Column 14
    real lambda;     // Column 15
    real beta;       // Column 16
} ParticleData;

// Structure to manage a dynamic collection of particles AND header info
typedef struct {
    ParticleData *particles; // Pointer to the array of particles
    size_t count;            // Number of particles currently stored
    size_t capacity;         // Allocated capacity of the array
    ParticleHeader header;   // Header information
} ParticleCollection;

// Function declarations
ParticleCollection* create_particle_collection(size_t initial_capacity);
void free_particle_collection(ParticleCollection *collection);
int add_particle(ParticleCollection *collection, ParticleData particle);
ParticleCollection* read_particle_file(const char *filename);

#endif // NBODY_PARTICLE_DATA_H