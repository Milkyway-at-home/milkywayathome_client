#ifndef PARTICLE_H
#define PARTICLE_H

// Defines the fundamental object in our simulation.
// Using 'double' for high precision, which is standard in physics sims.
typedef struct {
    double mass;
    double position[3]; // [x, y, z]
    double velocity[3]; // [vx, vy, vz]
    double acceleration[3]; // [ax, ay, az]
    int id;
} Particle;

#endif // PARTICLE_H