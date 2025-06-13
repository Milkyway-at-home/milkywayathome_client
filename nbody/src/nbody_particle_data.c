#include "nbody_particle_data.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

ParticleCollection* create_particle_collection(size_t initial_capacity) {
    ParticleCollection *collection = malloc(sizeof(ParticleCollection));
    if (!collection) {
        perror("Failed to allocate ParticleCollection");
        return NULL;
    }
    collection->capacity = (initial_capacity > 0) ? initial_capacity : 100;
    collection->count = 0;
    collection->particles = malloc(collection->capacity * sizeof(ParticleData));
    if (!collection->particles) {
        perror("Failed to allocate initial particle array");
        free(collection);
        return NULL;
    }

    // Initialize header fields
    collection->header.simple_output = 1;
    collection->header.has_milkyway = -1;
    collection->header.com_x = (real)0.0; 
    collection->header.com_y = (real)0.0; 
    collection->header.com_z = (real)0.0;
    collection->header.cmom_x = (real)0.0; 
    collection->header.cmom_y = (real)0.0; 
    collection->header.cmom_z = (real)0.0;

    return collection;
}

void free_particle_collection(ParticleCollection *collection) {
    if (!collection) return;
    free(collection->particles);
    free(collection);
}

int add_particle(ParticleCollection *collection, ParticleData particle) {
    if (collection->count >= collection->capacity) {
        // Resize needed
        size_t new_capacity = collection->capacity * 2;
        if (new_capacity <= collection->capacity) { // Handle potential overflow
            fprintf(stderr, "Error: Cannot resize particle collection further (capacity overflow).\n");
            return 0; // Indicate failure
        }
        ParticleData *new_particles = realloc(collection->particles, new_capacity * sizeof(ParticleData));
        if (!new_particles) {
            perror("Failed to reallocate particle array");
            return 0; // Indicate failure
        }
        collection->particles = new_particles;
        collection->capacity = new_capacity;
    }
    collection->particles[collection->count++] = particle;
    return 1; // Indicate success
}

ParticleCollection* read_particle_file(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }

    ParticleCollection *collection = create_particle_collection(20000);
    if (!collection) {
        fclose(file);
        return NULL;
    }

    char line_buffer[1024]; // Buffer for reading lines
    long line_number = 0;
    int header_lines_processed = 0; // Counter for processed header lines
    int in_header = 1; // Flag to indicate we are still processing header
    int simple_output = 1; // Default to simple output

    while (in_header && fgets(line_buffer, sizeof(line_buffer), file)) {
        line_number++;
        header_lines_processed++;

        // Trim leading/trailing whitespace
        char *start = line_buffer;
        while (*start == ' ' || *start == '\t') start++;
        char *end = start + strlen(start) - 1;
        while (end > start && (*end == ' ' || *end == '\t' || *end == '\n' || *end == '\r')) *end-- = '\0';

        // Parse header lines
        if (strstr(start, "simple_output") == start) {
            sscanf(start, "simple_output = %d", &simple_output);
            collection->header.simple_output = simple_output;
        } else if (strstr(start, "hasMilkyway") == start) {
            sscanf(start, "hasMilkyway = %d", &collection->header.has_milkyway);
        } else if (strstr(start, "centerOfMass") != NULL) {
            // Handle both combined and separate formats
            double com_x, com_y, com_z, cmom_x, cmom_y, cmom_z;
            char *com_equals = strstr(start, "centerOfMass =");
            char *cmom_equals = strstr(start, "centerOfMomentum =");
            
            if (com_equals) {
                com_equals += 13; // Skip "centerOfMass ="
                while (*com_equals == ' ' || *com_equals == '\t') com_equals++;
                // Skip the '=' character
                if (*com_equals == '=') {
                    com_equals++;
                    while (*com_equals == ' ' || *com_equals == '\t') com_equals++;
                }
                // Parse center of mass values
                if (sscanf(com_equals, "%lf, %lf, %lf", &com_x, &com_y, &com_z) == 3) {
                    collection->header.com_x = (real)com_x;
                    collection->header.com_y = (real)com_y;
                    collection->header.com_z = (real)com_z;
                }
            }
            
            if (cmom_equals) {
                cmom_equals += 16; // Skip "centerOfMomentum ="
                while (*cmom_equals == ' ' || *cmom_equals == '\t') cmom_equals++;
                // Skip any remaining characters from centerOfMass values
                while (*cmom_equals && *cmom_equals != '=') cmom_equals++;
                if (*cmom_equals == '=') {
                    cmom_equals++; // Skip the '='
                    while (*cmom_equals == ' ' || *cmom_equals == '\t') cmom_equals++;
                    if (sscanf(cmom_equals, "%lf, %lf, %lf", &cmom_x, &cmom_y, &cmom_z) == 3) {
                        collection->header.cmom_x = (real)cmom_x;
                        collection->header.cmom_y = (real)cmom_y;
                        collection->header.cmom_z = (real)cmom_z;
                    }
                }
            }
        } else if (start[0] == '#') {
            if (strstr(start, "ignore") && strstr(start, "id")) {
                in_header = 0;
            }
        } else {
            in_header = 0;
        }

        if (!in_header) {
            header_lines_processed--;
            break;
        }
    }

    if (feof(file) && in_header) {
        fclose(file);
        return collection;
    }

    int first_data_line_needs_processing = !in_header;

    do {
        if (!first_data_line_needs_processing) {
            if (fgets(line_buffer, sizeof(line_buffer), file) == NULL) {
                break;
            }
            line_number++;
        }
        first_data_line_needs_processing = 0;

        ParticleData current_particle = {0}; // Initialize all fields to 0
        char *token;
        char *rest = line_buffer;
        char *endptr;
        int col = 0;
        int parse_error = 0;
        long temp_long;
        double temp_double;
        int has_lambda_beta = 0; // Flag to track if lambda/beta columns are present

        while (*rest == ' ' || *rest == '\t') rest++;
        char *line_end = rest + strlen(rest) - 1;
        while (line_end > rest && (*line_end == ' ' || *line_end == '\t' || *line_end == '\n' || *line_end == '\r')) *line_end-- = '\0';

        // First pass: count columns to determine if lambda/beta are present
        char *count_rest = line_buffer;
        int total_cols = 0;
        while (strtok_r(count_rest, ",", &count_rest)) {
            total_cols++;
        }

        // Reset for actual parsing
        rest = line_buffer;
        while ((token = strtok_r(rest, ",", &rest)) != NULL) {
            while (*token == ' ' || *token == '\t') token++;
            char *end_token = token + strlen(token) - 1;
            while (end_token > token && (*end_token == ' ' || *end_token == '\t')) *end_token-- = '\0';

            if (*token == '\0') continue;

            errno = 0;

            if (col == 0 || col == 1) {
                temp_long = strtol(token, &endptr, 10);
                if (errno != 0 || endptr == token || *endptr != '\0' || temp_long > INT_MAX || temp_long < INT_MIN) {
                    parse_error = 1;
                } else {
                    if (col == 0) current_particle.type = (int)temp_long;
                    else current_particle.id = (int)temp_long;
                }
            } else {
                temp_double = strtod(token, &endptr);
                if (errno != 0 || endptr == token || *endptr != '\0') {
                    parse_error = 1;
                } else {
                    if (simple_output) {
                        // Simple output format: x,y,z,vx,vy,vz,mass
                        switch (col) {
                            case 2: current_particle.x = (real)temp_double; break;
                            case 3: current_particle.y = (real)temp_double; break;
                            case 4: current_particle.z = (real)temp_double; break;
                            case 5: current_particle.vx = (real)temp_double; break;
                            case 6: current_particle.vy = (real)temp_double; break;
                            case 7: current_particle.vz = (real)temp_double; break;
                            case 8: current_particle.mass = (real)temp_double; break;
                        }
                    } else {
                        // Non-simple output format: x,y,z,l,b,r,vx,vy,vz,mass,v_los,pmra,pmdec[,lambda,beta]
                        switch (col) {
                            case 2: current_particle.x = (real)temp_double; break;
                            case 3: current_particle.y = (real)temp_double; break;
                            case 4: current_particle.z = (real)temp_double; break;
                            case 5: current_particle.l = (real)temp_double; break;
                            case 6: current_particle.b = (real)temp_double; break;
                            case 7: current_particle.r = (real)temp_double; break;
                            case 8: current_particle.vx = (real)temp_double; break;
                            case 9: current_particle.vy = (real)temp_double; break;
                            case 10: current_particle.vz = (real)temp_double; break;
                            case 11: current_particle.mass = (real)temp_double; break;
                            case 12: current_particle.v_los = (real)temp_double; break;
                            case 13: current_particle.pm_ra = (real)temp_double; break;
                            case 14: current_particle.pm_dec = (real)temp_double; break;
                            case 15: 
                                if (total_cols >= 16) { // Only parse lambda if we have enough columns
                                    current_particle.lambda = (real)temp_double;
                                    has_lambda_beta = 1;
                                }
                                break;
                            case 16: 
                                if (has_lambda_beta) { // Only parse beta if we found lambda
                                    current_particle.beta = (real)temp_double;
                                }
                                break;
                        }
                    }
                }
            }

            if (parse_error) break;
            col++;
        }

        // Check if we have enough columns for the format
        int min_required_cols = simple_output ? 9 : 13; // 9 for simple (ignore,id,x,y,z,vx,vy,vz,mass), 13 for non-simple
        int max_expected_cols = simple_output ? 9 : (has_lambda_beta ? 16 : 14); // Account for optional lambda/beta
        
        if (!parse_error && col >= min_required_cols && col <= max_expected_cols) {
            if (!add_particle(collection, current_particle)) {
                fclose(file);
                return NULL;
            }
        } else if (!parse_error) {
            fprintf(stderr, "Warning: Line %ld has %d columns, expected %d-%d columns\n", 
                    line_number, col, min_required_cols, max_expected_cols);
        }
    } while (1);

    fclose(file);
    return collection;
}