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

        // Skip empty lines and comments
        char *start = line_buffer;
        while (*start == ' ' || *start == '\t') start++;
        if (*start == '\0' || *start == '\n' || *start == '\r' || *start == '#') {
            continue;
        }

        ParticleData current_particle = {0}; // Initialize all fields to 0
        int col = 0;
        int parse_error = 0;
        long temp_long;
        double temp_double;
        int has_lambda_beta = 0;

        // Process each field in the CSV line
        char *line = line_buffer;
        char *field_start = line;
        char *field_end;
        
        // First pass: count actual columns
        int total_cols = 0;
        char *count_line = line_buffer;
        char *count_start = count_line;
        char *count_end;
        
        while (*count_line) {
            while (*count_start == ' ' || *count_start == '\t') count_start++;
            if (*count_start == '\0' || *count_start == '\n' || *count_start == '\r') break;
            
            count_end = count_start;
            while (*count_end && *count_end != ',' && *count_end != '\n' && *count_end != '\r') count_end++;
            
            char *trim_end = count_end - 1;
            while (trim_end > count_start && (*trim_end == ' ' || *trim_end == '\t')) trim_end--;
            
            if (count_start != count_end) total_cols++;
            
            count_line = count_end + 1;
            count_start = count_line;
        }

        // Determine format based on column count and simple_output flag
        int is_full_format = !simple_output && total_cols >= 17;
        int is_partial_format = !simple_output && total_cols >= 15 && total_cols < 17;
        int is_simple_format = simple_output && total_cols == 9;

        // Reset for actual parsing
        while (*line) {
            // Find the start of the field (skip leading whitespace)
            while (*field_start == ' ' || *field_start == '\t') field_start++;
            if (*field_start == '\0' || *field_start == '\n' || *field_start == '\r') break;
            
            // Find the end of the field (either comma or end of line)
            field_end = field_start;
            while (*field_end && *field_end != ',' && *field_end != '\n' && *field_end != '\r') field_end++;
            
            // Trim trailing whitespace
            char *trim_end = field_end - 1;
            while (trim_end > field_start && (*trim_end == ' ' || *trim_end == '\t')) trim_end--;
            *(trim_end + 1) = '\0';
            
            // Skip empty fields
            if (field_start == field_end) {
                line = field_end + 1;
                field_start = line;
                continue;
            }

            errno = 0;
            if (col == 0 || col == 1) {
                temp_long = strtol(field_start, &field_end, 10);
                if (errno != 0 || field_end != trim_end + 1 || temp_long > INT_MAX || temp_long < INT_MIN) {
                    parse_error = 1;
                    break;
                }
                if (col == 0) current_particle.type = (int)temp_long;
                else current_particle.id = (int)temp_long;
            } else {
                temp_double = strtod(field_start, &field_end);
                if (errno != 0 || field_end != trim_end + 1) {
                    parse_error = 1;
                    break;
                }
                if (is_simple_format) {
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
                } else if (is_partial_format) {
                    // Partial format: x,y,z,l,b,r,vx,vy,vz,mass,v_los,pmra,pmdec
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
                    }
                } else if (is_full_format) {
                    // Full format: x,y,z,l,b,r,vx,vy,vz,mass,v_los,pmra,pmdec,lambda,beta
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
                        case 15: current_particle.lambda = (real)temp_double; break;
                        case 16: current_particle.beta = (real)temp_double; break;
                    }
                }
            }

            col++;
            line = field_end + 1;
            field_start = line;
        }

        // Validate format and column count
        int min_required_cols, max_expected_cols;
        const char* format_name;
        
        if (is_simple_format) {
            min_required_cols = max_expected_cols = 9;
            format_name = "simple";
        } else if (is_partial_format) {
            min_required_cols = max_expected_cols = 15;
            format_name = "partial";
        } else if (is_full_format) {
            min_required_cols = max_expected_cols = 17;
            format_name = "full";
        } else {
            fprintf(stderr, "Warning: Line %ld has %d columns, which doesn't match any known format\n", 
                    line_number, total_cols);
            continue;
        }
        
        if (!parse_error && col >= min_required_cols && col <= max_expected_cols) {
            if (!add_particle(collection, current_particle)) {
                fclose(file);
                return NULL;
            }
        } else if (!parse_error) {
            fprintf(stderr, "Warning: Line %ld has %d columns, expected %d-%d columns for %s format\n", 
                    line_number, col, min_required_cols, max_expected_cols, format_name);
        }
    } while (1);

    fclose(file);
    return collection;
}