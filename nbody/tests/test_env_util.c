/* This file holds utility functions for the testing suite */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include "nbody.h"
#include "nbody_particle_data.h"
#include "nbody_dwarf_potential.h"
#include "nbody_types.h"
#include "nbody_lua.h"
#include "nbody_lua_types.h"
#include "nbody_lua_models.h"
#include "milkyway_alloc.h"
#include "milkyway_math.h"

/* Search paths configuration */
#define MAX_SEARCH_PATHS 10        /* Maximum number of paths to search */
#define MAX_PATH_LENGTH 1024       /* Maximum length of a path */

/**
 * @brief Find the test lua file using multiple search paths
 * 
 * @param filename The filename to search for
 * @return char* Path to the file if found, NULL otherwise (caller must free)
 */
char* find_lua_file(const char* filename) {
    if (filename == NULL || strlen(filename) == 0) {
        fprintf(stderr, "Error: Invalid filename provided to find_lua_file\n");
        return NULL;
    }
    
    char* result = mwCallocA(MAX_PATH_LENGTH, sizeof(char));
    if (result == NULL) {
        fprintf(stderr, "Error: Memory allocation failed in find_lua_file\n");
        return NULL;
    }
    
    char cwd[MAX_PATH_LENGTH];
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        fprintf(stderr, "Error: Failed to get current working directory: %s\n", strerror(errno));
        free(result);
        return NULL;
    }
    
    // Try common relative paths for different possible execution locations
    const char* relative_paths[] = {
        "../../../nbody/sample_workunits/test_env_lua",
        "/../nbody/sample_workunits/test_env_lua",
        "../nbody/sample_workunits/test_env_lua",
        "../sample_workunits/test_env_lua",
        "../test_env_lua"
    };
    
    size_t num_paths = sizeof(relative_paths) / sizeof(relative_paths[0]);
        
    // Try each relative path
    for (size_t i = 0; i < num_paths && i < MAX_SEARCH_PATHS; i++) {
        if (strlen(relative_paths[i]) == 0) {
            // Current directory
            snprintf(result, MAX_PATH_LENGTH, "%s/%s", cwd, filename);
        } else {
            snprintf(result, MAX_PATH_LENGTH, "%s%s/%s", 
                    (relative_paths[i][0] == '/') ? cwd : "", 
                    relative_paths[i], filename);
        }
        
        printf("Trying path: %s\n", result);
        fflush(stdout);
        
        if (access(result, F_OK) != -1) {
            printf("Found file at: %s\n", result);
            fflush(stdout);
            return result;
        }
    }
    
    // File not found
    fprintf(stderr, "Error: Could not find file '%s' in any search path\n", filename);
    free(result);
    return NULL;
}

/**
 * @brief Find the milkyway_nbody executable
 * 
 * @return char* Path to the executable if found, NULL otherwise (caller must free)
 */
char* find_milkyway_nbody() {
    char* result = mwCallocA(MAX_PATH_LENGTH, sizeof(char));
    if (result == NULL) {
        fprintf(stderr, "Error: Memory allocation failed in find_milkyway_nbody\n");
        return NULL;
    }
    
    char cwd[MAX_PATH_LENGTH];
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        fprintf(stderr, "Error: Failed to get current working directory: %s\n", strerror(errno));
        free(result);
        return NULL;
    }
    
    // Try common relative paths for the executable
    const char* relative_paths[] = {
        "../../bin",
        "../bin",
        "bin",
        "..",
        "../.."
    };
    
    size_t num_paths = sizeof(relative_paths) / sizeof(relative_paths[0]);
    
    // Try each relative path
    for (size_t i = 0; i < num_paths && i < MAX_SEARCH_PATHS; i++) {
        snprintf(result, MAX_PATH_LENGTH, "%s/%s/milkyway_nbody", cwd, relative_paths[i]);
        
        // For the current directory option, don't add an extra slash
        if (strcmp(relative_paths[i], ".") == 0) {
            snprintf(result, MAX_PATH_LENGTH, "%s/milkyway_nbody", cwd);
        }
        
        printf("Trying executable path: %s\n", result);
        fflush(stdout);
        
        if (access(result, X_OK) != -1) {
            printf("Found executable at: %s\n", result);
            fflush(stdout);
            return result;
        }
    }
    
    // If not found with relative paths, try to find it in PATH
    snprintf(result, MAX_PATH_LENGTH, "milkyway_nbody");
    if (system("which milkyway_nbody > /dev/null 2>&1") == 0) {
        printf("Found executable in PATH: milkyway_nbody\n");
        fflush(stdout);
        return result;
    }
    
    // Executable not found
    fprintf(stderr, "Error: Could not find milkyway_nbody executable\n");
    free(result);
    return NULL;
}

/**
 * @brief Run the nbody simulation with the given parameters
 * 
 * @param dwarf_params Array of parameters for the dwarf galaxy
 * @param lua_file Path to the Lua file
 * @return int 0 on success, non-zero on failure
 */
int run_nbody(const char** dwarf_params, const char* lua_file) {
    if (!dwarf_params) {
        fprintf(stderr, "Error: Null dwarf parameters array in run_nbody\n");
        return 1;
    }
    
    if (!lua_file) {
        fprintf(stderr, "Error: Null Lua file path in run_nbody\n");
        return 1;
    }
    
    // Validate the lua file exists
    if (access(lua_file, F_OK) == -1) {
        fprintf(stderr, "Error: Lua file '%s' does not exist\n", lua_file);
        return 1;
    }

    char command[MAX_PATH_LENGTH * 2];  // Larger buffer for the command
    char cwd[MAX_PATH_LENGTH];
    
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        fprintf(stderr, "Error: Failed to get current working directory: %s\n", strerror(errno));
        return 1;
    }
    
    // Find the milkyway_nbody executable
    char* bin_path = find_milkyway_nbody();
    if (!bin_path) {
        fprintf(stderr, "Error: Could not find milkyway_nbody executable\n");
        return 1;
    }
    
    // Validate all parameters before building command
    for (int i = 0; i < 6; i++) {
        if (!dwarf_params[i] || strlen(dwarf_params[i]) == 0) {
            fprintf(stderr, "Error: Invalid dwarf parameter at index %d\n", i);
            free(bin_path);
            return 1;
        }
    }
    
    // Build the command with proper escaping
    snprintf(command, sizeof(command), 
             "%s "
             "-f \"%s\" "
             "-o \"%s/output.out\" "
             "-z \"%s/output.hist\" "
             "-n 8 -b -w 1 -P -e 54231651 "
             "-i %s %s %s %s %s %s",
             bin_path, lua_file, cwd, cwd,
             dwarf_params[0], dwarf_params[1], dwarf_params[2], 
             dwarf_params[3], dwarf_params[4], dwarf_params[5]);
    
    printf("Running command: %s\n", command);
    fflush(stdout);
    
    int result = system(command);
    if (result != 0) {
        fprintf(stderr, "Error: Command execution failed with code %d\n", result);
    }
    
    // Free memory
    free(bin_path);
    
    return result;
}

/**
 * @brief Read parameters from a Lua file
 * 
 * @param input_lua_file Path to the Lua file
 * @param dwarf_params Array of parameters for the dwarf galaxy
 * @param nbody Output: total number of bodies
 * @param nbody_baryon Output: number of baryon bodies
 * @param comp1 Output: first dwarf component
 * @param comp2 Output: second dwarf component
 * @param timestep Output: simulation timestep
 * @return int 0 on success, non-zero on failure
 */
int read_lua_parameters(const char* input_lua_file, const char** dwarf_params, real* nbody, real* nbody_baryon, Dwarf** comp1, Dwarf** comp2, real* timestep) {
    if (!input_lua_file || !dwarf_params || !nbody || !nbody_baryon || !comp1 || !comp2 || !timestep) {
        fprintf(stderr, "Error: Null parameter provided to read_lua_parameters\n");
        return 1;
    }
    
    printf("Opening Lua state with file: %s\n", input_lua_file);
    fflush(stdout);
    
    // Validate the lua file exists
    if (access(input_lua_file, F_OK) == -1) {
        fprintf(stderr, "Error: Lua file '%s' does not exist\n", input_lua_file);
        return 1;
    }
    
    NBodyFlags nbf = {
        .inputFile = input_lua_file,
        .debugLuaLibs = TRUE,  
        .forwardedArgs = dwarf_params,
        .numForwardedArgs = 6
    };

    lua_State* L = nbOpenLuaStateWithScript(&nbf, NULL);
    if (!L) {
        fprintf(stderr, "Failed to open Lua state with script %s\n", input_lua_file);
        fflush(stdout);
        return 1;
    }

    // Register all necessary types and models
    registerNBodyTypes(L);
    registerPredefinedModelGenerators(L);

    // Get nbody values from global variables 
    lua_getglobal(L, "totalBodies");
    if (!lua_isnumber(L, -1)) {
        fprintf(stderr, "Error: totalBodies is not a number in Lua file\n");
        fflush(stdout);
        lua_close(L);
        return 1;
    }
    *nbody = lua_tonumber(L, -1);
    printf("totalBodies from Lua: %f\n", *nbody);
    fflush(stdout);
    lua_pop(L, 1);

    lua_getglobal(L, "totalLightBodies");
    if (!lua_isnumber(L, -1)) {
        fprintf(stderr, "Error: totalLightBodies is not a number in Lua file\n");
        fflush(stdout);
        lua_close(L);
        return 1;
    }
    *nbody_baryon = lua_tonumber(L, -1);
    printf("totalLightBodies from Lua: %f\n", *nbody_baryon);
    fflush(stdout);
    lua_pop(L, 1);

    // Validate nbody values before proceeding
    if (*nbody <= 0 || *nbody_baryon <= 0 || *nbody_baryon > *nbody) {
        fprintf(stderr, "Error: Invalid nbody values - nbody: %f, nbody_baryon: %f\n", *nbody, *nbody_baryon);
        fflush(stdout);
        lua_close(L);
        return 1;
    }

    // Get the makeBodies function
    lua_getglobal(L, "makeBodies");
    if (!lua_isfunction(L, -1)) {
        fprintf(stderr, "Error: makeBodies function not found in Lua script\n");
        fflush(stdout);
        lua_close(L);
        return 1;
    }

    // Get the makeContext function
    lua_getglobal(L, "makeContext");
    if (!lua_isfunction(L, -1)) {
        fprintf(stderr, "Error: makeContext function not found in Lua script\n");
        fflush(stdout);
        lua_close(L);
        return 1;
    }

    // Call makeContext
    printf("Calling makeContext...\n");
    fflush(stdout);
    if (lua_pcall(L, 0, 1, 0) != 0) {
        fprintf(stderr, "Error calling makeContext: %s\n", lua_tostring(L, -1));
        fflush(stdout);
        lua_close(L);
        return 1;
    }

    // Get the timestep from the context
    lua_getfield(L, -1, "timestep");
    if (!lua_isnumber(L, -1)) {
        fprintf(stderr, "Error: timestep is not a number in Lua context\n");
        fflush(stdout);
        lua_close(L);
        return 1;
    }
    *timestep = lua_tonumber(L, -1);
    lua_pop(L, 1);

    // Validate timestep
    if (*timestep <= 0.0 || !isfinite(*timestep)) {
        fprintf(stderr, "Error: Invalid timestep value: %f\n", *timestep);
        fflush(stdout);
        lua_close(L);
        return 1;
    }
    printf("Valid timestep value: %f\n", *timestep);
    fflush(stdout);

    // Push nil for potential
    lua_pushnil(L);

    // Call makeBodies
    printf("Calling makeBodies...\n");
    fflush(stdout);
    if (lua_pcall(L, 2, 1, 0) != 0) {
        fprintf(stderr, "Error calling makeBodies: %s\n", lua_tostring(L, -1));
        fflush(stdout);
        lua_close(L);
        return 1;
    }

    // Get the components from the model's table
    lua_getfield(L, -1, "components");
    if (!lua_istable(L, -1)) {
        fprintf(stderr, "Error: components table not found in model\n");
        fflush(stdout);
        lua_close(L);
        return 1;
    }

    // Get comp1
    lua_getfield(L, -1, "comp1");
    if (!lua_isuserdata(L, -1)) {
        fprintf(stderr, "Error: comp1 is not a userdata in Lua model\n");
        fflush(stdout);
        lua_close(L);
        return 1;
    }
    *comp1 = (Dwarf*)lua_touserdata(L, -1);
    lua_pop(L, 1);

    // Get comp2
    lua_getfield(L, -1, "comp2");
    if (!lua_isuserdata(L, -1)) {
        fprintf(stderr, "Error: comp2 is not a userdata in Lua model\n");
        fflush(stdout);
        lua_close(L);
        return 1;
    }
    *comp2 = (Dwarf*)lua_touserdata(L, -1);
    lua_pop(L, 1);

    // Validate the components
    if (!*comp1 || !*comp2) {
        fprintf(stderr, "Error: Invalid dwarf components\n");
        fflush(stdout);
        lua_close(L);
        return 1;
    }

    printf("Read parameters from Lua file:\n");
    printf("nbody: %f\n", *nbody);
    printf("nbody_baryon: %f\n", *nbody_baryon);
    printf("comp1 mass: %f\n", (*comp1)->mass);
    printf("comp1 scale length: %f\n", (*comp1)->scaleLength);
    printf("comp2 mass: %f\n", (*comp2)->mass);
    printf("comp2 scale length: %f\n", (*comp2)->scaleLength);
    printf("timestep: %f\n", *timestep);
    fflush(stdout);

    // Validate component parameters
    if ((*comp1)->mass <= 0 || (*comp1)->scaleLength <= 0 || 
        (*comp2)->mass <= 0 || (*comp2)->scaleLength <= 0) {
        fprintf(stderr, "Error: Invalid component parameters\n");
        fflush(stdout);
        lua_close(L);
        return 1;
    }

    lua_close(L);
    return 0;
}