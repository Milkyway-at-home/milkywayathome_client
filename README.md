# Milkyway@Home Client

[![Linux Build](https://travis-ci.org/Milkyway-at-home/milkywayathome_client.svg?branch=master)](https://travis-ci.org/Milkyway-at-home/milkywayathome_client)
[![MinGW Build](https://travis-ci.org/Milkyway-at-home/milkywayathome_client.svg?branch=travis-xcompile)](https://travis-ci.org/Milkyway-at-home/milkywayathome_client)

> **Note:** CMake version 4.0 and later are **not currently supported**.

---

## Table of Contents

- [Nbody](#nbody)
- [Compiling Nbody](#instructions-for-compiling-nbody)
- [Running Nbody](#running-nbody-options)
- [Input Lua File Dwarf Model Options](#input-lua-file-dwarf-model-options)
- [N-Body CMAKE Flags](#n-body-cmake-flags)
- [Tests](#tests)
- [Separation](#separation)
- [TAO](#tao)
- [Random Notes](#random-notes)

Nbody
---
- Simulations are described with Lua input files which can be used
  to produce an arbitrary initial configuration of particles. 

- Number of particles can be indicated in the Lua input file as 
  a total number of bodies where half will be baryons and half 
  will be dark matter particles or as the total number of bodies
  with the number of baryons as an extra parameter  

- Various options are available for applying external potentials
  to a system.

- Graphics can be run separately and attach to existing simulations,
  or can be launched at the same time with the --visualizer argument
  to the main process.

- N-body videos can be produced by using a separate program to
  record OpenGL. A wrapper script that uses this can be used as
  the --visualizer-bin argument to record a video of the
  visualization. An example script is at tools/RecordNBodyVideo.sh

- Consistent N-body results between different systems require crlibm
  and SSE2 (at least on x86, not sure about other architectures)

- Returning nil from makePotential() for N-body will run the
  simulation without an external potential

- Device information is exposed to the workunit through the
  deviceInfo table if it is used.

- **Bar code currently does not pass all tests**


Instructions for Compiling Nbody
---
Step 0.  Ensure proper packages are installed

    (For Ubuntu) sudo apt-get install mingw-w64 cmake
    (OpenGL)     sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
    (NCurses)    sudo apt-get install libncurses5-dev libncursesw5-dev
    (OpenSSL)    sudo apt-get install libssl-dev

Step 1.  Download all necessary files (Only need to git submodule if cross compiling with BOINC)
```
git clone https://github.com/Milkyway-at-home/milkywayathome_client.git
cd milkywayathome_client
git submodule update --init --recursive
```
NOTE: If you are running on WSL (Windows Subsystem Linux), you may need to run the following commands
```
git submodule sync
git submodule init
git submodule update
```
Step 2.  Compile Nbody
```
./build_client
```
Step 3.  Run a Nbody Simulation
```
./run_nbody
```

Running N-Body Options
---
The type of run is set by setting one of the following flags to `true`:  
`run`, `run_compare`, `compare_only`, or `get_flag_list`.

### Command-Line Options

| Option | Description |
|--------|-------------|
| `-f`   | Path to input LUA file |
| `-o`   | Path to bodies output file |
| `-z`   | Path to histogram output file |
| `-h`   | Path to histogram input file (used with `run_compare` only) |
| `-e`   | Seed |
| `-n`   | Number of threads to use for simulation |
| `-P`   | Print the percentage of progress of the simulation to standard output |
| `-u`   | Runs the visualizer (may require additional packages and compilation with OpenGL) |
| `-p`   | Simulation paramters list (6, 7, 8, 12, 13 or 14 arguments)

#### `-p` Options

- **Required 6 arguments:**  
  `[1] Forward Time, [2] Time Ratio, [3] Baryon Scale Radius, [4] Radius Ratio, [5] Baryon Mass, [6] Mass Ratio`
- **If 7 arguments:**  
  - If `manual_bodies = true`: `[7] Manual Bodies Input File`  
  - Else: `[7] LMC_mass`
- **If 8 arguments:**  
  `[7] LMC Mass, [8] Manual Bodies Input File`
- **If 12 arguments:**  
  `[7] l, [8] b, [9] r, [10] vx, [11] vy, [12] vz`
- **If 13 arguments:**  
  - If `manual_bodies = true`: `[13] Manual Bodies Input File`  
  - Else: `[13] LMC_mass`
- **If 14 arguments:**  
  `[13] LMC Mass, [14] Manual Bodies Input File`

#### Likelihood Comparison Flags

| Flag | Description |
|------|-------------|
| `-s` | Compare using only EMD and cost component |
| `-S` | Use EMD, cost, beta dispersion |
| `-V` | Use EMD, cost, velocity dispersion |
| `-D` | Use EMD, cost, beta dispersion and velocity dispersion |
| `-A` | Compare all components of the likelihood |

---

## Input Lua File Dwarf Model Options

### Double Component Model

- **Plummer:** `{mass, scaleLength}`
- **NFW:** `{mass, scaleLength}`
- **General Hernquist:** `{mass, scaleLength}`
- **Cored:** `{mass, scaleLength, r1, rc}`  _***Not fully tested***_

### Single Component Model

- **Plummer:** `{nbody, mass, scaleRadius, position, velocity, ignore, prng}`
- **NFW:** `{nbody, mass, rho_0, scaleRadius, position, velocity, ignore, prng}`
- **Hernquist:** `{nbody, mass, radius, a, position, velocity, ignore, prng}`

---

## N-Body CMAKE Flags

| Flag | Values | Description |
|------|--------|-------------|
| `DCMAKE_BUILD_TYPE`      | Debug, Release, RelWithDebInfo, MinSizeRel | Set to `Release` for a normal build. Other options include debugging information. |
| `DNBODY_DEV_OPTIONS`     | ON, OFF | Set to `ON` for developer options. `OFF` to use client-side parameter files. |
| `DNBODY_GL`              | ON, OFF | Builds the visualizer. Requires additional OpenGL packages. |
| `DBOINC_APPLICATION`     | ON, OFF | Cross-compile with BOINC. |
| `DSEPARATION`            | ON, OFF | Option for building the Separation code. Defaults to `OFF`. |
| `DDOUBLEPREC`            | ON, OFF | Enable double-precision floating point calculation. |
| `DNBODY_OPENMP`          | ON, OFF | Build the algorithm single-threaded (`OFF`) or multithreaded (`ON`). |
| `DNBODY_OPENCL`          | ON, OFF | Build with OpenCL libraries to support running N-Body on GPUs. |

Tests
---
  Tests can be run by running:
  ```
  $ make test
  ```
  However this runs all of the tests, which takes forever. You can run
  (from the tests directory) some core functionality tests with:
  ```
  $ make check
  ```
  Other tests  can be run with a certain number of bodies depending on
  how long you want to wait with:
  ```
  $ make test_${n}
  ```
  Currently n = 100, 1024, 10000 are available.

  Single tests can be run with:
  ```
  $ ctest -R <Test_Name> 
  ```
  Get a more versbose output with:
  ```
  $ ctest -R <Test_Name> -VV
  ```
  If only 25 tests are running instead of 88 tests, you are missing libraries (check Step 0 for compiling N-body)

Separation
---
- separation will do a separation after the integration if given an
  output file. There is also an argument to set the random number seed.

TAO
---
TODO: update for latest tao version

- Maximum Likelihood Evaluation Code for running milkyway separation program

- Note: Lua files for TAO searches are different from those used by the separation code.

- The terminal output from this program appears confusing since it mixes the output of each separation run with that of TAO.  Using the linux ">" operator to port output to a file only takes the TAO output, making it much clearer.

To run call:
```
  $ ./TAO <options>

  General required options:
    --separation "<path/to/separation_binary>"
    --stars "<path/to/stars_file>"
    --params "<path/to/search_paramaters_file>"
    --search_type <options>
      search type options:
        de    - differential evolution
        ps    - particle swarm
        snm   - synchronous newton method
        gd    - gradient descent
        cgd   - conjugate gradient descent
        sweep - paramater sweep

  Search Specific Options:
    de:
      optional:
        --population_size <int>         (default:200)
        --maximum_iterations <int>        (default:will run forever - Ctrl-C to kill)
        --maximum_created <int>         (default:will run forever - Ctrl-C to kill)
        --maximum_reported <int>        (default:will run forever - Ctrl-C to kill)
        --parent_scaling_factor <float>     (default:1.0)
        --differential_scaling_factor <float> (default:1.0)
        --crossover_rate <float>        (default:0.5)
        --int_pairs <int>           (default:1)
        --parent_selection <option>       (defualt:best)
          options:
            best
            random
            current-to-best
            current-to-random
        --recombination_selection <option>    (default:binary)
          options:
            binary
            exponential
            sum
            none
    ps:
      optional:
        --population_size <int>       (default:200)
        --maximum_iterations <int>      (default:will run forever - Ctrl-C to kill)
        --maximum_created <int>       (default:will run forever - Ctrl-C to kill)
        --maximum_reported <int>      (default:will run forever - Ctrl-C to kill)
        --inertia <float>         (default:0.75)
        --global_best_weight <float>    (default:1.5)
        --local_best_weight <float>     (default:1.5)
        --initial_velocity_scale <float>  (default:0.25)
    snm:
      required:
        --iterations <int>
      optional:
        --rand <double>     (randomizes the search parameters by +- the given percent)
    gd:
      required:
        --iterations <int>
      optional:
        --loop1_max <int>   (default:300 iterations)
        --loop2_max <int>   (default:300 iterations)
        --nquad <int>     (default:4 iterations for loop 3)
        --tol <double>      (default:1e-6 for tolerance of dstar in loop 3)
        --min_threshold <double_1, double_2, ... , double_n>
                    (default:line search will not quit if the input direction is very small)
        --rand <double>     (randomizes the search parameters by +- the given percent)

    gd:
      required:
        --iterations <int>
        --cgd_reset <int> **roughly speaking this should be the number of paramaters...
      optional:
        --loop1_max <int>   (default:300 iterations)
        --loop2_max <int>   (default:300 iterations)
        --nquad <int>     (default:4 iterations for loop 3)
        --tol <double>      (default:1e-6 for tolerance of dstar in loop 3)
        --min_threshold <double_1, double_2, ... , double_n>
                    (default:line search will not quit if the input direction is very small)
        --rand <double>     (randomizes the search parameters by +- the given percent)
```

Random notes:
---

 - All give usage with --help/-? arguments

make nbody_release and make separation_release will produce release
tarballs if git and xz are installed and found.

- Make sure when building with MSVC to set built to use Multithreaded
  (/MT) for the builds of the various libraries
