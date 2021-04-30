Linux: [![Linux Build](https://travis-ci.org/Milkyway-at-home/milkywayathome_client.svg?branch=master)](https://travis-ci.org/Milkyway-at-home/milkywayathome_client)
MinGW: [![MinGW Build](https://travis-ci.org/Milkyway-at-home/milkywayathome_client.svg?branch=travis-xcompile)](https://travis-ci.org/Milkyway-at-home/milkywayathome_client)

Separation
--------------------------------------------------------------------------------
- separation will do a separation after the integration if given an
  output file. There is also an argument to set the random number seed.

Nbody
--------------------------------------------------------------------------------
- Simulations are described with Lua input files which can be used
  to produce an arbitrary initial configuration of particles.

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

Tests
-----------------
  Tests can be run by running:
  $ make test

  However this runs all of the tests, which takes forever. You can run
  (from the tests directory) some core functionality tests with:
  $ make check

  Other tests  can be run with a certain number of bodies depending on
  how long you want to wait with:

  $ make test_${n}

  Currently n = 100, 1024, 10000 are available.

TAO
--------------------------------------------------------------------------------
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

---------------------------------------------------------------------------------------------------

Random notes:

 - All give usage with --help/-? arguments

make nbody_release and make separation_release will produce release
tarballs if git and xz are installed and found.

- Make sure when building with MSVC to set built to use Multithreaded
  (/MT) for the builds of the various libraries

Instructions for Cross Compiling (With BOINC on)
---------------------------------------------------------------------------------------------------
Step 0.  Ensure proper packages are installed
    (For Ubuntu) `sudo apt-get install mingw-w64 cmake`

Step 1.  Download all necessary files.
```
git clone https://github.com/Milkyway-at-home/milkywayathome_client.git
cd milkywayathome_client
git submodule init  
git submodule update --recursive
```
Step 2.  Set up build directory
```
cd ../
mkdir build
cd build
```
Step 3.  Set up build files with CMake and Compile
(Only Confirmed to work with Separation)
```
cmake -DCMAKE_TOOLCHAIN_FILE="../milkywayathome_client/cmake_modules/MinGW32-Cross-Toolchain.cmake" -DBUILD_32=<ON/OFF> -       DSEPARATION_OPENCL=<ON/OFF> -DSEPARATION_STATIC=ON -DOPENCL_LIBRARIES=<Path to OpenCL.lib for 32 or 64 bit depending on build> -DOPENCL_INCLUDE_DIRS=<Path to OpenCL "include" files> ../milkywayathome_client/
make
```
* These instructions may need to be updated as BOINC updates their software.
