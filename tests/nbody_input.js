{  // I'm a C++ style comment
"nbody-parameters-file": {        /* I'm a C style comment */
    "nbody-context": {
        "outfile": "out",
        "headline" : "Captain Picard",
        "cheese count" : 246,
        "I'm another parameter that doesn't exist" : null,
        "criterion" : "sw93",        // I'm an enum!
        "accuracy-parameter" : 12.0,
        "potential-softening" : 0.0,
        "use-quadrupole-corrections" : false,    /* If true, use quad moments */
        "allow-incest" : false,
        "freq" : 1.0,
        "freqout" : 1.0,
        "nbody" : 300,       // I'm an int
        "seed" : 123,       /* Random number seed for test run */
        "tstop" : 10,       //evolution time
        "nbody" : 100,      /* Number of particles for test run */

        "tstop" : 2.0,      /* Time to stop integration */
        "dtout" : 0.25,     /* Data output interval */

        "potential" : {
            "disk" : {
                "exponential" : {   // I say what model I am
                    //We're parameters used by this model
                    "mass" : 3.4,        // I'm a double
                    "scale-length" : 5,  //I'm a double without a decimal
                    "scale-height" : 9.7
                }
            },

            "spherical" : {
                "sphere" : {
                    "mass" : 1.0,
                    "r0-scale" : 1.0,
                    "ITS A FAKE" : true  // I produce a warning
                }
            },

            "halo" : {
                "nfw" : {
                    "vhalo" : 1.0,
                    "scale-length" : 9000.0
                }
            }
        },

        "dwarf-model": {
            "plummer" : {
                "mass" : 42,
                "scale-radius" : 9,
                "time-orbit" : null, // I'm optional. I could also be omitted.
                "time-dwarf" : 1.2e9,
                "IT ISN'T REAL" : true
            }
        }
    },

    "tree" : {
        "rsize" : 4       // Size of initial t.root cell
    },

    "rings" : {
        "a" : 1.0,
        "M" : 1.0,
        "XC" : 1.0,
        "YC" : 1.0,
        "ZC" : 1.0
    }


    /*
      "initial-coordinates": {
      // From the vhalo = 73 model result from Newberg et al 2009
      // The l,b,r and vx, vy, vz are hard coded for the Orphan project
      "angle-use-radians" : false,
      "lstart" : 218.0,
      "bstart" : 53.5,
      "Rstart" : 28.6,

      "Xinit" : -156,
      "Yinit" : 79,
      "Zinit" : 107,

      "sunGCDist" : 8.0
      }
    */

}}







