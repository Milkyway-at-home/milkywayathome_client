{"nbody-parameters-file": {
    "nbody-parameters": {
        "PluMass": 10,      // mass
        "r0" : 10,          // radius
        "orbittstop" : 3,   // how far back the orbit goes
        "theta" : 1.0,      /* Cell subdivision tolerence */
        "eps" : 0.025,      /* Potential softening parameter */
        "dtime" : 0.03125,  /* Integration time-step */
        "tstop" : 2.0,      /* Time to stop integration */
        "dtout" : 0.25      /* Data output interval */
    },

    // Calculate starting galactic coordinates
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
    },

    "nbody-context": {
        "outfile": "out",
        "headline" : "headline",

        "model" : "bh86",

        "usequad" : false,    /* If true, use quad moments */
        "allow-incest" : false,

        "freq" : 1.0,
        "freqout" : 1.0,

        "seed" : 123,         /* Random number seed for test run */
        "tstop" : 10,      //evolution time
        "nbody" : 100      //number of bodies
    },

    "tree" :{
        "rsize" : 4       /* Size of initial t.root cell */
    }

}}

