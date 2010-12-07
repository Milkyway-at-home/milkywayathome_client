{  //test results are using double
"nbody-parameters-file": {
    "nbody-context": {
        "headline" : "disk=miyamoto-nagai_halo=nfw_quad=false_criterion=bh86",
        "criterion" : "bh86",
        "use-quadrupole-corrections" : false,
        "accuracy-parameter" : 1.0,
        "seed" : 0,

        "potential" : {
            "disk" : {
                "miyamoto-nagai" : {
                    "mass" : 6.0E5,
                    "scale-length" : 4.0,
                    "scale-height" : 0.4
                }
            },

            "spherical" : {
                "sphere" : {
                    "mass" : 1.0E5,
                    "r0-scale" : 1.0
                }
            },

            "halo" : {
                "nfw" : {
                    "vhalo" : 155,
                    "scale-length" : 15.0,
                    "z-flattening" : 1.1,
                    "x-flattening" : 0.8,
                    "y-flattening" : 1.0,
                    "triaxial-angle" : 0.1
                }
            }
        },

        "dwarf-model": {
            "plummer" : {
                "mass" : 10,
                "nbody" : 100,
                "scale-radius" : 0.3,
                "time-orbit" : 4,
                "time-dwarf" : 3.945
            }
        }
    },

    "initial-conditions": {
        "use-galactic-coordinates" : false,
        "angle-use-radians" : false,
        "velocity" : [ -156, 79, 107 ],
        "position" : [ 28.6, 218.0, 53.5 ]
    },


}}

