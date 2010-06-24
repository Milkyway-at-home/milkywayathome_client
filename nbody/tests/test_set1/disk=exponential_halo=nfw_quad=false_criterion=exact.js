{  //test results are using double
"nbody-parameters-file": {
    "nbody-context": {
        "headline" : "disk=exponential_halo=nfw_quad=false_criterion=exact",
        "criterion" : "exact",
        "use-quadrupole-corrections" : false,
        "accuracy-parameter" : 1.0,
        "seed" : 0,

        "potential" : {
            "disk" : {
                "exponential" : {
                    "mass" : 4.45865888E5,
                    "scale-length" : 6.5,
                    "scale-height" : 0.26
                }
            },

            "spherical" : {
                "sphere" : {
                    "mass" : 1.52954402E5,
                    "r0-scale" : 0.7
                }
            },

            "halo" : {
                "nfw" : {
                    "vhalo" : 73,
                    "scale-length" : 12.0,
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
        "useGalC" : false,
        "angle-use-radians" : false,
        "velocity" : [ -156, 79, 107 ],
        "position" : [ 28.6, 218.0, 53.5 ]
    },


}}

