{  //test results are using double
"nbody-parameters-file": {
    "nbody-context": {
        "headline" : "model-triaxial-disk=miyamoto-nagai_halo=triaxial_quad=true_criterion=sw93",
        "criterion" : "sw93",
        "use-quadrupole-corrections" : true,
        "accuracy-parameter" : 1.0,
        "seed" : 0,

        "potential" : {
            "disk" : {
                "miyamoto-nagai" : {
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
                "triaxial" : {
                    "vhalo" : 116,
                    "scale-length" : 16.3,
                    "z-flattening" : 1.43,
                    "x-flattening" : 1.26,
                    "y-flattening" : 1.0,
                    "triaxial-angle" : 96,
                }
            }
        },

        "dwarf-model": {
            "plummer" : {
            "mass" : 16,
            "nbody" : 32768,
            "scale-radius" : 0.2,
            "time-orbit" : 4,
            "time-dwarf" : 3.945
            }
        }
    },

    "initial-conditions": {
        "useGalC" : false,
        "angle-use-radians" : false,
        "velocity" : [ -183, 101, 107 ],
        "position" : [ 218, 53.5, 29.5 ]
    },


}}

