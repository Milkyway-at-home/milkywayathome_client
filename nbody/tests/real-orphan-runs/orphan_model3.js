{  //test results are using double
"nbody-parameters-file": {
    "nbody-context": {
        "headline" : "model3-disk=exponential_halo=logarithmic_quad=true_criterion=sw93",
        "criterion" : "sw93",
        "use-quadrupole-corrections" : true,
        "accuracy-parameter" : 1.0,
        "seed" : 0,

        "potential" : {
            "disk" : {
                "exponential" : {
                    "mass" : 224933,
                    "scale-length" : 4
                }
            },

            "spherical" : {
                "sphere" : {
                    "mass" : 67479.9,
                    "r0-scale" : 0.6
                }
            },

            "halo" : {
                "logarithmic" : {
                    "vhalo" : 81,
                    "scale-length" : 12.0,
                    "z-flattening" : 1.0,
                }
            }
        },

        "dwarf-model": {
            "plummer" : {
            "mass" : 16,
            "nbody" : 1000,
            "scale-radius" : 0.2,
            "time-orbit" : 4,
            "time-dwarf" : 3.945
            }
        }
    },

    "initial-conditions": {
        "useGalC" : false,
        "angle-use-radians" : false,
        "velocity" : [ -152, 72, 106 ],
        "position" : [ 218, 53.5, 28.4 ]
    },


}}

