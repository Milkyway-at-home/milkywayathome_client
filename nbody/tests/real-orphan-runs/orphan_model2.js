{  //test results are using double
"nbody-parameters-file": {
    "nbody-context": {
        "headline" : "model2-disk=exponential_halo=NFW_quad=true_criterion=sw93",
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
                "nfw" : {
                    "vhalo" : 120,
                    "scale-length" : 22.25
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
        "velocity" : [ -157, 78, 107 ],
        "position" : [ 218, 53.5, 28.5 ]
    },


}}

