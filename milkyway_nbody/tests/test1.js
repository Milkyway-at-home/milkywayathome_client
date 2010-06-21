{  //test results are using double
"nbody-parameters-file": {
    "nbody-context": {
        "headline" : "Test 1, uses double precision",
        "criterion" : "new-criterion",
        "use-quadrupole-corrections" : false,
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
                "logarithmic" : {
                    "vhalo" : 73,
                    "scale-length" : 12.0,
                    "z-flattening" : 1.0
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







