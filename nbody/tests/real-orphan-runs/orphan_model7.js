{  //test results are using double
"nbody-parameters-file": {
    "nbody-context": {
        "headline" : "model7-disk=miyamoto-nagai_halo=nfw_quad=true_criterion=sw93",
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
                "nfw" : {
                    "vhalo" : 109,
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
        "velocity" : [ -161, 85, 107 ],
        "position" : [ 218, 53.5, 28.7 ]
    },


}}

