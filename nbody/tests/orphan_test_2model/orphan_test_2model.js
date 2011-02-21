{
"nbody-parameters-file": {
    "nbody-context": {
        "headline" : "orphan test 2 model",
        "criterion" : "new-criterion",
        "use-quadrupole-corrections" : true,
        "accuracy-parameter" : 1.0,
        "seed" : 0,

        "time-orbit" : 4,
        "time-evolve" : 3.945,

        "potential" : {
            "disk" : {
                "type" : "miyamoto-nagai",
                "mass" : 4.45865888E5,
                "scale-length" : 6.5,
                "scale-height" : 0.26
            },

            "spherical" : {
                "type" : "sphere",
                "mass" : 1.52954402E5,
                "r0-scale" : 0.7
            },

            "halo" : {
                "type" : "logarithmic",
                "vhalo" : 73,
                "scale-length" : 12.0,
                "z-flattening" : 1.0
            }
        },

        "dwarf-model": [
            {
                "type" : "plummer",
                "mass" : 12,
                "nbody" : 10000,
                "scale-radius" : 0.2,

                "reverse-orbit" : true,
                "ignore-final" : false,

                "initial-conditions": {
                    "use-galactic-coordinates" : false,
                    "angle-use-radians" : false,
                    "velocity" : [ -156, 79, 107 ],
                    "position" : [ 218, 53.5, 28.6 ]
                }
            },

            {
                "type" : "plummer",
                "mass" : 5000,
                "nbody" : 1024,
                "scale-radius" : 0.9,

                "reverse-orbit" : true,
                "ignore-final" : true,

                "initial-conditions": {
                    "use-galactic-coordinates" : false,
                    "angle-use-radians" : false,
                    "velocity" : [ -156, 79, 107 ],
                    "position" : [ 218, 53.5, 28.6 ]
                }
            }
        ]
    },

    "histogram" : {
        "phi"     : 128.79,
        "theta"   : 54.39,
        "psi"     : 90.70,
        "start"   : -50,
        "end"     : 50,
        "binsize" : 2.9411764705882355,
        "center"  : 0.0
    }

}}

