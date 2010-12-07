{
"nbody-parameters-file": {
    "nbody-context": {
        "headline" : "orphan test tiny",
        "criterion" : "sw93",
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
                "z-flattening" : 1.0,
            }
        },

        "dwarf-model": [
            {
                "type" : "plummer",
                "mass" : 11.25,
                "nbody" : 256,
                "scale-radius" : 0.2,
                "initial-conditions": {
                    "use-galactic-coordinates" : false,
                    "angle-use-radians" : false,
                    "velocity" : [ -156, 79, 107 ],
                    "position" : [ 28.6, 218, 53.5 ]
                }
            }
        ]
    },

    "histogram" : { }

}}

