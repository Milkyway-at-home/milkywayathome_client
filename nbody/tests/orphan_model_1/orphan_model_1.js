{
"nbody-parameters-file": {
    "nbody-context": {
        "headline" : "orphan model 1",
        "criterion" : "sw93",
        "use-quadrupole-corrections" : true,
        "accuracy-parameter" : 1.0,
        "seed" : 0,

        "time-orbit" : 4,
        "time-evolve" : 3.945,

        "potential" : {
            "disk" : {
                "type" : "exponential",
                "mass" : 224933,
                "scale-length" : 4
            },

            "spherical" : {
                "type" : "sphere",
                "mass" : 67479.9,
                "r0-scale" : 0.6
            },

            "halo" : {
                "type" : "nfw",
                "vhalo" : 155,
                "scale-length" : 22.25,
            }
        },

        "dwarf-model": [
            {
                "type" : "plummer",
                "mass" : 16,
                "nbody" : 1024,
                "scale-radius" : 0.2,
                "initial-conditions": {
                    "use-galactic-coordinates" : false,
                    "angle-use-radians" : false,
                    "velocity" : [ -170, 94, 108 ],
                    "position" : [ 218, 53.5, 28.8 ]
                }
            }
        ]
    },

    "histogram" : { }

}}

