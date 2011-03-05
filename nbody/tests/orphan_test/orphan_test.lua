

   -- Sample values
dwarfMass = 16
dwarfRadius = 0.2
reverseTime = 4.0
evolveTime = 3.945

-- Make sure to use the seed from the server
prng = DSFMT.create(argSeed)
nbody = 100


-- This is a required function. You get the defaults with no arguments
function makeHistogram()
   return HistogramParams.create()
end

-- This is a required function
function makePotential()
   local disk = Disk.miyamotoNagai{ mass = 4.45865888e5,
                                    scaleLength = 6.5,
                                    scaleHeight = 0.26
                                  }

   local halo = Halo.logarithmic{ vhalo = 73,
                                  scaleLength = 12.0,
                                  flattenZ = 1.0
                                }

   local spherical = Spherical.spherical{ mass = 1.52954402e5,
                                          scale = 0.7
                                        }

   local pot = Potential.create{ disk = disk,
                                 halo = halo,
                                 spherical = spherical
                               }
   return pot
end

-- This is also required
function makeContext()
   local r0 = dwarfRadius
   local eps = r0 / (10.0 * sqrt(nbody))
   local dt = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / dwarfMass)
   local ctx = NBodyCtx.create{ criterion = "sw93",
                                useQuad = true,
                                theta = 1.0,
                                timeEvolve = evolveTime,
                                eps2 = eps * eps,
                                timestep = dt
                              }
   return ctx
end

-- Also required
function makeBodies(ctx, potential)
   local iniVel = Vector.create(-156, 79, 107)
   local iniPos = Vector.create(218, 53.5, 28.6)

   local ic = InitialConditions.create{ context = ctx,
                                        velocity = iniVel,
                                        position = iniPos
                                      }

   local orbitTimestep = ctx.timestep / 10.0
   local finalPosition = reverseOrbit{ potential = potential,
                                       initialConditions = ic,
                                       tstop = reverseTime,
                                       dt = orbitTimestep
                                     }

   local plummerArgs = { prng = prng,
                         initialConditions = finalPosition,
                         mass = dwarfMass,
                         ignore = false,
                         scaleRadius =dwarfRadius
                       }
   local bodies = predefinedModels.plummer(nbody, plummerArgs)

   return bodies
end

