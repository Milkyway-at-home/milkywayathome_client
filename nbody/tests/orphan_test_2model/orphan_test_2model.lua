
prng = DSFMT.create(3845024)

function makeHistogram()
   return HistogramParams.create()
end

function makePotential()
   local disk, halo, spherical
   disk = Disk.miyamotoNagai{
      mass = 4.45865888e5,
      scaleLength = 6.5,
      scaleHeight = 0.26
   }

   halo = Halo.logarithmic{
      vhalo = 73,
      scaleLength = 12.0,
      flattenZ = 1.0
   }

   spherical = Spherical.spherical{
      mass = 1.52954402e5,
      scale = 0.7
   }

   return Potential.create{ disk = disk,
                            halo = halo,
                            spherical = spherical
                         }
end

model1Bodies = 50000
model2Bodies = 10000
totalBodies = model1Bodies + model2Bodies

r0 = 0.2
dwarfMass = 12

-- This is also required
function makeContext()
   -- You can change the eps and timestep calculation however you want
   local eps = r0 / (10.0 * sqrt(totalBodies))
   local dt = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / dwarfMass)
   local ctx = NBodyCtx.create{ criterion = "NewCriterion",
                                useQuad = true,
                                theta = 1.0,
                                timeEvolve = 5.945,
                                eps2 = eps * eps,
                                timestep = dt
                              }
   return ctx
end

-- Also required
function makeBodies(ctx, potential)
   local firstModel, secondModel
   local orbitTimestep = ctx.timestep / 10.0
   local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
      velocity = Vector.create(-156, 79, 107),
      tstop = 6.0,
      dt = orbitTimestep
   }

   firstModel = predefinedModels.plummer(model1Bodies, {
                                            prng = prng,
                                            position = finalPosition,
                                            velocity = finalVelocity,
                                            mass = dwarfMass,
                                            ignore = false,
                                            scaleRadius = r0
                                         })

   secondModel = predefinedModels.plummer(model2Bodies, {
                                             prng = prng,
                                             position = finalPosition,
                                             velocity = finalVelocity,

                                             mass = 5000,
                                             scaleRadius = 0.9,
                                             ignore = true
                                          })

   return firstModel, secondModel
end

