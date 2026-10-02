
arg = { ... }

--prng = DSFMT.create(3845024)
prng = DSFMT.create(argSeed)
nbody = arg[1]

assert(nbody ~= nil, "Number of bodies not set for test unit")


function makePotential()
   return Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
   }
end

model1Bodies = nbody / 2
model2Bodies = nbody / 2

totalBodies = model1Bodies + model2Bodies

r0, r02 = 0.2, 0.5
dwarfMass, dwarfMass2 = 12, 190

encMass = plummerTimestepIntegral(r0, r02, dwarfMass2, 1e-7)

-- This is also required
function makeContext()
   -- secondModel below is generated with ignore=true (particle type -1),
   -- so eps2_index must cover both types, not just 1 -- reusing the same
   -- single computed value for all four matrix entries, same as this
   -- project's gpu_test_checkpoint.lua fix.
   local eps2_val = calculateEps2(totalBodies, r0)
   return NBodyCtx.create{
      timeEvolve = 5.945,
      timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (encMass + dwarfMass)),
      eps2       = {eps2_val, eps2_val, eps2_val, eps2_val},
      eps2_index = {1, -1},
      eps2_size  = 2,
      criterion  = "TreeCode",
      useQuad    = true,
      theta      = 1.0
   }
end

-- Also required
function makeBodies(ctx, potential)
   local firstModel, secondModel
   local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
      velocity  = Vector.create(-156, 79, 107),
      tstop     = 6.0,
      dt        = ctx.timestep / 10.0
   }

   firstModel = predefinedModels.plummer{
      nbody       = model1Bodies,
      prng        = prng,
      position    = finalPosition,
      velocity    = finalVelocity,
      mass        = dwarfMass,
      scaleRadius = r0,
      ignore      = false
   }

   secondModel = predefinedModels.plummer{
      nbody       = model2Bodies,
      prng        = prng,
      position    = finalPosition,
      velocity    = finalVelocity,
      mass        = dwarfMass2,
      scaleRadius = r02,
      ignore      = true
   }
   return firstModel, secondModel
end

function makeHistogram()
   return HistogramParams.create()
end

