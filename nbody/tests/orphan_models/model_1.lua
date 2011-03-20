
arg = {...}

seed = argSeed
nbody = arg[1]

assert(seed ~= nil, "Seed argument not set for test unit")
assert(nbody ~= nil, "Number of bodies not set for test unit")

prng = DSFMT.create(seed)

dwarfMass = 16
dwarfRadius = 0.2

function makePotential()
   return Potential.create{
      spherical = Spherical.spherical{ mass = 67479.9, scale = 0.6 },
      disk      = Disk.exponential{ mass = 224933, scaleLength = 4 },
      halo      = Halo.nfw{ vhalo = 155, scaleLength = 22.25 }
   }
end

function makeContext()
   return NBodyCtx.create{
      timestep   = calculateTimestep(dwarfMass, dwarfRadius),
      timeEvolve = 3.945,
      eps2       = calculateEps2(nbody, dwarfRadius),
      criterion  = "sw93",
      useQuad    = true,
      theta      = 1.0
   }
end

function makeBodies(ctx, potential)
   local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.8)),
      velocity  = Vector.create(-170, 94, 108),
      tstop     = 4.0,
      dt        = ctx.timestep / 10.0
   }

   return predefinedModels.plummer{
      nbody       = nbody,
      prng        = prng,
      position    = finalPosition,
      velocity    = finalVelocity,
      mass        = dwarfMass,
      scaleRadius = dwarfRadius,
      ignore      = false
   }
end

function makeHistogram()
   return HistogramParams.create()
end


