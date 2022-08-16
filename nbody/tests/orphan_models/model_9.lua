
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
      spherical = Spherical.plummer{ mass = 14852.322299938, scale = 0.25 },
      disk      = Disk.miyamotoNagai{ mass = 285750.313545277, scaleLength = 5.16, scaleHeight = 0.311 },
      disk2     = Disk.none{ mass = 3.0e5 },
      halo      = Halo.wilkinsonEvans{ mass = 2570079.320803279, scaleLength = 64.3 }
   }
end

function makeContext()
   return NBodyCtx.create{
      timestep   = calculateTimestep(dwarfMass, dwarfRadius),
      timeEvolve = 3.945,
      eps2       = calculateEps2(nbody, dwarfRadius,1),
      criterion  = "sw93",
      useQuad    = true,
      theta      = 1.0,
      BestLikeStart = 0.95,
      BetaSigma     = 2.5,
      VelSigma      = 2.5,
      DistSigma      = 2.5,
      BetaCorrect   = 1.111,
      VelCorrect    = 1.111,
      DistCorrect    = 1.111,
      IterMax       = 6
   }
end

function makeBodies(ctx, potential)
   local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.9)),
      velocity  = Vector.create(-179, 106, 109),
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
   return HistogramParams.create{
     --Orphan Stream coordinate transformation angles
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     
     -- ANGULAR RANGE AND NUMBER OF BINS
     lambdaStart = -50,
     lambdaEnd   = 50,
     lambdaBins  = 34,
     
     betaStart = -15,
     betaEnd   = 15,
     betaBins  = 1
}
end


