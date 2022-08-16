
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
      spherical = Spherical.hernquist{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.18972480e5, scaleLength = 6.5, scaleHeight = 0.26 },
      disk2     = Disk.orbitingBar{ mass = 2.689340798e4, scaleLength = 5.4, patternSpeed = 250.61, startAngle = 0.488692},
      halo      = Halo.logarithmic{ vhalo = 74.61, scaleLength = 12.0, flattenZ = 1.0 }
   }
end

evolveTime = 3.945
best_like_start = 0.95

evolveTime = (2.0 - best_like_start) * evolveTime --making it evolve to end of best-likelihood window
eff_best_like_start = best_like_start / (2.0 - best_like_start) --correct for changed evolve time

function makeContext()
   return NBodyCtx.create{
      timestep   = calculateTimestep(dwarfMass, dwarfRadius),
      timeEvolve = evolveTime,
      timeBack    = 3.945,
      eps2       = calculateEps2(nbody, dwarfRadius,0),
      criterion  = "sw93",
      useQuad    = true,
      theta      = 1.0,
      useBestLike = true,
      BestLikeStart = eff_best_like_start,
      BetaSigma     = 2.5,
      VelSigma      = 2.5,
      DistSigma     = 2.5,
      BetaCorrect   = 1.111,
      VelCorrect    = 1.111,
      DistCorrect   = 1.111,
      IterMax       = 6,
      calibrationRuns = 2
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


