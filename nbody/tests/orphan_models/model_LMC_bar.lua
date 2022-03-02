
arg = {...}

seed = argSeed
nbody = arg[1]

--assert(seed ~= nil, "Seed argument not set for test unit")
--assert(nbody ~= nil, "Number of bodies not set for test unit")

prng = DSFMT.create(seed)

dwarfMass = 16
dwarfRadius = 0.2
LMCMASS = 449865.888
LMCSCALE = 15.0

function makePotential()
   return Potential.create{
      spherical = Spherical.hernquist{ mass = 67479.9, scale = 0.6 },
      disk      = Disk.miyamotoNagai{ mass = 198040, scaleLength = 4, scaleHeight = 0.26 },
      --disk      = Disk.miyamotoNagai{ mass = 224933, scaleLength = 4, scaleHeight = 0.26 },
      disk2     = Disk.orbitingBar{ mass = 2.689340798e4, scaleLength = 5.4, patternSpeed = 250.61, startAngle = 0.488692},
      --disk2     = Disk.none{mass = 2},
      halo      = Halo.logarithmic{ vhalo = 155, scaleLength = 22.25, flattenZ = 1.1 }
   }
end

function makeContext()
   return NBodyCtx.create{
      timestep   = calculateTimestep(dwarfMass, dwarfRadius),
      timeEvolve = 3.945,
      timeBack   = 3.945,
      eps2       = calculateEps2(nbody, dwarfRadius),
      criterion  = "sw93",
      useQuad    = true,
      theta      = 1.0,
      BetaSigma     = 2.5,
      VelSigma      = 2.5,
      DistSigma     = 2.5,
      BetaCorrect   = 1.111,
      VelCorrect    = 1.111,
      DistCorrect   = 1.111,
      IterMax       = 6,
      LMC           = true,
      LMCmass       = LMCMASS,
      LMCscale      = LMCSCALE,
      LMCDynaFric   = true
   }
end

function makeBodies(ctx, potential)
   local finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity = reverseOrbit_LMC{
      potential   = potential,
      position    = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.8)),
      velocity    = Vector.create(-170, 94, 108),
      LMCposition = Vector.create(-1.1, -41.1, -27.9),
      LMCvelocity = Vector.create(-57, -226, 221), 
      LMCmass     = LMCMASS,
      LMCscale    = LMCSCALE,
      LMCDynaFric = true,
      ftime       = 3.945,
      tstop       = 3.945,
      dt          = ctx.timestep / 10.0,
      sunGCDist = 8.0
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


