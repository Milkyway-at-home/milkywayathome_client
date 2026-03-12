
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
      spherical = Spherical.hernquist{ mass = 67479.9, scale = 0.6 },
      disk      = Disk.doubleExponential{ mass = 224933, scaleLength = 6, scaleHeight = 0.3 },
      disk2     = Disk.none{ mass = 3.0e5 },
      halo      = Halo.nfw{ vhalo = 120, scaleLength = 22.25 }
   }
end

-- rdh edited
dec = 9.0   -- -- number of decimals to round to (default: 9.0)
function round(num, places)
    local mult = 10.0^(places)
    return floor(num * mult + 0.5) / mult
  end

function makeContext()
   return NBodyCtx.create{
      -- rdh edited
      dwarfn = 1,
      b           = {round( -44.328, dec)},
      r           = {round( 62.4,    dec)},
      vx          = {round( 21.99,   dec)},
      vy          = {round( -201.36, dec)},
      vz          = {round( 171.25,  dec)},
      sunGCDist   = SunGCDist,
      sunVelx     = SunVelx,
      sunVely     = SunVely,
      sunVelz     = SunVelz,
      timestep   = calculateTimestep(dwarfMass, dwarfRadius),
      timeEvolve = 3.945,
      eps2       = calculateEps2(nbody, dwarfRadius),
      criterion  = "sw93",
      useQuad    = true,
      theta      = 1.0,
      BestLikeStart = 0.95,
      BetaSigma     = 2.5,
      VelSigma      = 2.5,
      DistSigma     = 2.5,
      PMSigma       = 2.5,
      MomentumSigma = 2.5,
      BetaCorrect   = 1.111,
      VelCorrect    = 1.111,
      DistCorrect   = 1.111,
      PMCorrect     = 1.111,
      MomentumCorrect = 1.111,
      IterMax       = 6
   }
end

function makeBodies(ctx, potential)
   local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.5)),
      velocity  = Vector.create(-157, 78, 107),
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


