
arg = {...}

seed = argSeed
nbody = arg[1]

assert(seed ~= nil, "Seed argument not set for test unit")
assert(nbody ~= nil, "Number of bodies not set for test unit")

prng = DSFMT.create(seed)

dwarfMass = 16
dwarfRadius = 0.2

dec = 9.0   -- -- number of decimals to round to (default: 9.0)
function round(num, places)
    local mult = 10.0^(places)
    return floor(num * mult + 0.5) / mult
  end
rscale_l            = {round( 2.9,     dec)}  -- Baryonic Radius (kpc)
light_r_ratio       = {round( 0.2,     dec)}  -- Baryonic Radius / (Baryonic Radius + Dark Matter Radius)
mass_l              = {round( 2429.198,dec)}  -- Baryonic Mass (Structure Mass Units)
light_mass_ratio    = {round( 0.0830,  dec)}  -- Baryonic Mass / (Baryonic Mass + Dark Matter Mass)
orbit_parameter_l   = {round( 302.801, dec)}  -- Galactocentric l
orbit_parameter_b   = {round( -44.328, dec)}  -- Galactocentric b
orbit_parameter_r   = {round( 62.4,    dec)}  -- Galactocentric r
orbit_parameter_vx  = {round( 21.99,   dec)}  -- Galactocentric vx
orbit_parameter_vy  = {round( -201.36, dec)}  -- Galactocentric vy
orbit_parameter_vz  = {round( 171.25,  dec)}  -- Galactocentric vz

dwarf = Dwarf.plummer{mass = dwarfMass, scaleLength = dwarfRadius}

function makePotential()
   return Potential.create{
      spherical = Spherical.hernquist{ mass = 67479.9, scale = 0.6 },
      disk      = Disk.freeman{ mass = 224933, scaleLength = 4 },
      disk2     = Disk.none{ mass = 3.0e5 },
      halo      = Halo.nfw{ vhalo = 155, scaleLength = 22.25 }
   }
end

function makeContext()
   return NBodyCtx.create{
      dwarfn = 1,
      b           = orbit_parameter_b,
      r           = orbit_parameter_r,
      vx          = orbit_parameter_vx,
      vy          = orbit_parameter_vy,
      vz          = orbit_parameter_vz,
      sunGCDist   = SunGCDist,
      sunVelx     = SunVelx,
      sunVely     = SunVely,
      sunVelz     = SunVelz,
      timestep   = calculateTimestep(dwarfMass, dwarfRadius),
      timeEvolve = 3.945,
      eps2       = calculateEps2Dwarf(dwarf, nbody),
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


