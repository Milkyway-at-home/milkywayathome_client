arg = { ... }
-- /* Copyright (c) 2016 Siddhartha Shelton */
assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")
argSeed = 34086709
prng = DSFMT.create(argSeed)

evolveTime       = tonumber(arg[1])
rev_ratio        = tonumber(arg[2])
rscale_l         = tonumber(arg[3])
light_r_ratio    = tonumber(arg[4])
mass_l           = tonumber(arg[5])
light_mass_ratio = tonumber(arg[6])

function round(num, places)
  local mult = 10.0^(places)
  return floor(num * mult + 0.5) / mult
end

dec = 9.0
evolveTime       = round( evolveTime,       dec )
rev_ratio        = round( rev_ratio,        dec )
rscale_l         = round( rscale_l,         dec )
light_r_ratio    = round( light_r_ratio,    dec )
mass_l           = round( mass_l,           dec )
light_mass_ratio = round( light_mass_ratio, dec )
model1Bodies = 2000
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.62"

revOrbTime = evolveTime / rev_ratio
dwarfMass = mass_l / light_mass_ratio
rscale_t  = rscale_l / light_r_ratio
rscale_d  = rscale_t *  (1.0 - light_r_ratio)
mass_d    = dwarfMass * (1.0 - light_mass_ratio)

--masses in order:
--spherical bulge: 3.4e10
--MN disk: 10e10
--log halo: 14e10

function makePotential()
   return  Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
   }
end


function get_timestep()
    --Mass of a single dark matter sphere enclosed within light rscale
    mass_enc_d = mass_d * (rscale_l)^3 * ( (rscale_l)^2 + (rscale_d)^2  )^(-3.0/2.0)

    --Mass of a single light matter sphere enclosed within dark rscale
    mass_enc_l = mass_l * (rscale_d)^3 * ( (rscale_l)^2 + (rscale_d)^2  )^(-3.0/2.0)

    s1 = (rscale_l)^3 / (mass_enc_d + mass_l)
    s2 = (rscale_d)^3 / (mass_enc_l + mass_d)
    
    --return the smaller time step
    if(s1 < s2) then
        s = s1
    else
        s = s2
    end
    
    -- I did it this way so there was only one place to change the time step. 
    t = (1 / 100.0) * ( pi_4_3 * s)^(1.0/2.0)
    return t
end


function makeContext()
   soften_length  = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
   return NBodyCtx.create{
      timeEvolve  = evolveTime,
      timestep    = get_timestep(),
      eps2        = calculateEps2(totalBodies, soften_length ),
      criterion   = "NewCriterion",
      useQuad     = true,
      theta       = 1.0
   }
end
--       position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
--       velocity  = Vector.create(-156, 79, 107),

-- Also required
-- for orphan: lbr in sun centered
-- position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6))
-- velocity  = Vector.create(-156, 79, 107),
      
-- dwarf starting positions
l  = 218
b  = 53.5
r  = 28.6
vx = -156 
vy = 79 
vz = 107
function makeBodies(ctx, potential)
  local firstModel
  local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(l, b, r)),
      velocity  = Vector.create(vx, vy, vz),
      tstop     = revOrbTime,
      dt        = ctx.timestep / 10.0
    }

  
-- print(lbrToCartesian(ctx, Vector.create(l, b, r)), Vector.create(vx, vy, vz))
print(finalPosition, finalVelocity)
  firstModel = predefinedModels.isotropic{
      nbody       = model1Bodies,
      prng        = prng,
      position    = finalPosition,
      velocity    = finalVelocity,
      mass1       = mass_l,
      mass2       = mass_d,
      scaleRadius1 = rscale_l,
      scaleRadius2 = rscale_d,
      ignore      = true
  }

--   print(finalPosition, finalVelocity)
  
  return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     --Orphan Stream coordinate transformation angles
     phi = 128.79,-- these are in degrees
     theta = 54.39,
     psi = 90.70,
     lambdaStart = -75,
     lambdaEnd = 50,
     lambdaBins = 50,
     betaStart = -40,
     betaEnd = 40,
     betaBins = 1
}
end


