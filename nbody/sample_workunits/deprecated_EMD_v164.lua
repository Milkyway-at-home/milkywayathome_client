-- /* Copyright (c) 2016 Siddhartha Shelton */

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- Client side Lua file
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
                
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
totalBodies           = 20000   -- -- NUMBER OF BODIES           -- --
nbodyLikelihoodMethod = "EMD"   -- -- HIST COMPARE METHOD        -- --
nbodyMinVersion       = "1.64"  -- -- MINIMUM APP VERSION        -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 


-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- PARAMETER SETTINGS   -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- -- -- -- -- -- -- -- -- HISTOGRAM   -- -- -- -- -- -- -- -- -- -- -- -- --
lda_bins        = 50      -- number of bins in lamdba direction
lda_lower_range = -150    -- lower range for lambda
lda_upper_range = 150     -- upepr range for lamdba

bta_bins        = 1       -- number of beta bins. normally use 1 for 1D hist
bta_lower_range = -15     -- lower range for beta
bta_upper_range = 15      -- upper range for beta
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

-- -- -- -- -- -- -- -- -- AlGORITHM OPTIONS -- -- -- -- -- -- -- --
use_best_likelihood  = true    -- use the best likelihood return code
best_like_start      = 0.98    -- what percent of sim to start
use_vel_disps        = true    -- use velocity dispersions in likelihood
        

-- -- -- -- -- -- -- -- -- DWARF STARTING LOCATION   -- -- -- -- -- -- -- --
l  = 218
b  = 53.5
r  = 28.6
vx = -156 
vy = 79 
vz = 107
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

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
    
    tmp = sqr(1/10.0) * sqrt((pi_4_3 * cube(rscale_d)) / (mass_l + mass_d))
--     print('timestep ', t, tmp)
    
    return t
end


function makeContext()
   soften_length  = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
   return NBodyCtx.create{
      timeEvolve  = evolveTime,
      timestep    = get_timestep(),
      eps2        = calculateEps2(totalBodies, soften_length ),
      criterion   = "TreeCode",
      useQuad     = true,
      useBestLike = use_best_likelihood,
      useVelDisp  = use_vel_disps,
      BestLikeStart = best_like_start,
      theta       = 1.0
   }
end


function makeBodies(ctx, potential)
    local firstModel
    local finalPosition, finalVelocity
    
    finalPosition, finalVelocity = reverseOrbit{
        potential = potential,
        position  = lbrToCartesian(ctx, Vector.create(l, b, r)),
        velocity  = Vector.create(vx, vy, vz),
        tstop     = revOrbTime,
        dt        = ctx.timestep / 10.0
    }
    
    firstModel = predefinedModels.mixeddwarf{
        nbody       = totalBodies,
        prng        = prng,
        position    = finalPosition,
        velocity    = finalVelocity,
        comp1       = Dwarf.plummer{mass = mass_l, scaleLength = rscale_l}, -- Dwarf Options: plummer, nfw, general_hernquist
        comp2       = Dwarf.plummer{mass = mass_d, scaleLength = rscale_d}, -- Dwarf Options: plummer, nfw, general_hernquist
        ignore      = true
    }

  return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     --Orphan Stream coordinate transformation angles
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     
     -- ANGULAR RANGE AND NUMBER OF BINS
     lambdaStart = lda_lower_range,
     lambdaEnd   = lda_upper_range,
     lambdaBins  = lda_bins,
     
     betaStart = bta_lower_range,
     betaEnd   = bta_upper_range,
     betaBins  = bta_bins
}
end


arg = { ... } -- -- TAKING USER INPUT
assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed") -- STILL EXPECTING SEED AS INPUT FOR THE FUTURE
argSeed = 34086709 -- -- SETTING SEED TO FIXED VALUE
prng = DSFMT.create(argSeed)

-- -- -- -- -- -- -- -- -- ROUNDING USER INPUT -- -- -- -- -- -- -- --
function round(num, places)
  local mult = 10.0^(places)
  return floor(num * mult + 0.5) / mult
end

-- -- -- -- -- -- ROUNDING TO AVOID DIFFERENT COMPUTER TERMINAL PRECISION -- -- -- -- -- --
dec = 9.0
evolveTime       = round( tonumber(arg[1]), dec )
rev_ratio        = round( tonumber(arg[2]), dec )
rscale_l         = round( tonumber(arg[3]), dec )
light_r_ratio    = round( tonumber(arg[4]), dec )
mass_l           = round( tonumber(arg[5]), dec )
light_mass_ratio = round( tonumber(arg[6]), dec )

-- -- -- -- -- -- -- -- -- DWARF PARAMETERS   -- -- -- -- -- -- -- --
revOrbTime = evolveTime
dwarfMass = mass_l / light_mass_ratio
rscale_t  = rscale_l / light_r_ratio
rscale_d  = rscale_t *  (1.0 - light_r_ratio)
mass_d    = dwarfMass * (1.0 - light_mass_ratio)
