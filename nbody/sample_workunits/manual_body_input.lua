-- /* Copyright (c) 2016 Siddhartha Shelton */

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- DEAR LUA USER:
-- This is the developer version of the lua parameter file. 
-- It gives all the options you can have. 
-- Many of these the client will not need.

-- NOTE --
-- if you are using single component plummer model, it will take the baryonic
-- matter component parameters. meaning you input should look like
-- ft, bt, rscale_baryon, radius_ratio, baryon mass, mass ratio
-- typical parameters: 4.0, 1.0, 0.2, 0.2, 12, 0.2
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        
        
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- STANDARD  SETTINGS   -- -- -- -- -- -- -- -- -- --        
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
nbodyLikelihoodMethod = "EMD"   -- -- HIST COMPARE METHOD        -- --
nbodyMinVersion       = "1.66"  -- -- MINIMUM APP VERSION        -- --

run_null_potential    = true   -- -- NULL POTENTIAL SWITCH      -- --
use_tree_code         = false   -- -- USE TREE CODE NOT EXACT    -- --
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
use_best_likelihood  = false    -- use the best likelihood return code
best_like_start      = 0.98    -- what percent of sim to start
use_vel_disps        = true    -- use velocity dispersions in likelihood
        
timestep_control     = false -- -- control number of steps    -- --
Ntime_steps          = 0    -- -- number of timesteps to run -- --

-- -- -- -- -- -- -- -- -- DWARF STARTING LOCATION   -- -- -- -- -- -- -- --
orbit_parameter_l  = 218
orbit_parameter_b  = 53.5
orbit_parameter_r  = 28.6
orbit_parameter_vx = -156 
orbit_parameter_vy = 79 
orbit_parameter_vz = 107
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        

        
function makePotential()
   if(run_null_potential == true) then
       print("running in null potential")
       return nil
   else
        return  Potential.create{
            spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
            disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
            halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
        }
   end
end


function get_timestep()
    -- this is a hard coded value. We may need to revist this at a later date.
    t = 1e-5
--     print(t)
    return t
end


function makeContext()
   soften_length  = 2e-5
   
   return NBodyCtx.create{
      timeEvolve  = evolveTime,
      timestep    = get_timestep(),
      eps2        = soften_length,
      criterion   = criterion,
      useQuad     = true,
      useBestLike = use_best_likelihood,
      useVelDisp  = use_vel_disps,
      BestLikeStart = best_like_start,
      Nstep_control = timestep_control,
      Ntsteps       = Ntime_steps,
      theta       = 1.0
   }
end

-- soften_length = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
-- print('soften_length ', calculateEps2(totalBodies, soften_length ))

function makeBodies(ctx, potential)
  local firstModel
  
    firstModel = predefinedModels.manual_bodies{
        body_file   = file,
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
assert(#arg == 2, "Expected 2 arguments")
assert(argSeed ~= nil, "Expected seed") -- STILL EXPECTING SEED AS INPUT FOR THE FUTURE
-- argSeed = 34086709 -- -- SETTING SEED TO FIXED VALUE
argSeed = 7854614814 -- -- SETTING SEED TO FIXED VALUE
prng = DSFMT.create(argSeed)

evolveTime       = tonumber(arg[1])
file             = arg[2]

-- print(evolveTime, file)

if(use_tree_code) then
    criterion = "TreeCode"
--     criterion = "NewCriterion"
else
    criterion = "Exact"
end

