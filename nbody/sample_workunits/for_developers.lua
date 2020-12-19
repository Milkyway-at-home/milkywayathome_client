-- /* Copyright (c) 2016-2018 Siddhartha Shelton */

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- DEAR LUA USER:
-- This is the developer version of the lua parameter file. 
-- It gives all the options you can have. 
-- Many of these the client will not need.

-- NOTE --
-- to fully utilize this lua, need to compile with -DNBODY_DEV_OPTIONS=ON
-- if you are using single component plummer model, it will take the baryonic
-- matter component parameters. meaning you input should look like
-- ft, time_ratio, rscale_baryon, radius_ratio, baryon mass, mass ratio
-- typical parameters: 4.0, 1.0, 0.2, 0.2, 12, 0.2 (52.5, 28.6, -156, 79, 107)

-- available option: using a user inputted list of bodies. Sent in as an 
-- optional arguement after dwarf parameter list
-- MUST still include dwarf parameter list
-- can control what model to use below
-- simulation time still taken as the first parameter in the list
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        
        
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- STANDARD  SETTINGS   -- -- -- -- -- -- -- -- -- --        
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
totalBodies           = 20000   -- -- NUMBER OF BODIES           -- --
nbodyLikelihoodMethod = "EMD"   -- -- HIST COMPARE METHOD        -- --
nbodyMinVersion       = "1.76"  -- -- MINIMUM APP VERSION        -- --

run_null_potential    = false   -- -- NULL POTENTIAL SWITCH      -- --
use_tree_code         = true    -- -- USE TREE CODE NOT EXACT    -- --
print_reverse_orbit   = false   -- -- PRINT REVERSE ORBIT SWITCH -- --
print_out_parameters  = false   -- -- PRINT OUT ALL PARAMETERS   -- --

LMC_body              = true    -- -- PRESENCE OF LMC            -- --
LMCtotalBodies        = 2000    -- -- DON'T SET TO VALUES SMALLER THAN 100 
LMC_scaleRadius       = 15
LMC_Mass              = 449865.888
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 



-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- MODEL SETTINGS -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- ModelComponent Options: 
-- --       2 - TWO COMPONENT MODEL     -- -- -- -- -- -- -- -- -- -- 
-- --       1 - SINGLE COMPONENT MODEL  -- -- -- -- -- -- -- -- -- -- 
-- --       0 - NO DWARF MODEL          -- -- -- -- -- -- -- -- -- -- 
ModelComponents   = 2         -- -- TWO COMPONENTS SWITCH      -- --
manual_bodies     = false     -- -- USE THE MANUAL BODY LIST   -- --
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

SigmaCutoff          = 2.5     -- -- sigma cutoff for outlier rejection DO NOT CHANGE -- --
SigmaIter            = 6       -- -- number of times to apply outlier rejection DO NOT CHANGE -- --
Correction           = 1.111   -- -- correction for outlier rejection   DO NOT CHANGE -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

-- -- -- -- -- -- -- -- -- AlGORITHM OPTIONS -- -- -- -- -- -- -- --
use_best_likelihood  = true    -- use the best likelihood return code
best_like_start      = 0.98    -- what percent of sim to start

use_beta_disps       = true    -- use beta dispersions in likelihood
use_vel_disps        = false    -- use velocity dispersions in likelihood

-- if one of these is true, will get output for all 3 of the new histograms
-- if not computing likelihood scores, still need one of these to be true if want them computed/output
use_beta_comp        = false  -- calculate average beta, use in likelihood
use_vlos_comp        = false  -- calculate average los velocity, use in likelihood
use_avg_dist         = false  -- calculate average distance, use in likelihood
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- ADVANCED DEVELOPER OPTIONS -- -- -- -- -- -- -- --        
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- These options only work if you compile nbody with  -- -- --
-- -- -- -- -- -- the -DNBODY_DEV_OPTIONS set to on                  -- -- --   

useMultiOutputs       = false        -- -- WRITE MULTIPLE OUTPUTS       -- --
freqOfOutputs         = 6            -- -- FREQUENCY OF WRITING OUTPUTS -- --

timestep_control     = false         -- -- control number of steps      -- --
Ntime_steps          = 10            -- -- number of timesteps to run   -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        




-- -- -- -- -- -- -- -- -- DWARF STARTING LOCATION   -- -- -- -- -- -- -- --
orbit_parameter_l  = 218

-- these only get used if only 6 parameters are input from shell script
-- otherwise they get reset later with the inputs (if 11 given)
orbit_parameter_b  = 53.5
orbit_parameter_r  = 28.6
orbit_parameter_vx = -156 
orbit_parameter_vy = 79 
orbit_parameter_vz = 107
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        
-- -- -- -- -- -- -- -- -- CHECK TIMESTEPS -- -- -- -- -- -- -- -- 
TooManyTimesteps = 0
        
function makePotential()
   if(run_null_potential == true) then
       print("running in null potential")
       return nil
   else
        --NOTE: To exculde a component from the potential, set component to "<component_name>.none" and include only an arbitrary "mass" argument
        return  Potential.create{
            spherical = Spherical.hernquist{ mass  = 1.52954402e5, scale = 0.7 },
            disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
            disk2     = Disk.orbitingBar{ mass = 2.429275796e4, scaleLength = 5.4, patternSpeed = 39, startAngle = 0.488692},
            halo      = Halo.logarithmic{ vhalo = 74.61, scaleLength = 12.0, flattenZ = 1.0 }
        }--vhalo = 74.61 kpc/gy = 73 km/s
   end
end

function get_timestep()
    if(timestep_control) then
      t = (evolveTime) / (Ntime_step)
    elseif(ModelComponents == 2) then

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
        t = (1.0 / 100.0) * ( pi_4_3 * s)^(1.0/2.0)
        
    else 
        t = sqr(1.0 / 10.0) * sqrt((pi_4_3 * cube(rscale_l)) / (mass_l))
    end

    if (evolveTime/t > 150000 or t ~= t) then
        TooManyTimesteps = 1
        t = evolveTime/4.0
    end

    return t
end


function makeContext()
   soften_length  = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
   return NBodyCtx.create{
      timeEvolve  = evolveTime,
      timeBack    = revOrbTime,
      timestep    = get_timestep(),
      eps2        = calculateEps2(totalBodies, soften_length),
      b           = orbit_parameter_b,
      r           = orbit_parameter_r,
      vx          = orbit_parameter_vx,
      vy          = orbit_parameter_vy,
      vz          = orbit_parameter_vz,
      criterion   = criterion,
      useQuad     = true,
      useBestLike   = use_best_likelihood,
      BestLikeStart = eff_best_like_start,
      useVelDisp    = use_vel_disps,
      useBetaDisp   = use_beta_disps,
      useBetaComp   = use_beta_comp,
      useVlos       = use_vlos_comp,
      useDist       = use_avg_dist,
      Nstep_control = timestep_control,
      Ntsteps       = Ntime_steps,
      BetaSigma     = SigmaCutoff,
      VelSigma      = SigmaCutoff,
      DistSigma     = SigmaCutoff,
      IterMax       = SigmaIter,
      BetaCorrect   = Correction,
      VelCorrect    = Correction,
      DistCorrect   = Correction,
      MultiOutput   = useMultiOutputs,
      OutputFreq    = freqOfOutputs,
      theta         = 1.0,
      LMC           = LMC_body
   }
end



function makeBodies(ctx, potential)
  local firstModel, LMCModel
  local finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity
    if TooManyTimesteps == 1 then
        totalBodies = 1
    end

    if(run_null_potential == true) then
        print("placing dwarf at origin")
        finalPosition, finalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)
        LMCfinalPosition, LMCfinalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)
    else 
    	if (LMC_body) then
    		finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity = reverseOrbit_LMC{
	            potential   = potential,
	            position    = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
	            velocity    = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
	            LMCposition = Vector.create(-1.1, -41.1, -27.9),
	            LMCvelocity = Vector.create(-57, -226, 221), 
              LMCmass     = LMC_Mass,
	            tstop       = revOrbTime,
	            dt          = ctx.timestep / 10.0
	            }

              
	    else
	        finalPosition, finalVelocity = reverseOrbit{
	            potential = potential,
	            position  = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
	            velocity  = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
	            tstop     = revOrbTime,
	            dt        = ctx.timestep / 10.0
	            }
         end
    end
    
    if(print_reverse_orbit == true) then
        local placeholderPos, placeholderVel = PrintReverseOrbit{
            potential = potential,
            position  = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
            velocity  = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
            tstop     = .14,
            tstopf    = .20,
            dt        = ctx.timestep / 10.0
        }
        print('Printing reverse orbit')
    end
    
  	if(LMC_body) then
  		LMCModel = predefinedModels.plummer{
            nbody       = LMCtotalBodies,
            prng        = prng,
            position    = LMCfinalPosition,
            velocity    = LMCfinalVelocity,
            mass        = LMC_Mass,
            scaleRadius = LMC_scaleRadius,
            ignore      = false
        }
    end

    if(ModelComponents == 2) then 
        firstModel = predefinedModels.mixeddwarf{
            nbody       = totalBodies,
            prng        = prng,
            position    = finalPosition,
            velocity    = finalVelocity,
            comp1       = Dwarf.plummer{mass = mass_l, scaleLength = rscale_l}, -- Dwarf Options: plummer, nfw, general_hernquist
            comp2       = Dwarf.plummer{mass = mass_d, scaleLength = rscale_d}, -- Dwarf Options: plummer, nfw, general_hernquist
            ignore      = true
        }
        
    elseif(ModelComponents == 1) then
        firstModel = predefinedModels.plummer{
            nbody       = totalBodies,
            prng        = prng,
            position    = finalPosition,
            velocity    = finalVelocity,
            mass        = mass_l,
            scaleRadius = rscale_l,
            ignore      = false
        }
  
    end
  
    if(manual_bodies) then
        manualModel = predefinedModels.manual_bodies{
        body_file   = manual_body_file,
    }
         
    end
    
    if(ModelComponents > 0 and manual_bodies) then 
        return firstModel, manualModel
    elseif(ModelComponents == 0 and manual_bodies) then
        return manualModel
    elseif(ModelComponents > 0 and not manual_bodies and LMC_body) then
        return firstModel, LMCModel
    elseif(ModelComponents > 0 and not manual_bodies) then
        return firstModel
    else    
        print("Don't you want to simulate something?")
    end
    
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
assert(#arg >= 6, "Expects either 6 or 11 arguments")
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
evolveTime       = round( tonumber(arg[1]), dec )    -- Forward Time
time_ratio       = round( tonumber(arg[2]), dec )    -- Forward Time / Backward Time
rscale_l         = round( tonumber(arg[3]), dec )    -- Baryonic Radius
light_r_ratio    = round( tonumber(arg[4]), dec )    -- Baryonic Radius / (Baryonic Radius + Dark Matter Radius)
mass_l           = round( tonumber(arg[5]), dec )    -- Baryonic Mass (Structure Mass Units)
light_mass_ratio = round( tonumber(arg[6]), dec )    -- Baryonic Mass / (Baryonic Mass + Dark Matter Mass)
if((#arg == 11)or(#arg == 12)) then
  orbit_parameter_b   = round( tonumber(arg[7]), dec )
  orbit_parameter_r   = round( tonumber(arg[8]), dec )
  orbit_parameter_vx  = round( tonumber(arg[9]), dec )
  orbit_parameter_vy  = round( tonumber(arg[10]), dec )
  orbit_parameter_vz  = round( tonumber(arg[11]), dec )
  manual_body_file = arg[12]
else
  manual_body_file = arg[7] -- File with Individual Particles (.out file)
end

-- -- -- -- -- -- -- -- -- DWARF PARAMETERS   -- -- -- -- -- -- -- --
revOrbTime = evolveTime / time_ratio
if use_best_likelihood then
    evolveTime = (2.0 - best_like_start) * evolveTime --making it evolve slightly longer
    eff_best_like_start = best_like_start / (2.0 - best_like_start)
else
    eff_best_like_start = best_like_start
end


dwarfMass = mass_l / light_mass_ratio
rscale_t  = rscale_l / light_r_ratio
rscale_d  = rscale_t *  (1.0 - light_r_ratio)
mass_d    = dwarfMass * (1.0 - light_mass_ratio)

if(manual_bodies and manual_body_file == nil) then 
    print 'WARNING: No body list given. Manual body input turn off'
    manual_bodies = false  --optional body list was not included
elseif(manual_bodies and ModelComponents == 0) then
    print 'Using user inputted body list only' 
    print( manual_body_file)
elseif(manual_bodies and ModelComponents ~= 0) then
    print 'Using dwarf model and user inputted body list'
end


if(use_tree_code) then
    criterion = "TreeCode"
else
    criterion = "Exact"
end

if(print_out_parameters) then
    print('forward time=', evolveTime, '\nrev time=',  revOrbTime)
    print('mass_l sim=', mass_l, '\nmass_d sim=', mass_d)
    print('light mass solar=', mass_l * 222288.47, '\ndark mass solar=', mass_d * 222288.47)
    print('total mass solar= ', (mass_d + mass_l) * 222288.47)
    print('rl = ', rscale_l, 'rd = ', rscale_d)
end
