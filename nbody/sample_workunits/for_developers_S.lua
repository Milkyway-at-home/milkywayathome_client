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
totalBodies           = 40000   -- -- NUMBER OF BODIES                                  -- --
nbodyLikelihoodMethod = "EMD"   -- -- HIST COMPARE METHOD                               -- --
nbodyMinVersion       = "1.86"  -- -- MINIMUM APP VERSION                               -- --

run_null_potential    = false   -- -- NULL POTENTIAL SWITCH                             -- --
use_tree_code         = true    -- -- USE TREE CODE NOT EXACT                           -- --
print_reverse_orbit   = false   -- -- PRINT REVERSE ORBIT SWITCH                        -- --
print_out_parameters  = false   -- -- PRINT OUT ALL PARAMETERS                          -- --

LMC_body              = true    -- -- PRESENCE OF LMC                                   -- --
LMC_scaleRadius       = 8.76 --15
LMC_Mass              = 218634.8191    --449865.888
LMC_DynamicalFriction = true    -- -- LMC DYNAMICAL FRICTION SWITCH (IGNORED IF NO LMC) -- --
CoulombLogarithm      = 0.470003629 -- -- (ln(1.6)) COULOMB LOGARITHM USED IN DYNAMICAL FRACTION CALCULATION -- --

SunGCDist             = 8.0       -- -- Distance between Sun and Galactic Center -- --

UseOldSofteningLength = 1         -- -- Uses old softening length formula from v1.76 and eariler -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 



-- -- -- -- MULTIPLE INPUT SWITCH -- -- -- --
n=2
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- MODEL SETTINGS -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- ModelComponent Options: 
-- --       2 - TWO COMPONENT MODEL     -- -- -- -- -- -- -- -- -- -- 
-- --       1 - SINGLE COMPONENT MODEL  -- -- -- -- -- -- -- -- -- -- 
-- --       0 - NO DWARF MODEL          -- -- -- -- -- -- -- -- -- -- 
---ModelComponents   = 1         -- -- TWO COMPONENTS SWITCH      -- --
componentList = {1,1}
-- componentList should be a table (list) with each entry being either 1 or 2.
-- 1 indicates a single component model, 2 indicates a two-component model.
if #componentList ~= n then
    error("Error: The length of the component list does not match the value of n.")
end

local hasNonZeroComponent = false
for _, component in ipairs(componentList) do
    if component ~= 0 then
        hasNonZeroComponent = true
        break
    end
end


manual_bodies     = false     -- -- USE THE MANUAL BODY LIST   -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 


-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- PARAMETER SETTINGS   -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- -- -- -- -- -- -- -- -- HISTOGRAM   -- -- -- -- -- -- -- -- -- -- -- -- --
Output_LB_coord = false    -- include Lambda-Beta coordinates in output file

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
use_best_likelihood  = false    -- use the best likelihood return code (ONLY SET TO TRUE FOR RUN-COMPARE)
best_like_start      = 0.98    -- what percent of sim to start

use_beta_disps       = true    -- use beta dispersions in likelihood
use_vel_disps        = false    -- use velocity dispersions in likelihood

-- if one of these is true, will get output for all 3 of the new histograms
-- if not computing likelihood scores, still need one of these to be true if want them computed/output
use_beta_comp        = true  -- calculate average beta, use in likelihood
use_vlos_comp        = true  -- calculate average los velocity, use in likelihood
use_avg_dist         = true  -- calculate average distance, use in likelihood

-- number of additional forward evolutions to do to calibrate the rotation of the bar
-- numCalibrationRuns + 1 additional forward evolutions will be done
-- if no bar potential is being used, this variable will be ignored
numCalibrationRuns = 0
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- ADVANCED DEVELOPER OPTIONS -- -- -- -- -- -- -- --        
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- These options only work if you compile nbody with  -- -- --
-- -- -- -- -- -- the -DNBODY_DEV_OPTIONS set to on                  -- -- --   

useMultiOutputs       = true       -- -- WRITE MULTIPLE OUTPUTS       -- --
freqOfOutputs         = 1000         -- -- FREQUENCY OF WRITING OUTPUTS -- --

timestep_control      = true       -- -- control number of steps      -- --
Ntime_steps           = 3000        -- -- number of timesteps to run   -- --

use_max_soft_par      = false       -- -- limit the softening parameter value to a max value
max_soft_par          = 0.8         -- -- kpc, if switch above is turned on, use this as the max softening parameter
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        




arg = { ... } -- -- TAKING USER INPUT
assert(#arg >= 6, "Expects either 6 or 12 arguments, and optional manual body list")
assert(argSeed ~= nil, "Expected seed") -- STILL EXPECTING SEED AS INPUT FOR THE FUTURE
argSeed = 7854614814 -- -- SETTING SEED TO FIXED VALUE
--argSeed = 34086710 -- -- SETTING SEED TO FIXED VALUE
prng = DSFMT.create(argSeed)

-- -- -- -- -- -- -- -- -- ROUNDING USER INPUT -- -- -- -- -- -- -- --
function round(num, places)
  local mult = 10.0^(places)
  return floor(num * mult + 0.5) / mult
end

-- -- -- -- -- -- ROUNDING TO AVOID DIFFERENT COMPUTER TERMINAL PRECISION -- -- -- -- -- --
dec = 9.0
evolveTime       = round( 3.63330 , dec )    -- Forward Time (Gyrs)
time_ratio       = round( 1, dec )    -- Forward Time / Backward Time

-- rscale_l         = {round( 2.89, dec )}    -- Baryonic Radius (kpc)
-- light_r_ratio    = {round( 0.2, dec )}    -- Baryonic Radius / (Baryonic Radius + Dark Matter Radius)
-- mass_l           = {round( 6298.125, dec )}    -- Baryonic Mass (Structure Mass Units)
-- light_mass_ratio = {round( 0.4, dec )}    -- Baryonic Mass / (Baryonic Mass + Dark Matter Mass)
-- orbit_parameter_l   = {round( 302.801, dec )}
-- orbit_parameter_b   = {round( -44.328, dec )}
-- orbit_parameter_r   = {round( 62.4, dec )}
-- orbit_parameter_vx  = {round( 21.99, dec )}
-- orbit_parameter_vy  = {round( -201.36, dec )}
-- orbit_parameter_vz  = {round( 171.25, dec )}

-- File with Individual Particles (.out file)
-- Dwarf series: 1.SMC 2.Sag 3.Fornax 4.LeoI 5.Sculptor 6.LeoII 7.Sextans 8.Carina 9.Draco 10.Umi 11.CvnI
-- Orphan Sagittarius

-- rscale_l         = {round( 1.0,dec)}    -- Baryonic Radius (kpc)
-- light_r_ratio    = {round( 0.2, dec )}    -- Baryonic Radius / (Baryonic Radius + Dark Matter Radius)
-- mass_l           = {round(  44.98658882 ,dec)}    -- Baryonic Mass (Structure Mass Units)
-- light_mass_ratio = {round( 0.05,dec)}    -- Baryonic Mass / (Baryonic Mass + Dark Matter Mass)
-- orbit_parameter_l   = {round( 299 ,dec)}
-- orbit_parameter_b   = {round( 5.75 ,dec)}
-- orbit_parameter_r   = {round( 17.8,dec)}
-- orbit_parameter_vx  = {round(  -237.4,dec)}
-- orbit_parameter_vy  = {round(  4.4,dec)}
-- orbit_parameter_vz  = {round( 233.1 ,dec)}

rscale_l         = {round( 1.5895,dec),round( 1.0, dec )}    -- Baryonic Radius (kpc) 2.89,1.445, 1.5895, 1.734, 2.023, 
light_r_ratio    = {round( 0.2, dec ),round( 0.2,dec)}    -- Baryonic Radius / (Baryonic Radius + Dark Matter Radius)
mass_l           = {round(  6298.125,dec),round( 44.98658882, dec )}    -- Baryonic Mass (Structure Mass Units) 6298.125 ，
light_mass_ratio = {round( 0.4,dec),round( 0.05, dec )}    -- Baryonic Mass / (Baryonic Mass + Dark Matter Mass)
orbit_parameter_l   = {round( 302.801 ,dec),round( 299, dec )}
orbit_parameter_b   = {round( -44.328 ,dec),round( 5.75, dec )}
orbit_parameter_r   = {round( 62.4,dec),round( 17.8, dec )}
orbit_parameter_vx  = {round(  21.99,dec),round( -237.4, dec )}
orbit_parameter_vy  = {round(  -201.36,dec),round(4.4, dec )}
orbit_parameter_vz  = {round( 171.25 ,dec),round( 233.1, dec )}

-- -- -- -- -- -- -- -- -- DWARF STARTING LOCATION   -- -- -- -- -- -- -- --
-- these only get used if only 6 parameters are input from shell script
-- otherwise they get reset later with the inputs (if 11 given)
--[[
preset_orbit_parameter_l  = 258
preset_orbit_parameter_b  = 45.8
preset_orbit_parameter_r  = 21.5
preset_orbit_parameter_vx = -185.5
preset_orbit_parameter_vy = 54.7
preset_orbit_parameter_vz = 147.4
]]
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        
-- -- -- -- -- -- -- -- -- CHECK TIMESTEPS -- -- -- -- -- -- -- -- 
-- manual_body_file = arg[13]
TooManyTimesteps = 0
        
function makePotential()
   if(run_null_potential == true) then
       print("running in null potential")
       return nil
   else
        --NOTE: To exclude a component from the potential, set component to "<component_name>.none" and include only an arbitrary "mass" argument
        -- return  Potential.create{
        --     spherical = Spherical.hernquist{ mass  = 1.52954402e5, scale = 0.7 },
        --     disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
        --     disk2     = Disk.none{ mass = 3.0e5 },
        --     halo      = Halo.logarithmic{ vhalo = 74.61, scaleLength = 12.0, flattenZ = 1.0 }
        -- }--vhalo = 74.61 kpc/gy = 73 km/s
        return  Potential.create{
            spherical = Spherical.hernquist{ mass  = 20243.9650, scale = 0.442 },
            disk  = Disk.miyamotoNagai{ mass = 305908.804, scaleLength = 3.0, scaleHeight = 0.28 },
            disk2 = Disk.none{ mass = 0.0 },
            halo  = Halo.nfwmass{ scaleLength = 16.0, mass = 1.96591393e6 }
            }--vhalo = 74.61 kpc/gy = 73 km/s
   end
end

function get_timestep()
    if(timestep_control) then
        t = (evolveTime) / (Ntime_steps)
    elseif(ModelComponents == 2) then--disable now for multidwarfs

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

    if ((evolveTime/t > 150000 or t ~= t) and not timestep_control) then
        TooManyTimesteps = 1
        t = evolveTime/4.0
    end

    return t
end


function get_soft_par()
    --softening parameter only calculated based on dwarf,
    --so if manual bodies is turned on the calculated s.p. may be too large
    sp = calculateEps2(totalBodies, rscale_l[1], rscale_d[1], mass_l[1], mass_d[1], UseOldSofteningLength)

    if ((manual_bodies or use_max_soft_par) and (sp > max_soft_par^2)) then --dealing with softening parameter squared
        print("Using maximum softening parameter value of " .. tostring(max_soft_par) .. " kpc")
        return max_soft_par^2
    else
        return sp
    end
end


function makeContext()
   return NBodyCtx.create{
      dwarfn      = n,
      timeEvolve  = evolveTime,
      timeBack    = revOrbTime,
      timestep    = get_timestep(),
      eps2        = get_soft_par(),
      b           = orbit_parameter_b,
      r           = orbit_parameter_r,
      vx          = orbit_parameter_vx,
      vy          = orbit_parameter_vy,
      vz          = orbit_parameter_vz,
      sunGCDist   = SunGCDist,
      criterion   = criterion,
      OutputLB    = Output_LB_coord,
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
      LMC           = LMC_body,
      LMCmass       = LMC_Mass,
      LMCscale      = LMC_scaleRadius,
      LMCDynaFric   = LMC_DynamicalFriction,
      coulomb_log   = CoulombLogarithm,
      calibrationRuns = numCalibrationRuns
   }
end



function makeBodies(ctx, potential)
  local firstModel = {}
  local finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity = {}, {}
  --Setting finalPosition, finalVelocity as empty list, LMC value will be nil
    if TooManyTimesteps == 1 then
        totalBodies = 1
    end

    if(run_null_potential == true and manual_bodies == true) then
        for i = 1, n do
            table.insert(finalPosition, lbrToCartesian(ctx, Vector.create(orbit_parameter_l[i], orbit_parameter_b[i], orbit_parameter_r[i])))
            table.insert(finalVelocity, Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz))
        end
    elseif(run_null_potential == true) then
        print("placing dwarf at origin")
        for i = 1, n do
            table.insert(finalPosition, Vector.create(0, 0, 0))
            table.insert(finalVelocity, Vector.create(0, 0, 0))
        end
        print(finalPosition[1])
    else 
    	if (LMC_body) then
            -- Old Single Body function with LMC and Method:

    		-- finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity = reverseOrbit_LMC{
	        --     potential   = potential,
	        --     position    = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
	        --     velocity    = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
	        --     LMCposition = Vector.create(-1.1, -41.1, -27.9),
	        --     LMCvelocity = Vector.create(-57, -226, 221), 
            --         LMCmass     = LMC_Mass,
            --         LMCscale    = LMC_scaleRadius,
            --         LMCDynaFric = LMC_DynamicalFriction,
            --         coulomb_log = CoulombLogarithm,
            --         ftime       = evolveTime,
	        --     tstop       = revOrbTime,
	        --     dt          = ctx.timestep / 10.0
	        --     }

            local potential = potential
            local position = lbrToCartesianTable(ctx, Vector.creates(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r))
            local velocity = Vector.creates(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz)
            local LMCposition = Vector.create(-0.52, -40.8, -26.5)
            local LMCvelocity = Vector.create(-58.2, -231, 226)
            local LMCmass = LMC_Mass
            local LMCscale = LMC_scaleRadius
            local LMCDynaFric = LMC_DynamicalFriction and 1 or 0
            local coulomb_log = CoulombLogarithm
            local ftime = evolveTime
            local tstop = revOrbTime
            local dt = ctx.timestep / 10.0      
            local masses    = dwarfMass  
            local rscales = rscale_t

            finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity = reverseOrbitS_LMC(potential, position, velocity, LMCposition, LMCvelocity, LMCmass, LMCscale, LMCDynaFric, coulomb_log, ftime, tstop, dt, masses, rscales)      
	    else
            local potential = potential
            local position  = lbrToCartesianTable(ctx, Vector.creates(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r))
            local velocity  = Vector.creates(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz)
            local tstop     = revOrbTime
            local dt        = ctx.timestep / 10.0 
            local masses    = dwarfMass

            finalPosition, finalVelocity = reverseOrbitS(potential, position, velocity, tstop, dt, masses)
            
            -- Old Single Body function and Method:
            
            -- finalPosition, finalVelocity = reverseOrbit{
	        --     potential = potential,
	        --     position  = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
	        --     velocity  = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
	        --     tstop     = revOrbTime,
	        --     dt        = ctx.timestep / 10.0
	        --     }
         end
    end

    if(print_reverse_orbit == true) then
        local placeholderPos, placeholderVel = {}, {}
            for i = 1, n do
                local phPos, phVel = PrintReverseOrbit{
                potential = potential,
                position  = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
                velocity  = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
                tstop     = .14,
                tstopf    = .20,
                dt        = ctx.timestep / 10.0
                }
            table.insert(placeholderPos, phPos)
            table.insert(placeholderVel, phVel)
            end
        print('Printing reverse orbit')
    end

    -- if(ModelComponents == 2) then 
    --     for i = 1, n do
    --         local Model = predefinedModels.mixeddwarf{
    --             nbody       = totalBodies,
    --             prng        = prng,
    --             position    = finalPosition[i],
    --             velocity    = finalVelocity[i],
    --             comp1       = Dwarf.plummer{mass = mass_l[i], scaleLength = rscale_l[i]}, -- Dwarf Options: plummer, nfw, general_hernquist
    --             comp2       = Dwarf.plummer{mass = mass_d[i], scaleLength = rscale_d[i]}, -- Dwarf Options: plummer, nfw, general_hernquist
    --             ignore      = true
    --             }
    --         for _, row in ipairs(Model) do
    --             table.insert(firstModel, row)
    --         end
    --         print(string.format("Dwarf %d bodies generation finished", i))
    --     end
        
    -- elseif(ModelComponents == 1) then
    --     for i = 1, n do
    --         local Model = predefinedModels.plummer{
    --             nbody       = totalBodies,
    --             prng        = prng,
    --             position    = finalPosition[i],
    --             velocity    = finalVelocity[i],
    --             mass        = mass_l[i],
    --             scaleRadius = rscale_l[i],
    --             ignore      = false
    --         }
    --         for _, row in ipairs(Model) do
    --             table.insert(firstModel, row)
    --         end
    --         print(string.format("Dwarf %d bodies generation finished", i))
    --     end
    -- end
    for i = 1, n do
        local componentType = componentList[i]
        
        if componentType ~= 1 and componentType ~= 2 and componentType ~= 0 then
            error(string.format("Error: Invalid component type %d at index %d. Only 1 or 2 are allowed.", componentType, i))
        end
    
        if componentType == 2 then 
            local Model = predefinedModels.mixeddwarf{
                nbody       = totalBodies,
                prng        = prng,
                position    = finalPosition[i],
                velocity    = finalVelocity[i],
                comp1       = Dwarf.plummer{mass = mass_l[i], scaleLength = rscale_l[i]}, -- Dwarf Options: plummer, nfw, general_hernquist
                comp2       = Dwarf.plummer{mass = mass_d[i], scaleLength = rscale_d[i]}, -- Dwarf Options: plummer, nfw, general_hernquist
                ignore      = true
            }
            for _, row in ipairs(Model) do
                table.insert(firstModel, row)
            end
            print(string.format("Dwarf %d bodies generation finished with two components", i))
    
        elseif componentType == 1 then
            local Model = predefinedModels.plummer{
                nbody       = totalBodies,
                prng        = prng,
                position    = finalPosition[i],
                velocity    = finalVelocity[i],
                mass        = mass_l[i],
                scaleRadius = rscale_l[i],
                ignore      = false
            }
            for _, row in ipairs(Model) do
                table.insert(firstModel, row)
            end
            print(string.format("Dwarf %d bodies generation finished with a single component", i))
        end
    end
    if(manual_bodies) then
        manualModel = predefinedModels.manual_bodies{
        body_file   = manual_body_file,
    }
         
    end
    if hasNonZeroComponent and manual_bodies then
        print("FirstModel & ManualModel")
        return firstModel, manualModel
    elseif hasNonZeroComponent and not manual_bodies then
        print("FirstModel")
        return firstModel
    elseif not hasNonZeroComponent and manual_bodies then        
        print("ManualModel")
        return manualModel
    else
        error("Don't you want to simulate something?")
    end    
    print("finished makebodies")
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


-- -- -- -- -- -- -- -- -- DWARF PARAMETERS   -- -- -- -- -- -- -- --
revOrbTime = evolveTime / time_ratio
if use_best_likelihood then
    evolveTime = (2.0 - best_like_start) * evolveTime --making it evolve slightly longer
    eff_best_like_start = best_like_start / (2.0 - best_like_start)
else
    eff_best_like_start = best_like_start
end


dwarfMass = {}
rscale_t = {}
rscale_d = {}
mass_d = {}
if(ModelComponents == 1) then
    for i = 1, n do
        dwarfMass[i]  = mass_l[i]
        rscale_t[i]   = rscale_l[i]
        rscale_d[i]  = 1.0
        mass_d[i]     = 0.0
    end
else    
    for i = 1, n do
        dwarfMass[i] = mass_l[i] / light_mass_ratio[i]
        rscale_t[i]  = rscale_l[i] / light_r_ratio[i]
        rscale_d[i]  = rscale_t[i] * (1.0 - light_r_ratio[i])
        mass_d[i]    = dwarfMass[i] * (1.0 - light_mass_ratio[i])
    end
end
   

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
    print('forward time=', evolveTime, '\nreverse time=',  revOrbTime)
    print('mass_l sim=', mass_l, '\nmass_d sim=', mass_d)
    print('light mass solar=', mass_l * 222288.47, '\ndark mass solar=', mass_d * 222288.47)
    print('total mass solar= ', (mass_d + mass_l) * 222288.47)
    print('rl = ', rscale_l, 'rd = ', rscale_d)
end
