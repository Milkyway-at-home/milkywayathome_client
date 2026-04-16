-- /* Copyright (c) 2016-2018 Siddhartha Shelton */

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
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
-- 222288.47 solar masses = 1 Structure Mass Unit (SMU)

-- available option: using a user inputted list of bodies. Sent in as an 
-- optional arguement after dwarf parameter list
-- MUST still include dwarf parameter list
-- can control what model to use below
-- simulation time still taken as the first parameter in the list

-- IMPORTANT -- IMPORTANT -- IMPORTANT -- IMPORTANT -- IMPORTANT -- 
-- Structural changes to this file also need to be changed in the 
-- lua files in the tests directory (nbody/tests/mixeddwarf_models/) and (nbody/tests/orphan_models/)
-- especially if the changes are not backwards compatible with the previous format
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- --  BASIC  SETTINGS  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --      
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
totalBodies           = 5000   -- -- NUMBER OF BODIES PER DWARF                                                -- --
totalLightBodies      = 2500   -- -- NUMBER OF LIGHT MATTER BODIES                                            -- --

nbodyLikelihoodMethod = "EMD"       -- -- HIST COMPARE METHOD                                                  -- --
nbodyMinVersion       = "1.95"      -- -- MINIMUM APP VERSION                                                  -- --

run_null_potential    = false       -- -- NULL POTENTIAL SWITCH                                                -- --
use_tree_code         = true        -- -- USE TREE CODE NOT EXACT                                              -- --
print_reverse_orbit   = false       -- -- PRINT REVERSE ORBIT SWITCH (WORKS FOR LMC_body = false)              -- --
print_out_parameters  = false       -- -- PRINT OUT ALL PARAMETERS                                             -- --

LMC_body              = false        -- -- PRESENCE OF LMC (TURN OFF FOR NULL POTENTIAL)                        -- --
LMC_function          = 1           -- -- 1: Plummer 2: Henrquist 3: Hernquist with cutoff                     -- --
LMC_scaleRadius       = 15          -- --  kpc                                                                 -- --
LMC_cutoff            = 16          -- --  kpc  This is used only for Hernquist with cutoff                    -- --
LMC_Mass       = 449865.888  -- -- SMU (used unless specified in arguments)                             -- --
LMC_DynamicalFriction = true    -- -- LMC DYNAMICAL FRICTION SWITCH (IGNORED IF NO LMC)                        -- --
CoulombLogarithm      = 15      -- -- ln(r/1.22*CoulombLogarithm) (Patel et al. 2020) COULOMB LOGARITHM USED   -- --
                                -- -- IN DYNAMICAL FRACTION CALCULATION                                        -- --

SunGCDist             = 8.0       -- -- Distance between Sun and Galactic Center                               -- --
SunVelx               = 10.3      -- -- Sun's x-velocity (kpc/Gyr) (Hogg et al. (2005))                        -- --
SunVely               = 229.2     -- -- Sun's y-velocity (kpc/Gyr)                                             -- --
SunVelz               = 6.9       -- -- Sun's z-velocity (kpc/Gyr)                                             -- --

UseOldSofteningLength = 0         -- -- Uses old softening length formula from v1.76 and eariler               -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- -- -- -- NOTE: USER INPUT AT RUNTIME IS CURRENTLY NOT FUNCTIONAL -- -- -- --
arg = { ... } -- -- TAKING USER INPUT
assert((#arg == 6 or #arg == 7 or #arg == 8 or #arg == 11 or #arg == 12 or #arg == 13), "Expects either 6, 7, 8, 11, 12, or 13 arguments")
assert(argSeed ~= nil, "Expected seed") -- STILL EXPECTING SEED AS INPUT FOR THE FUTURE
argSeed = 34086709 -- -- SETTING SEED TO FIXED VALUE
--argSeed = 34086710 -- -- SETTING SEED TO FIXED VALUE
prng = DSFMT.create(argSeed)

-- -- -- -- -- -- -- -- -- INPUT ROUNDING -- -- -- -- -- -- -- --
function round(num, places)
  local mult = 10.0^(places)
  return floor(num * mult + 0.5) / mult
end

dec = 9.0   -- -- number of decimals to round to (default: 9.0)

-- -- -- -- MISC. INPUTS -- -- -- --

n = 11                                  -- number of simulated dwarfs
evolveTime       = round( 3.0, dec )    -- Forward Time (Gyrs)
time_ratio       = round( 1, dec )      -- Forward Time / Backward Time

-- -- -- -- -- --  DWARF PARAMETER INPUTS  -- -- -- -- -- --
-- -- -- make sure arrays are of length n !!! -- -- -- -- --

-- vv will likely not apply if settings are changed
-- default index/name :   00  SMC          |   01  Sagittarius  |   02  Fornax       |   03  Leo I        |   04  Sculptor     |   05  Leo II       |   06  Sextans      |   07  Carina       |   08  Draco        |   09  Ursa Minor   |   10  C.Venatici I
rscale_l            = {round( 2.9,     dec),round( 1.53,    dec),round( 1.425,   dec),round( 0.43,    dec),round( 0.725,   dec),round( 0.96,    dec),round( 1.25,    dec),round( 0.465,   dec),round( 0.59,    dec),round( 0.42,    dec),round( 0.505,   dec)}  -- Baryonic Radius (kpc)
light_r_ratio       = {round( 0.2,     dec),round( 0.2,     dec),round( 0.2,     dec),round( 0.2,     dec),round( 0.2,     dec),round( 0.2,     dec),round( 0.2,     dec),round( 0.2,     dec),round( 0.2,     dec),round( 0.2,     dec),round( 0.2,     dec)}  -- Baryonic Radius / (Baryonic Radius + Dark Matter Radius)
mass_l              = {round( 2429.198,dec),round( 107.041, dec),round( 80.159,  dec),round( 20.968,  dec),round( 9.384,   dec),round( 2.918,   dec),round( 1.892,   dec),round( 1.647,   dec),round( 1.134,   dec),round( 0.899,   dec),round( 1.061,   dec)}  -- Baryonic Mass (Structure Mass Units)
light_mass_ratio    = {round( 0.0830,  dec),round( 0.0594,  dec),round( 0.1429,  dec),round( 0.0067,  dec),round( 0.0674,  dec),round( 0.0240,  dec),round( 0.0100,  dec),round( 0.0159,  dec),round( 0.0115,  dec),round( 0.0038,  dec),round( 0.0087,  dec)}  -- Baryonic Mass / (Baryonic Mass + Dark Matter Mass)
orbit_parameter_l   = {round( 302.801, dec),round( 5.569,   dec),round( 237.104, dec),round( 225.985, dec),round( 287.535, dec),round( 220.164, dec),round( 243.498, dec),round( 260.112, dec),round( 86.368,  dec),round( 104.9,   dec),round( 74.305,  dec)}  -- Galactocentric l
orbit_parameter_b   = {round( -44.328, dec),round( -14.166, dec),round( -65.651, dec),round( 49.112,  dec),round( -83.157, dec),round( 67.229,  dec),round( 42.272,  dec),round( -22.223, dec),round( 34.722,  dec),round( 44.8,    dec),round( 79.823,  dec)}  -- Galactocentric b
orbit_parameter_r   = {round( 62.4,    dec),round( 25,      dec),round( 143,     dec),round( 250,     dec),round( 88.91,   dec),round( 220,     dec),round( 90,      dec),round( 100,     dec),round( 80,      dec),round( 60,      dec),round( 220,     dec)}  -- Galactocentric r
orbit_parameter_vx  = {round( 21.99,   dec),round( 223.97,  dec),round( -27.04,  dec),round( 48.17,   dec),round( -22.11,  dec),round( 94.87,   dec),round( -194.39, dec),round( -28.48,  dec),round( -59.22,  dec),round( 19.12,   dec),round( 23.95,   dec)}  -- Galactocentric vx
orbit_parameter_vy  = {round( -201.36, dec),round( -5.34,   dec),round( -172.14, dec),round( -16.36,  dec),round( 197.28,  dec),round( 209.73,  dec),round( 30.33,   dec),round( -79.13,  dec),round( 60.33,   dec),round( 38.13,   dec),round( 47.45,   dec)}  -- Galactocentric vy
orbit_parameter_vz  = {round( 171.25,  dec),round( 185.78,  dec),round( 101.21,  dec),round( 254.15,  dec),round( -102.1,  dec),round( 114.61,  dec),round( 49.13,   dec),round( 164.44,  dec),round( -263.33, dec),round( -160.51, dec),round( 68.05,   dec)}  -- Galactocentric vz


-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- MODEL SETTINGS -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- --       ModelComponent Options:    -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- --       2 - TWO COMPONENT MODEL    -- -- -- -- -- -- -- -- -- -- -- -- --
-- --       1 - SINGLE COMPONENT MODEL  -- -- -- - -- -- -- -- -- -- -- -- -- 
-- --       0 - NO DWARF MODEL         -- -- -- -- -- -- -- -- -- -- -- -- --
ModelComponents   = 2         -- -- TWO COMPONENTS SWITCH   -- -- -- -- -- --
manual_bodies     = false     -- -- USE THE MANUAL BODY LIST   -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --



-- -- -- -- -- -- -- -- --  MANUAL INPUT CODE  -- -- -- -- -- -- -- -- --
-- -- -- requires more work to make operational with multidwarfs -- -- --
-- -- -- line-by-line prompt??  -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- make use as need be o7 -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- make use as need be o7 -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

-- evolveTime       = round( tonumber(arg[1]), dec )    -- Forward Time (Gyrs)
-- time_ratio       = round( tonumber(arg[2]), dec )    -- Forward Time / Backward Time
-- rscale_l         = round( tonumber(arg[3]), dec )    -- Baryonic Radius (kpc)
-- light_r_ratio    = round( tonumber(arg[4]), dec )    -- Baryonic Radius / (Baryonic Radius + Dark Matter Radius)
-- mass_l           = round( tonumber(arg[5]), dec )    -- Baryonic Mass (Structure Mass Units)
-- light_mass_ratio = round( tonumber(arg[6]), dec )    -- Baryonic Mass / (Baryonic Mass + Dark Matter Mass)
-- if (#arg == 7) then
--     if manual_bodies then
--         manual_body_file = arg[7]
--     else 
--         LMC_Mass = round( tonumber(arg[7]), dec )
--     end
-- elseif (#arg == 8) then
--     LMC_Mass = round( tonumber(arg[7]), dec )
--     manual_body_file = arg[8]
-- elseif (#arg == 12) then
--     orbit_parameter_l   = round( tonumber(arg[7]), dec )
--     orbit_parameter_b   = round( tonumber(arg[8]), dec )
--     orbit_parameter_r   = round( tonumber(arg[9]), dec )
--     orbit_parameter_vx  = round( tonumber(arg[10]), dec )
--     orbit_parameter_vy  = round( tonumber(arg[11]), dec )
--     orbit_parameter_vz  = round( tonumber(arg[12]), dec )
-- elseif (#arg == 13) then
--     orbit_parameter_l   = round( tonumber(arg[7]), dec )
--     orbit_parameter_b   = round( tonumber(arg[8]), dec )
--     orbit_parameter_r   = round( tonumber(arg[9]), dec )
--     orbit_parameter_vx  = round( tonumber(arg[10]), dec )
--     orbit_parameter_vy  = round( tonumber(arg[11]), dec )
--     orbit_parameter_vz  = round( tonumber(arg[12]), dec )
--     if manual_bodies then
--         manual_body_file = arg[13]
--     else
--         LMC_Mass = round( tonumber(arg[13]), dec )
--     end
-- elseif (#arg == 14) then
--     orbit_parameter_l   = round( tonumber(arg[7]), dec )
--     orbit_parameter_b   = round( tonumber(arg[8]), dec )
--     orbit_parameter_r   = round( tonumber(arg[9]), dec )
--     orbit_parameter_vx  = round( tonumber(arg[10]), dec )
--     orbit_parameter_vy  = round( tonumber(arg[11]), dec )
--     orbit_parameter_vz  = round( tonumber(arg[12]), dec )
--     LMC_Mass = round( tonumber(arg[13]), dec )
--     manual_body_file = arg[14]
-- else
--     -- fallback to preset orbit parameters and LMC mass if not enough args
--     orbit_parameter_l   = preset_orbit_parameter_l
--     orbit_parameter_b   = preset_orbit_parameter_b
--     orbit_parameter_r   = preset_orbit_parameter_r
--     orbit_parameter_vx  = preset_orbit_parameter_vx
--     orbit_parameter_vy  = preset_orbit_parameter_vy
--     orbit_parameter_vz  = preset_orbit_parameter_vz
--     LMC_Mass = preset_LMC_Mass
-- end


-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- PARAMETER SETTINGS   -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- -- -- -- -- -- -- --  OUTPUT SETTINGS  -- -- -- -- -- -- -- -- -- -- -- --
generateSimpleOutput = true       -- Simple output file includes: x, y, z, vx, vy, vz, mass
-- Full output file includes: x, y, z, l, b, r, vx, vy, vz, mass, vlos, pmra, pmdec, [lambda, beta]
-- NOTE: Lambda and Beta are optional and will only be included if the histogram parameters are set in makeHistogram()

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
use_best_likelihood  = false    -- use the best likelihood return code (ONLY SET TO TRUE FOR RUN-COMPARE)
best_like_start      = 0.98    -- what percent of sim to start

use_beta_disps       = true    -- use beta dispersions in likelihood
use_vel_disps        = false    -- use velocity dispersions in likelihood

-- if one of these is true, will get output for all 3 of the new histograms
-- if not computing likelihood scores, still need one of these to be true if want them computed/output
use_beta_comp        = true  -- calculate average beta, use in likelihood
use_vlos_comp        = true  -- calculate average los velocity, use in likelihood
use_avg_dist         = true  -- calculate average distance, use in likelihood
use_pm_comp          = true  -- calculate proper motion, use in likelihood

-- if using momentum likelihood, include momentum information in the parameters of the input
-- histogram (after <histogram> )with the following lines:
    -- L = {Lx, Ly, Lz} (angular momentum vector)
    -- LErr = {Err_Lx, Err_Ly, Err_Lz} (uncertainty in angular momentum vector)
-- These are in units of kpc^2/Gyr (no mass included)
use_momentum         = true  -- calculate angular momentum, use in likelihood

-- number of additional forward evolutions to do to calibrate the rotation of the bar
-- numCalibrationRuns + 1 additional forward evolutions will be done
-- if no bar potential is being used, this variable will be ignored
numCalibrationRuns = 0
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- ADVANCED DEVELOPER OPTIONS -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- These options only work if you compile nbody with  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- the -DNBODY_DEV_OPTIONS set to on -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- - -- -- -- -- -- -- --  

useMultiOutputs       = true     -- -- WRITE MULTIPLE OUTPUTS                                                            -- --
freqOfOutputs         = 30         -- -- FREQUENCY OF WRITING OUTPUTS                                                     -- --

timestep_control      = true       -- -- control number of steps                                                          -- --
Ntime_steps           = 3000        -- -- number of timesteps to run                                                       -- --

use_max_soft_par      = false       -- -- limit the softening parameter value to a max value                               -- --
max_soft_par          = 0.8         -- -- kpc, if switch above is turned on, use this as the max softening parameter       -- --

generateInitialOutput = true       -- -- save initial dwarf galaxy state to initial.out before evolution                   -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- -- -- -- -- -- -- -- -- CHECK TIMESTEPS -- -- -- -- -- -- -- -- 
TooManyTimesteps = 0
        
function makePotential()
   if(run_null_potential == true) then
       print("running in null potential")
       return nil
   else
        --NOTE: To exclude a component from the potential, set component to "<component_name>.none" and include only an arbitrary "mass" argument
        return  Potential.create{
            spherical = Spherical.hernquist{ mass  = 20243.9650, scale = 0.442 },
            disk      = Disk.miyamotoNagai{ mass = 305908.804, scaleLength = 3.0, scaleHeight = 0.28 },
            disk2     = Disk.none{ mass = 0.0 },
            halo      = Halo.nfwmass{ scaleLength = 16.0, mass = 1.96591393e6 }
        }
   end
end

function get_timestep()
    if(timestep_control) then
        t = (evolveTime) / (Ntime_steps)
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

    if ((evolveTime/t > 150000 or t ~= t) and not timestep_control) then
        -- We could throw an error here, but instead let it run fast and return a poor likelihood
        -- This way users won't see errors in their workunit logs
        TooManyTimesteps = 1
        t = evolveTime/4.0
    end

    return t
end


function get_soft_par()
    --softening parameter only calculated based on dwarf,
    --so if manual bodies is turned on the calculated s.p. may be too large
    if (UseOldSofteningLength == 1) then
        sp = calculateEps2(totalBodies, rscale_l[1], rscale_d[1], mass_l[1], mass_d[1])
    else
        sp = calculateEps2Dwarf(Dwarf.plummer{mass = mass_l[1], scaleLength = rscale_l[1]}, totalLightBodies)
    end
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
      sunVelx     = SunVelx,
      sunVely     = SunVely,
      sunVelz     = SunVelz,
      criterion   = criterion,
      useQuad     = true,
      useBestLike   = use_best_likelihood,
      BestLikeStart = eff_best_like_start,
      useVelDisp    = use_vel_disps,
      useBetaDisp   = use_beta_disps,
      useBetaComp   = use_beta_comp,
      useVlos       = use_vlos_comp,
      useDist       = use_avg_dist,
      usePropMot    = use_pm_comp,
      useMomentum   = use_momentum,
      Nstep_control = timestep_control,
      Ntsteps       = Ntime_steps,
      BetaSigma     = SigmaCutoff,
      VelSigma      = SigmaCutoff,
      DistSigma     = SigmaCutoff,
      PMSigma       = SigmaCutoff,
      MomentumSigma = SigmaCutoff,
      IterMax       = SigmaIter,
      BetaCorrect   = Correction,
      VelCorrect    = Correction,
      DistCorrect   = Correction,
      PMCorrect     = Correction,
      MomentumCorrect = Correction,
      SimpleOutput  = generateSimpleOutput,
      MultiOutput   = useMultiOutputs,
      OutputFreq    = freqOfOutputs,
      InitialOutput = generateInitialOutput,
      theta         = 1.0,
      LMC           = LMC_body,
      LMCfunction   = LMC_function,
      LMCmass       = LMC_Mass,
      LMCscale      = LMC_scaleRadius,
      LMCscale2     = LMC_cutoff,
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
        -- Setting bodies to 1 ensures worst case likelihood
        totalBodies = 1
        totalLightBodies = 1
    end

    if(run_null_potential == true and manual_bodies == true) then
        for i = 1, n do
            table.insert(finalPosition, lbrToCartesian(ctx, Vector.create(orbit_parameter_l[i], orbit_parameter_b[i], orbit_parameter_r[i])))
        end
        for i = 1, n do
            table.insert(finalVelocity, Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz))
        end
    elseif(run_null_potential == true) then
        print("placing dwarf at origin")
        finalPosition, finalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)
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
            local LMCposition = Vector.create(-1.1, -41.1, -27.9)
            local LMCvelocity = Vector.create(-57, -226, 221)
            local LMCmass = LMC_Mass
            local LMCscale = LMC_scaleRadius
            local LMCscale2 = LMC_scaleRadius/4 -- << PLACEHOLDER VALUE, PLEASE UPDATE
            local LMCfunction = 1   -- 1 = Plummer, 2 = hernquist 
            local LMCDynaFric = LMC_DynamicalFriction and 1 or 0
            local coulomb_log = CoulombLogarithm
            local ftime = evolveTime
            local tstop = revOrbTime
            local dt = ctx.timestep / 10.0      
            local masses    = dwarfMass  
            local rscales = rscale_t

            finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity = reverseOrbitS_LMC(potential, position, velocity, LMCposition, LMCvelocity, LMCmass, LMCfunction, LMCscale, LMCscale2, LMCDynaFric, coulomb_log, ftime, tstop, dt, masses, rscales)      
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

    if(ModelComponents == 2) then 
        for i = 1, n do
            local Model = predefinedModels.mixeddwarf{
                nbody       = totalBodies,
                nbody_baryon  = totalLightBodies,
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
            print(string.format("Dwarf %d bodies generation finished", i))
        end
        
    elseif(ModelComponents == 1) then
        for i = 1, n do
            local Model = predefinedModels.plummer{
                nbody       = totalBodies,
                prng        = prng,
                position    = finalPosition[i],
                velocity    = finalVelocity[i],
                mass        = mass_l[i],
                scaleRadius = rscale_l[i],
                ignore      = true
                }
            for _, row in ipairs(Model) do
                table.insert(firstModel, row)
                -- print(row)  -- < is this a debug function? double check...

            end
            print(string.format("Dwarf %d bodies generation finished", i))
        end

        -- firstModel = predefinedModels.plummer{
        --     nbody       = totalBodies,
        --     prng        = prng,
        --     position    = finalPosition,
        --     velocity    = finalVelocity,
        --     mass        = mass_l,
        --     scaleRadius = rscale_l,
        --     ignore      = false
        -- }
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
    elseif(ModelComponents > 0 and not manual_bodies) then        
        print("done")
        return firstModel
    else    
        print("Don't you want to simulate something?")
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
     betaBins  = bta_bins,

     -- Optional params
     L = {0.0, 0.0, 0.0}, --If any L components are nonzero, will use this L and LErr for momentum likelihood
     LErr = {0.0, 0.0, 0.0}, --This will overwrite any momentum values passed in through histogram. Input these as lua tables

     nRange = 0, --If non-zero, will use EMDRange values below. Overwrites values given in input histogram
     EMDRange = {} --Make sure this has an even number of elements and matches nRange. Input as a lua table
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