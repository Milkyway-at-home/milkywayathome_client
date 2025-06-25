-- /* Copyright (c) 2016-2018 Siddhartha Shelton */
-- /*               2019-2021    Eric Mendelsohn */
-- /*               2021              Tom Donlon */
-- /*      Rensselaer Polytechnic Institute      */

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- DEAR USER:
-- This is the GUI version of the lua parameter file.
-- It is meant to store the options for the GUI. 
-- The user should not need to make any changes to this file: 
-- rather, all changes should be made in the GUI.
-- As a result, things may not be easily human-readable. 
-- The GUI parser is very rudimentary, and changes to the format 
-- of this file are likely to break the GUI program.
-- 
-- If you want to make changes to this file, it is recommended that you
-- use the non-GUI version of nbody lite (settings.lua) instead.
-- 
-- The first half of the file is settings the Nbody simulation.
-- The second half of the file is code for generating the initial
-- conditions of the simulation, and is not something that
-- the typical end user needs to change.
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- FOR NOSY USERS:
-- The settings in this file have the postfix
-- "-- -- comment $ type | default value ^ top * bottom"
-- and add spaces before the comment to align it
--
-- If comment is "IGNORE", then the parameter will not be included in the GUI settings list
-- 
-- Limits by default set to [0,0]
--
-- Interface types:
-- entry: basic text box
-- q-entry: basic text box where the value is surrounded by quotes when outputed
-- l-entry or l-q-entry: wider version of the respective entry type
-- button: a true/false toggle button
-- short: displayed value can be incremented or decremented with buttons
--
-- Some (like bools) naturally do not need tops and bottoms and use filler
-- values.
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- STANDARD SETTINGS -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
nbodyMinVersion       = "1.92"   -- IGNORE

run_null_potential    = false    -- NULL POTENTIAL SWITCH $ button | 0 ^ 1 * 0
use_tree_code         = true     -- USE TREE CODE (NOT EXACT) $ button | 1 ^ 1 * 0
print_reverse_orbit   = false    -- PRINT REVERSE ORBIT SWITCH $ button | 0 ^ 1 * 0
print_out_parameters  = false    -- PRINT OUT ALL PARAMETERS $ button | 0 ^ 1 * 0
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --



-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- COORDINATE OPTIONS   -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

SunGCDist             = 8.0         -- Distance between Sun and Galactic Center (kpc) $ entry | 8.0 ^ 0 * 0
SunVelx               = 10.3        -- Sun's x-velocity (kpc/Gyr) (Hogg et al. (2005)) $ entry | 10.3 ^ 0 * 0
SunVely               = 229.2       -- Sun's y-velocity (kpc/Gyr) $ entry | 229.2 ^ 0 * 0
SunVelz               = 6.9         -- Sun's z-velocity (kpc/Gyr) $ entry | 6.9 ^ 0 * 0
LeftHandedCoords      = false       -- If true, work in left-handed galactocentric cartesian coordinates $ button | 0 ^ 0 * 0
                                    -- (e.g. the Sun is located at positive X) 

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --



-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- LMC  SETTINGS  -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
LMC_body              = true        -- RUN WITH LMC $ button | 1 ^ 1 * 0
LMC_scaleRadius       = 15          -- NO COMMENT $ entry | 15 ^ 0 * 0
LMC_Mass              = 449865.888  -- NO COMMENT $ l-entry | 449865.888 ^ 0 * 0
LMC_DynamicalFriction = true        -- LMC DYNAMICAL FRICTION SWITCH (IGNORED IF NO LMC) $ button | 1 ^ 1 * 0
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --



-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- DWARF PARAMETERS  -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- Switch for whether to use a single component or -- -- -- -- --
-- -- -- -- -- two component (baryon & dark matter) model for  -- -- -- -- --
-- -- -- -- -- the generated dwarf galaxy -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- ModelComponent Options:
-- --                2 - TWO COMPONENT MODEL    -- -- -- -- -- -- -- -- -- --
-- --                1 - SINGLE COMPONENT MODEL -- -- -- -- -- -- -- -- -- --
-- --                0 - NO DWARF MODEL         -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
ModelComponents   = 2         -- TWO COMPONENTS SWITCH $ short | 2 ^ 2 * 0
manual_bodies     = false     -- USE THE MANUAL BODY LIST $ button | 0 ^ 1 * 0
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- Change the parameters below to specify the settings   -- -- --
-- -- -- -- -- for the generation of the dwarf bodies -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
totalBodies      = 40000        -- Number of Bodies $ entry | 40000 ^ 0 * 0
evolveTime       = 3.0          -- Forward Time (Gyr) $ entry | 3.0 ^ 0 * 0
revOrbTime       = 3.0          -- Reverse Orbit Time (Gyr) $ entry | 3.0 ^ 0 * 0
rscale_l         = 0.3          -- Baryonic Radius (kpc) $ entry | 0.3 ^ 0 * 0
light_r_ratio    = 0.2          -- Baryonic Radius / (Baryonic Radius + Dark Matter Radius) $ entry | 0.2 ^ 0 * 0
mass_l           = 45.0         -- Baryonic Mass (Structure Mass Units) $ entry | 45.0 ^ 0 * 0
light_mass_ratio = 0.1          -- Baryonic Mass / (Baryonic Mass + Dark Matter Mass) $ entry | 0.1 ^ 0 * 0
orbit_parameter_l  = 258        -- Galactic coordinates of dwarf position (deg) $ entry | 258 ^ 0 * 0
orbit_parameter_b  = 45.8       -- NO COMMENT $ entry | 45.8 ^ 0 * 0
orbit_parameter_r  = 21.5       -- Distance from Sun to dwarf (kpc) $ entry | 21.5 ^ 0 * 0
orbit_parameter_vx = -185.5     -- Galactocentric (no Solar motion) velocities of dwarf (km/s) $ entry | -185.5 ^ 0 * 0
orbit_parameter_vy = 54.7       -- NO COMMENT $ entry | 54.7 ^ 0 * 0
orbit_parameter_vz = 147.4      -- NO COMMENT $ entry | 147.4 ^ 0 * 0
manual_body_file = "manual_bodies_example.in" -- (Optional) Manual bodies list. Can be nil. $ l-q-entry | manual_bodies_example.in ^ 0 * 0
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --



-- -- -- -- -- -- -- --  OUTPUT SETTINGS  -- -- -- -- -- -- -- -- -- -- -- -- 
generateInitialOutput = false     -- Outputs the initial bodies file right after dwarf generation $ button | 0 ^ 1 * 0

generateSimpleOutput = true       -- Simple output file includes: x, y, z, vx, vy, vz, mass $ button | 1 ^ 1 * 0
-- Full output file includes: x, y, z, l, b, r, vx, vy, vz, mass, vlos, pmra, pmdec, [lambda, beta]
-- NOTE: Lambda and Beta are optional and will only be included if the histogram parameters are set in makeHistogram()
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 



-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- TIME CONTROL OPTIONS -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- Control how often to output data and   -- -- -- -- -- -- -- --
-- -- -- -- -- whether to manually set timestep size  -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
useMultiOutputs      = false  -- WRITE MULTIPLE OUTPUTS $ button | 0 ^ 1 * 0
freqOfOutputs        = 100    -- FREQUENCY OF WRITING OUTPUTS $ entry | 100 ^ 0 * 0

timestep_control     = false  -- control number of steps $ button | 0 ^ 1 * 0
Ntime_steps          = 3000   -- number of timesteps to run (ignored if timestep_control == false) $ entry | 3000 ^ 0 * 0
                              -- (<1 Myr timestep size strongly recommended)
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --



-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- --  MISC OPTIONS  -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- Any other options for the simulation   -- -- -- -- -- -- -- --
-- -- -- -- -- are provided here                      -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
use_max_soft_par      = false       -- limit the softening parameter value to a max value $ button | 0 ^ 1 * 0
                                    -- (NOTE: This is turned on automatically if manual_bodies is true,
                                    -- since the softening parameter is determined only by the dwarf bodies)
max_soft_par          = 0.8         -- kpc, if switch above is turned on, use this as the max softening parameter $ entry | 0.8 ^ 0 * 0
UseOldSofteningLength = false       -- If true, uses old softening length formula from v1.76 and earlier $ button | 0 ^ 1 * 0
                                    -- (this is only useful to compare with/match simulations 
                                    --  that were run before v1.80)
CoulombLogarithm      = 0.470003629 -- (ln(1.6)) COULOMB LOGARITHM USED IN DYNAMICAL FRICTION $ entry | 0.470003629 ^ 0 * 0
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- END GUI

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- -- POTENTIAL OPTIONS -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- Control the shape of the external gravitational -- -- -- -- --
-- -- -- -- -- potential of the simulation.  -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- Add/remove/change potentials in the -- -- -- -- -- -- -- -- --
-- -- -- -- -- Potential.create{ } brackets to edit the potential.   -- -- --
-- -- -- -- -- More information about avaiable potentials   -- -- -- -- -- --
-- -- -- -- -- coming soon.   -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
function makePotential()
   if(run_null_potential == true) then
       print("running in null potential")
       return nil
   else
        --NOTE: To exclude a component from the potential, set component to "<component_name>.none" and include only an arbitrary "mass" argument
        return  Potential.create{
            spherical = Spherical.hernquist{ mass  = 1.52954402e5, scale = 0.7 },
            disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
            disk2     = Disk.none{ mass = 3.0e5 },
            halo      = Halo.logarithmic{ vhalo = 74.61, scaleLength = 12.0, flattenZ = 1.0 }
        }--vhalo = 74.61 kpc/gy = 73 km/s
   end
end
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 



-- XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX --
-- XX -- -- -- -- -- -- -- -- DEVELOPER OPTIONS -- -- -- -- -- -- -- -- XX --
-- XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX --
-- XX -- -- -- The actual setup for the simulation is done below. -- -- XX --
-- XX -- -- -- The typical user should not have to change anything   -- XX --
-- XX -- -- -- in this file beyond this point.  -- -- -- -- -- -- -- -- XX --
-- XX -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- XX --
-- XX -- -- -- -- -- -- -- -- -- WARNING  -- -- -- -- -- -- -- -- -- -- XX --
-- XX -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- XX --
-- XX -- -- -- Changing anything below this point may result in   -- -- XX --
-- XX -- -- -- undesired behavior of the Nbody Lite software and  -- -- XX --
-- XX -- -- -- may generate inaccurate/unphysical results   -- -- -- -- XX --
-- XX -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- XX --
-- XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX XX --

-- -- -- -- -- -- -- -- -- CHECK TIMESTEPS -- -- -- -- -- -- -- -- 
TooManyTimesteps = 0

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
        TooManyTimesteps = 1
        t = evolveTime/4.0
    end

    return t
end


function get_soft_par()
    --softening parameter only calculated based on dwarf,
    --so if manual bodies is turned on the calculated s.p. may be too large
    sp = calculateEps2(totalBodies, rscale_l, rscale_d, mass_l, mass_d, UseOldSofteningLength)

    if ((manual_bodies or use_max_soft_par) and (sp > max_soft_par^2)) then --dealing with s.p. squared
        print("Using maximum softening parameter value of " .. tostring(max_soft_par) .. " kpc")
        return max_soft_par^2
    else
        return sp
    end
end

-- A lot of this is hard-coded to the default values.
-- This is because lite users don't need to optimize histograms, 
-- and having all the options for histogram optimization in the 
-- settings file just makes things bloated.
function makeContext()
   return NBodyCtx.create{
      timeEvolve    = evolveTime,
      timeBack      = revOrbTime,
      timestep      = get_timestep(),
      eps2          = get_soft_par(),
      b             = orbit_parameter_b,
      r             = orbit_parameter_r,
      vx            = orbit_parameter_vx,
      vy            = orbit_parameter_vy,
      vz            = orbit_parameter_vz,
      criterion     = criterion,
      useQuad       = true,
      useBestLike   = false,
      BestLikeStart = 0.98,
      useVelDisp    = false,
      useBetaDisp   = false,
      useBetaComp   = false,
      useVlos       = false,
      useDist       = false,
      Nstep_control = timestep_control,
      Ntsteps       = Ntime_steps,
      BetaSigma     = 2.5,
      VelSigma      = 2.5,
      DistSigma     = 2.5,
      PMSigma       = 2.5,
      IterMax       = 6,
      BetaCorrect   = 1.111,
      VelCorrect    = 1.111,
      DistCorrect   = 1.111,
      PMCorrect     = 1.111,
      leftHanded    = LeftHandedCoords,
      useContBins   = false,
      bleedInRange  = 1,
      SimpleOutput  = generateSimpleOutput,
      MultiOutput   = useMultiOutputs,
      OutputFreq    = freqOfOutputs,
      InitialOutput = generateInitialOutput,
      theta         = 1.0,
      LMC           = LMC_body,
      LMCmass       = LMC_Mass,
      LMCscale      = LMC_scaleRadius,
      LMCDynaFric   = LMC_DynamicalFriction,
      coulomb_log   = CoulombLogarithm,
      calibrationRuns = 0
   }
end



function makeBodies(ctx, potential)
  local firstModel
  local finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity
    if TooManyTimesteps == 1 then
        totalBodies = 1
    end

    if(run_null_potential == true and manual_bodies == true) then
        finalPosition = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r))
        finalVelocity = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz)
    elseif(run_null_potential == true) then
        print("placing dwarf at origin")
        finalPosition, finalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)
    else 
    	if (LMC_body) then
    		finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity = reverseOrbit_LMC{
	            potential   = potential,
	            position    = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
	            velocity    = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
	            LMCposition = Vector.create(-1.1, -41.1, -27.9),
	            LMCvelocity = Vector.create(-57, -226, 221), 
                    LMCmass     = LMC_Mass,
                    LMCscale    = LMC_scaleRadius,
                    LMCDynaFric = LMC_DynamicalFriction,
                    coulomb_log = CoulombLogarithm,
                    ftime       = evolveTime,
	            tstop       = revOrbTime,
	            dt          = ctx.timestep / 10.0,
	            sunGCDist   = SunGCDist
	            }

              
	    else
	        finalPosition, finalVelocity = reverseOrbit{
	            potential = potential,
	            position  = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
	            velocity  = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
	            tstop     = revOrbTime,
	            dt        = ctx.timestep / 10.0,
	            sunGCDist = SunGCDist
	            }
         end
    end
    
     if(print_reverse_orbit == true) then
        print('Printing reverse orbit')
        if (LMC_body) then
            local placeholderPos, placeholderVel, LMCplaceholderPos, LMCplaceholderVel = PrintReverseOrbit_LMC{
                potential = potential,
                position  = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
                velocity  = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
	        LMCposition = Vector.create(-1.1, -41.1, -27.9),
	        LMCvelocity = Vector.create(-57, -226, 221), 
                LMCmass     = LMC_Mass,
                LMCscale    = LMC_scaleRadius,
                LMCDynaFric = LMC_DynamicalFriction,
                tstop     = .14,
                tstopf    = .20,
                dt        = ctx.timestep / 10.0,
                coulomb_log = CoulombLogarithm
            }
        else
            local placeholderPos, placeholderVel = PrintReverseOrbit{
                potential = potential,
                position  = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
                velocity  = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
                tstop     = .14,
                tstopf    = .20,
                dt        = ctx.timestep / 10.0
            }
        end
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
    elseif(ModelComponents > 0 and not manual_bodies) then
        return firstModel
    else    
        print("Don't you want to simulate something?")
    end
    
end

assert(argSeed ~= nil, "Expected seed") -- STILL EXPECTING SEED AS INPUT FOR THE FUTURE
argSeed = 7854614814 -- -- SETTING SEED TO FIXED VALUE
prng = DSFMT.create(argSeed)


if(ModelComponents == 1) then
   dwarfMass = mass_l
   rscale_t  = rscale_l
   rscale_d  = 1.0
   mass_d    = 0.0
else
   dwarfMass = mass_l / light_mass_ratio
   rscale_t  = rscale_l / light_r_ratio
   rscale_d  = rscale_t *  (1.0 - light_r_ratio)
   mass_d    = dwarfMass * (1.0 - light_mass_ratio)
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
