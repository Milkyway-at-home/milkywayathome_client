arg = { ... } -- -- TAKING USER INPUT
argSeed = 34086709 -- -- SETTING SEED TO FIXED VALUE
prng = DSFMT.create(argSeed) 

evolveTime       = 0.1    -- Forward Time
time_ratio       = 1.0    -- Forward Time / Backward Time
rscale_l         = 0.2    -- Baryonic Radius
light_r_ratio    = 0.2   -- Baryonic Radius / (Baryonic Radius + Dark Matter Radius)
mass_l           = 12.0    -- Baryonic Mass (Structure Mass Units)
light_mass_ratio = 0.2
dwarfMass = mass_l / light_mass_ratio
rscale_t  = rscale_l / light_r_ratio
rscale_d  = rscale_t *  (1.0 - light_r_ratio)
mass_d    = dwarfMass * (1.0 - light_mass_ratio)

totalBodies           = arg[1]   -- -- NUMBER OF BODIES           -- --
nbodyLikelihoodMethod = "EMD"   -- -- HIST COMPARE METHOD        -- --
nbodyMinVersion       = "1.76"  -- -- MINIMUM APP VERSION        -- --
run_null_potential    = false   -- -- NULL POTENTIAL SWITCH      -- --
use_tree_code         = true    -- -- USE TREE CODE NOT EXACT    -- --
LMC_body              = false    -- -- PRESENCE OF LMC            -- --
lda_bins        = 50      -- number of bins in lamdba direction
lda_lower_range = -150    -- lower range for lambda
lda_upper_range = 150     -- upepr range for lamdba
bta_bins        = 1       -- number of beta bins. normally use 1 for 1D hist
bta_lower_range = -15     -- lower range for beta
bta_upper_range = 15      -- upper range for beta
SigmaCutoff          = 2.5     -- -- sigma cutoff for outlier rejection DO NOT CHANGE -- --
SigmaIter            = 6       -- -- number of times to apply outlier rejection DO NOT CHANGE -- --
Correction           = 1.111   -- -- correction for outlier rejection   DO NOT CHANGE -- --
use_best_likelihood  = true    -- use the best likelihood return code
best_like_start      = 0.98    -- what percent of sim to start
use_beta_disps       = true    -- use beta dispersions in likelihood
use_vel_disps        = false    -- use velocity dispersions in likelihood
use_beta_comp        = false  -- calculate average beta, use in likelihood
use_vlos_comp        = false  -- calculate average los velocity, use in likelihood
use_avg_dist         = false  -- calculate average distance, use in likelihood
useMultiOutputs       = false        -- -- WRITE MULTIPLE OUTPUTS       -- --
freqOfOutputs         = 6            -- -- FREQUENCY OF WRITING OUTPUTS -- --
timestep_control     = false         -- -- control number of steps      -- --
Ntime_steps          = 3000            -- -- number of timesteps to run   -- --
        
orbit_parameter_l  = 218
orbit_parameter_b  = 53.5
orbit_parameter_r  = 28.6
orbit_parameter_vx = -156 
orbit_parameter_vy = 79 
orbit_parameter_vz = 107

dwarfMass = 16
dwarfRadius = 0.2
        
function makePotential()
   if(run_null_potential == true) then
       print("running in null potential")
       return nil
   else
        --NOTE: To exculde a component from the potential, set component to "<component_name>.none" and include only an arbitrary "mass" argument
   	return Potential.create{
      spherical = Spherical.hernquist{ mass = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      disk2     = Disk.none{ mass = 3.0e5 },
      halo      = Halo.triaxial{ vhalo = 116,
                                 scaleLength = 16.3,
                                 flattenZ = 1.43,
                                 flattenX = 1.26,
                                 flattenY = 1.0,
                                 triaxAngle = 96
                              }
   }
   end
end

function makeContext()
   soften_length  = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
   return NBodyCtx.create{
      timeEvolve  = 3.945,
      timeBack    = revOrbTime,
      timestep    = calculateTimestep(dwarfMass, dwarfRadius),
      eps2        = calculateEps2(totalBodies, soften_length),
      b           = orbit_parameter_b,
      r           = orbit_parameter_r,
      vx          = orbit_parameter_vx,
      vy          = orbit_parameter_vy,
      vz          = orbit_parameter_vz, 
      criterion   = "sw93",
      useQuad     = true,
      useBestLike   = use_best_likelihood,
      BestLikeStart = .95,
      useVelDisp    = use_vel_disps,
      useBetaDisp   = use_beta_disps,
      useBetaComp   = use_beta_comp,
      useVlos       = use_vlos_comp,
      useDist       = use_avg_dist,
      Nstep_control = timestep_control,
      Ntsteps       = Ntime_steps,
      BetaSigma     = 2.5,
      VelSigma      = 2.5,
      DistSigma     = 2.5,
      IterMax       = 6,
      BetaCorrect   = 1.111,
      VelCorrect    = 1.111,
      DistCorrect   = 1.111,
      MultiOutput   = useMultiOutputs,
      OutputFreq    = freqOfOutputs,
      theta         = 1.0,
      LMC           = LMC_body
   }
end

function makeBodies(ctx, potential)
   local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 29.5)),
      velocity  = Vector.create(-183, 101, 107),
      tstop     = 4.0,
      dt        = ctx.timestep / 10.0
   }

   return predefinedModels.plummer{
      nbody       = totalBodies,
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
