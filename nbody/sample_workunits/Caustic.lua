arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")

prng = DSFMT.create(argSeed)

evolveTime       = arg[1]
reverseOrbitTime = arg[1] / arg[2]

r0  = arg[3]
light_r_ratio = arg[4]

dwarfMass  = arg[5]
light_mass_ratio = arg[6]

model1Bodies = 10000
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.32"

function makePotential()
   return  Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.caustic{ vhalo = 73, scaleLength = 12.0 }
   }
end

encMass = plummerTimestepIntegral(r0*light_r_ratio, sqr(r0) + sqr(r0/arg[4]) , dwarfMass, 1e-7)

function makeContext()
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (encMass + dwarfMass)),
      eps2       = calculateEps2(totalBodies, r0),
      criterion  = "TreeCode",
      useQuad    = true,
      theta      = 1.0
   }
end

-- Also required
function makeBodies(ctx, potential)
    local firstModel
    local finalPosition, finalVelocity = reverseOrbit{
        potential = potential,
       position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
       velocity  = Vector.create(-156, 79, 107),
        tstop     = reverseOrbitTime,
        dt        = ctx.timestep / 10.0
    }

    firstModel = predefinedModels.isotropic{
        nbody       = model1Bodies,
        prng        = prng,
        position    = finalPosition,
        velocity    = finalVelocity,
        mass1        = dwarfMass * arg[6],
	mass2       = dwarfMass - (dwarfMass * arg[6]),
        scaleRadius1 = r0,
	scaleRadius2 = r0/arg[4],
        ignore      = true
    }

    if (tonumber(arg[6]) == 1.0) then
      for i,v in ipairs(firstModel)
      do
     	v.ignore = false
      end
    

    else 
      count = 0
      while (count < arg[6] * totalBodies)
      do
    
      for i,v in ipairs(firstModel)
      do
         --Figure out properties
         r_2 = (finalPosition.x - v.position.x)^2 + (finalPosition.y - v.position.y)^2 + (finalPosition.z - v.position.z)^2
         r = r_2 ^ (0.5)
         scale1 = r0 
         lightDensityToMaxRatio = 1/((1 + r^2/scale1^2)^(5/2))

         --Get light sphere properties

         --Chance object is light_mass = light/dwarf
         chance = lightDensityToMaxRatio
         
         ---Do random calculation
         if(prng:random() > chance and v.ignore==true and count < arg[6] * totalBodies)
          then
              v.ignore=false
              count = count + 1
          end
       end
     end     
   end

return firstModel
end

function makeHistogram()
   return HistogramParams.create()
end


