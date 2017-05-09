arg = { ... }
-- /* Copyright (c) 2016 Siddhartha Shelton */
assert(#arg == 4, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")

prng = DSFMT.create(argSeed)

evolveTime       = tonumber(arg[1])
reverseOrbitTime = evolveTime / tonumber(arg[2])

print(evolveTime)

r0  = tonumber(arg[3])

dwarfMass  = tonumber(arg[4])

model1Bodies = 20000
totalBodies = model1Bodies

function makePotential()
    
return  nil
end



encMass = plummerTimestepIntegral(r0, dwarfMass, 1e-7)

-- This is also required
function makeContext()
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (dwarfMass)),
      eps2       = calculateEps2(totalBodies, r0),
      criterion  =  "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end

-- Also required
function makeBodies(ctx, potential)
   local firstModel
    local finalPosition, finalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)

   firstModel = predefinedModels.plummer{
      nbody       = model1Bodies,
      prng        = prng,
      position    = finalPosition,
      velocity    = finalVelocity,
      mass        = dwarfMass,
      scaleRadius = r0,
      ignore      = false
   }

   return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = -100,
     lambdaEnd = 100,
     lambdaBins = 50,
     betaStart = -100,
     betaEnd = 100,
     betaBins = 1
}
end