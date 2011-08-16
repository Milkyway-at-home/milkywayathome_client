
prng = DSFMT.create(argSeed)
nbody = 1000

dwarfMass = 16
dwarfRadius = 0.2
reverseTime = 4.0
evolveTime = 0.1


function makeHistogram()
   return HistogramParams.create()
end

function makePotential()
   return nil
end

function makeContext()
   return NBodyCtx.create{
      timestep   = calculateTimestep(dwarfMass, dwarfRadius),
      timeEvolve = evolveTime,
      eps2       = calculateEps2(nbody, dwarfRadius),
      criterion  = "BH86",
      useQuad    = false,
      theta      = 0.5
   }
end

function makeBodies(ctx, potential)
   return predefinedModels.plummer{
      nbody       = nbody,
      prng        = prng,
      position    = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
      velocity    = Vector.create(20, 0, 0),
      mass        = dwarfMass,
      scaleRadius = dwarfRadius,
      ignore      = false
   }
end

