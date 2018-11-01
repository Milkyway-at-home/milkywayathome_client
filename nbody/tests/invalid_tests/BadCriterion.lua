
prng = DSFMT.create(argSeed)
nbody = 4096


dwarfMass = 16
dwarfRadius = 0.2
reverseTime = 4.0
evolveTime = 3.945


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
      criterion  = "IDontExist",
      useQuad    = true,
      theta      = 1.0,
      BestLikeStart = 0.95,
      BetaSigma     = 2.5,
      VelSigma      = 2.5,
      BetaCorrect   = 1.111,
      VelCorrect    = 1.111,
      IterMax       = 6
   }
end

function makeBodies(ctx, potential)
   return predefinedModels.plummer{
      nbody       = nbody,
      prng        = prng,
      position    = Vector.create(0, 0, 0),
      velocity    = Vector.create(0, 0, 0),
      mass        = dwarfMass,
      scaleRadius = dwarfRadius,
      ignore      = false
   }
end

