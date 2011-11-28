
function makeHistogram()
   return HistogramParams.create()
end

function makePotential()
   return nil
end

function makeContext()
   return NBodyCtx.create{
      timestep   = calculateTimestep(16, 0.2),
      timeEvolve = 4.0,
      eps2       = calculateEps2(4096, 0.2),
      criterion  = "sw93",
      useQuad    = true,
      theta      = 1.0
   }
end

function makeBodies(ctx, potential)
   local bodies = predefinedModels.plummer{
      nbody       = 4096,
      prng        = DSFMT.create(argSeed),
      position    = Vector.create(0, 0, 0),
      velocity    = Vector.create(0, 0, 0),
      mass        = 16,
      scaleRadius = 0.2
   }

   return { { bodies } }
end

