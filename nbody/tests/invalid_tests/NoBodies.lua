
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
   return { }
end

