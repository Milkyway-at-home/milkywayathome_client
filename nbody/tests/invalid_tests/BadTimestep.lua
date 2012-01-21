

local nbody = 4096

local dwarfMass = 16
local dwarfRadius = 0.2
local reverseTime = 4.0
local evolveTime = 3.945


function makeHistogram()
   return HistogramParams.create()
end

function makePotential()
   return nil
end

function makeContext()
   return NBodyCtx.create{
      timestep   = -1.0,
      timeEvolve = evolveTime,
      eps2       = calculateEps2(nbody, dwarfRadius),
      criterion  = "sw93",
      useQuad    = true,
      theta      = 1.0
   }
end

function makeBodies(ctx, potential)
   return predefinedModels.plummer{
      nbody       = nbody,
      prng        = DSFMT.create(argSeed),
      position    = Vector.create(0, 0, 0),
      velocity    = Vector.create(0, 0, 0),
      mass        = dwarfMass,
      scaleRadius = dwarfRadius,
      ignore      = false
   }
end

