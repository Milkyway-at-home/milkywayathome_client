--
-- Run a single Plummer sphere with no external potential for
-- benchmarking
--

args = {...}

nbody = args[1]
nTimestep = args[2]

criterion = args[3]
theta = args[4]
useQuad = (args[5] == "true")

mass = args[6]
radius = args[7]


assert(nbody, "Nbody not set")
assert(nTimestep, "nTimestep not set")

assert(criterion, "criterion not set")
assert(theta, "theta not set")
assert(args[5] == "true" or args[5] == "false", "useQuad must be \"true\" or \"false\"")

assert(mass, "mass not set")
assert(radius, "radius not set")



dt = calculateTimestep(mass, radius)


function makeHistogram()
   return HistogramParams.create()
end

function makePotential()
   return nil
end

function makeContext()
   return NBodyCtx.create{
      timestep   = dt,
      timeEvolve = nTimestep * dt,
      eps2       = calculateEps2(nbody, radius),
      criterion  = criterion,
      useQuad    = useQuad,
      theta      = theta
   }
end

function makeBodies(ctx, potential)
   return predefinedModels.plummer{
      nbody       = nbody,
      prng        = DSFMT.create(argSeed),
      position    = Vector.create(0, 0, 0),
      velocity    = Vector.create(0, 0, 0),
      mass        = mass,
      scaleRadius = radius
   }
end

