-- A copy of orphan_test with a plummer sphere implementation in Lua

if isBOINCApplication and argv == nil then
  error("Running as BOINC application, but no server arguments")
end

function myPlummer(prng, nbody, mass, ignore, rShift, vShift, radiusScale)
   local massPerParticle = mass / nbody
   local velScale = sqrt(mass / radiusScale)
   local function randomR()
      return 1.0 / sqrt(pow(prng:genrandCloseOpen(), -2.0 / 3.0) - 1.0)
   end

  local function selectG()
     local x, y
     repeat
        x = prng:random(0.0, 1.0)
        y = prng:random(0.0, 0.1)
     until y <= sqr(x) * pow(1.0 - sqr(x), 3.5)
     return x
  end

  local function randomVel(r)
     local x = selectG()
     return sqrt2 * x / pow(1.0 + sqr(r), 0.25);
  end

  local function pickShell(r)
     local vec, rsq, rsc
     repeat
        vec = prng:randomVector()
        rsq = vec * vec;
     until rsq <= 1.0

     rsc = r / sqrt(rsq)
     return rsc * vec
  end

  local bodies = { }
  for i = 1, nbody do
      local r = randomR()
      bodies[i] = Body.create{ mass = massPerParticle,
                               ignore = ignore,
                               position = pickShell(radiusScale * r) + rShift,
                               velocity = pickShell(velScale * randomVel(r)) + vShift,
                            }
   end

  return bodies
end

if argv ~= nil then
   table.foreach(argv, print)

   -- Use the arguments from the server. Tables are indexed from 1
   dwarfMass = argv[1]
   dwarfRadius= argv[2]
   reverseTime = argv[3]
   evolveTime = argv[4]
else
  dwarfMass = 16
  dwarfRadius = 0.2
  reverseTime = 4.0
  evolveTime = 3.945
end

-- Make sure to use the seed from the server
prng = DSFMT.create(argSeed)
nbody = 100


-- This is a required function. You get the defaults with no arguments
function makeHistogram()
   return HistogramParams.create()
end

-- This is a required function
function makePotential()
   local disk = Disk.miyamotoNagai{ mass = 4.45865888e5,
                                    scaleLength = 6.5,
                                    scaleHeight = 0.26
                                 }

   local halo = Halo.logarithmic{ vhalo = 73,
                                  scaleLength = 12.0,
                                  flattenZ = 1.0
                               }

   local spherical = Spherical.spherical{ mass = 1.52954402e5,
                                          scale = 0.7
                                       }

  local pot = Potential.create{ disk = disk,
                                halo = halo,
                                spherical = spherical
                             }
  return pot
end

-- This is also required
function makeContext()
    local r0 = dwarfRadius
   local eps = r0 / (10.0 * sqrt(nbody))
   local ctx = NBodyCtx.create
   {
   criterion = "sw93",
   useQuad = true,
   theta = 1.0,
   timeEvolve = evolveTime,
    eps2 = eps * eps,
   timestep = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / dwarfMass)
   }

   return ctx
end

-- Also required
function makeBodies(ctx, potential)
   local iniVel = Vector.create(-156, 79, 107)
   local iniPos = Vector.create(218, 53.5, 28.6)

   local ic = InitialConditions.create{ context = ctx,
                                                          velocity = iniVel,
                                                          position = iniPos }

  local orbitTimestep = ctx.timestep / 10.0
  local finalPosition = reverseOrbit{ potential = potential,
                                                        initialConditions = ic,
                                                        tstop = reverseTime,
                                                        dt = orbitTimestep
                                                      }
  local bodies = myPlummer(prng,
                           nbody,
                           dwarfMass,
                           false,
                           finalPosition.position,
                           finalPosition.velocity,
                           dwarfRadius)
   return bodies
end

