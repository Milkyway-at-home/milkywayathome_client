
args = {...}

assert(#args == 5, "5 arguments required")

assert(args[1] == "hello",  "string argument failed")
assert(args[2] == "-4",     "Negative single digit argument failed")
assert(args[3] == "-3.14",  "Negative single digit argument failed")
assert(args[4] == "1.23",  "Positive single digit argument failed")
assert(args[5] == "13434", "Multi digit integer argument failed")

-- Just whatever so that the file pasess
function makePotential()
   return nil
end

function makeHistogram()
   return HistogramParams.create()
end

-- rdh edited
dec = 9.0   -- -- number of decimals to round to (default: 9.0)
function round(num, places)
    local mult = 10.0^(places)
    return floor(num * mult + 0.5) / mult
  end

function makeContext()
   return NBodyCtx.create{
      -- rdh edited
      dwarfn = 1,
      b           = {round( -44.328, dec)},
      r           = {round( 62.4,    dec)},
      vx          = {round( 21.99,   dec)},
      vy          = {round( -201.36, dec)},
      vz          = {round( 171.25,  dec)},
      timestep      = 0.1,
      timeEvolve    = 1.0,
      theta         = 0.5,
      eps2          = 0.01,
      BestLikeStart = 0.98,
      BetaSigma     = 2.5,
      VelSigma      = 2.5,
      DistSigma     = 2.5,
      PMSigma       = 2.5,
      MomentumSigma = 2.5,
      BetaCorrect   = 1.111,
      VelCorrect    = 1.111,
      DistCorrect   = 1.111,
      PMCorrect     = 1.111,
      MomentumCorrect = 1.111,
      IterMax       = 6
   }
end

function makeBodies2()
   return { { }, { } }
end

function makeBodies()
   return { Body.create{
               mass = 0.1,
               position = Vector.create(0, 0, 0),
               velocity = Vector.create(0, 0, 0)
            }
         }
end

