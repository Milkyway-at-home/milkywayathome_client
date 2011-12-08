
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

function makeContext()
   return NBodyCtx.create{
      timestep = 0.1,
      timeEvolve = 1.0,
      theta = 0.5,
      eps2 = 0.01
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

