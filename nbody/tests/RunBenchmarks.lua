
require "NBodyTesting"


local arg = {...}

local nbodyBinary = arg[1]

local nbodyFlags = getExtraNBodyFlags()
print("NBODY_FLAGS = ", nbodyFlags)



local sampleSeeds = { "670828913", "886885833", "715144259", "430281807", "543966758" }

local lineSep = string.rep("-", 80) .. "\n"

function runBenchmark(testName, seed, ...)
   local ret

   io.stdout:write(lineSep)
   print("Running benchmark " .. testName .. " with seed " .. seed)
   ret = os.readProcess(nbodyBinary,
                         "--ignore-checkpoint",
                         "--checkpoint-interval=-1", -- Disable checkpointing
                         "--debug-boinc", -- Prevent stderr from getting consumed with BOINC
                         "--timing",
                         "--input-file", "benchmark.lua",
                         "--output-cartesian",
                         "--output-file", testName .. ".out",
                         "--seed", seed,
                         nbodyFlags,
                         table.concat({...}, " ")
                      )
   io.stdout:write(ret)
   io.stdout:write(lineSep)
end

local nTimestep = 25


function runStandardPlummer(n, quad, crit, theta)
   local name
   name = string.format("Plummer_m=16_r=0.2__n=%d__quad=%s_%s=%.2f",
                        n, quad, crit, theta)
   runBenchmark(name, sampleSeeds[1], n, nTimestep, crit, theta, quad, 16, 0.2)
end

runStandardPlummer(100000, "true", "NewCriterion", 1.0)
runStandardPlummer(100000, "false", "NewCriterion", 1.0)

runStandardPlummer(100000, "true", "SW93", 1.0)
runStandardPlummer(100000, "false", "SW93", 1.0)

runStandardPlummer(100000, "true", "BH86", 0.6)
runStandardPlummer(100000, "false", "BH86", 0.6)

runStandardPlummer(100000, "true", "BH86", 0.5)
runStandardPlummer(100000, "false", "BH86", 0.5)


