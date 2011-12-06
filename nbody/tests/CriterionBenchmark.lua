
require "NBodyTesting"
require "persistence"

local bin = "milkyway_nbody"
local output = "/tmp/arst.out"

local nAvg = 1
local nTimestep = 40


local nbodies = { 1024, 10000, 20000, 32768, 50000, 75000, 100000 }
local thetas = { 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 }

local criteria = { "BH86", "SW93", "NewCriterion" }
local quads = { true, false }


function runSamples(deviceFlag)
   local timings = { }
   local mass = 16
   local radius = 0.2
   local args = {
      nbodyBin  = bin,
      input     = "benchmark.lua",
      output    = "/tmp/arst.out"
   }

   for _, nbody in ipairs(nbodies) do
      timings[nbody] = { }
      for _, crit in ipairs(criteria) do
         timings[nbody][crit] = { }
         for _, quad in ipairs(quads) do
            timings[nbody][crit][quad] = { }
            for _, theta in ipairs(thetas) do
               timings[nbody][crit][quad][theta] = { }
               local extraArgs = {
                  deviceFlag,
                  nbody,
                  nTimestep,
                  crit,
                  theta,
                  tostring(quad),
                  mass,
                  radius
               }

               local samples = { }
               args.extraArgs = extraArgs

               eprintf("Running test %s/%f/%d/%s...", crit, theta, nbody, tostring(quad))

               -- Take the best of n
               for i = 1, nAvg do
                  local output = runSimple(args)
                  samples[i] = findNumber(output, "run_time") or -1.0
               end

               local tMin = findMin(samples) / nTimestep
               timings[nbody][crit][quad][theta] = tMin
               eprintf("%f\n", tMin)
            end
         end
      end
   end

   -- exact ones
   for _, nbody in ipairs(nbodies) do
      timings[nbody]["Exact"] = { }
      timings[nbody]["Exact"][false] = { }
      timings[nbody]["Exact"][false][0.0] = { }

      local exactArgs = {
         deviceFlag,
         nbody,
         nTimestep,
         "Exact",
         0.0,
         "false",

         mass,
         radius
      }

      local samples = { }
      eprintf("Running test Exact/%d/false...", nbody)
      args.extraArgs = exactArgs
      for i = 1, nAvg do
         local output = runSimple(args)
         samples[i] = findNumber(output, "run_time") or -1.0
      end

      local tMin = findMin(samples) / nTimestep
      timings[nbody]["Exact"][false][0.0] = tMin
      eprintf("%f\n", tMin)
   end
   eprintf("\n")

   return timings
end

function boolToInt(b)
   if b then
      return 1
   else
      return 0
   end
end

function printCSV(timings)
   for _, nbody in ipairs(nbodies) do
      for _, crit in ipairs(criteria) do
         for _, quad in ipairs(quads) do
            for _, theta in ipairs(thetas) do
               printf("%d, \"%s\", %d, %f, %f\n",
                      nbody,
                      crit,
                      boolToInt(quad),
                      theta,
                      timings[nbody][crit][quad][theta])
            end
         end
      end
   end
end

local timings = runSamples("--disable-opencl")
--local timings = runSample("--device 0")
--local timings = runSample("--nthreads 1")
persistence.store("CriterionBenchmarkResults.lua", timings)
printCSV(timings)


