
require "NBodyTesting"
require "persistence"

local arg = {...}
local nbodyBinary = assert(arg[1])

local nbodyFlags = getExtraNBodyFlags()
eprintf("NBODY_FLAGS = %s\n", nbodyFlags)

local testDir = "orphan_models"

math.randomseed(os.time())

-- When generating results, write results to stdout and everything
-- else to stderr to make things easier

function generateSampleUnits(outputName, testName, histograms, nbodies, iterations)
   local file, results, newResults = false
   local allResults = persistence.load(outputName)
   if allResults == nil then
      eprintf("Could not open file '%s', creating new\n", outputName)
      allResults = { }
      newResults = true
   end

   local results
   if allResults[testName] == nil then
      allResults[testName] = { }
   end
   results = allResults[testName]

   for _, nbody in ipairs(nbodies) do
      if results[nbody] == nil then
         results[nbody] = { }
      end
      local histTable = results[nbody]
      local nSamples = { } -- we need to know how many are there originally to avoid clobbering old ones

      for _, histogram in ipairs(histograms) do
         if histTable[histogram] == nil then
            histTable[histogram] = { }
            histTable[histogram]["samples"] = { }
         end
         nSamples[histogram] = #histTable[histogram]["samples"]
      end


      for i = 1, iterations do
         local seed = randomSeed()

         local cached = false -- we only need to run the sim once and match against all histograms
         for _, histogram in ipairs(histograms) do
            local this = histTable[histogram]["samples"]
            eprintf("Running test '%s' with %d bodies. cached = %s\n", testName, nbody, tostring(cached))
            local output = runFullTest{
               nbodyBin  = nbodyBinary,
               testDir   = testDir,
               testName  = testName,
               histogram = histogram,
               seed      = seed,
               cached    = cached,
               extraArgs = { nbody }
            }

            local result = findLikelihood(output, false)
            assert(result, "Did not find result in test output")
            eprintf("  iteration %4d (seed %6d) = %f\n", i, seed, result)
            this[i + nSamples[histogram]] = result
            cached = true
         end
      end

      for _, histogram in pairs(histograms) do
         local histItem = histTable[histogram]
         local histSamples = histItem["samples"]
         histItem["mean"], histItem["stddev"] = calcStats(histSamples)
         histItem["min"], histItem["max"] =  findMinMax(histSamples)
      end
   end

   if not newResults then
      os.rename(outputName, outputName .. ".bak")
   end

   local tmpName = outputName .. ".tmp"
   persistence.store(tmpName, allResults)
   os.rename(tmpName, outputName)
end

local testName = assert(arg[2])
local outputFile = arg[3] or "GeneratedResults.lua"
local nbody = arg[4] or 10000
local nIter = arg[5] or 10

local histograms = {"orphan_model_histogram", "model_1_self.hist"}
local nbodies = { nbody }



for i = 1, nIter do
   generateSampleUnits(outputFile, testName, histograms, nbodies, 1)
end




