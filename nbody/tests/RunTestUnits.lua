
require "NBodyTesting"
require "persistence"

local arg = {...}

assert(#arg == 5, "Test driver expected 5 arguments got " .. #arg)

local nbodyBinary = arg[1]
local testDir = arg[2]
local testName = arg[3]
local histogramName = arg[4]
local testBodies = arg[5]

local nbodyFlags = getExtraNBodyFlags()
eprintf("NBODY_FLAGS = %s\n", nbodyFlags)

math.randomseed(os.time())

-- Pick one of the random seeds used in generating these tests
local testSeeds = {"670828913"} --no longer testing multiple seeds
local testSeed = testSeeds[math.random(1, #testSeeds)]
--local testSeed = testSeeds[5]


refResults = {
   ["model_1"] = {
      ["10000"] = {
         ["670828913"] = 771.495957702001306,
      }
   },

   ["model_2"] = {
      ["10000"] = {
         ["670828913"] = 1599.2064279689132,
      }
   },

   ["model_3"] = {
      ["10000"] = {
         ["670828913"] = 2668.371681712826103,
      }
   },

   ["model_4"] = {
      ["10000"] = {
         ["670828913"] = 471.0690774528726,
      }
   },

   ["model_5"] = {
      ["10000"] = {
         ["670828913"] = 442.85633731928453,
      }
   },

   ["model_5_bounds_test"] = {
      ["10000"] = {
         ["670828913"] = 3483.051322529342,
      }
   },

   ["model_6"] = {
      ["10000"] = {
         ["670828913"] = 427.723338242557702,
      }
   },

   ["model_7"] = {
      ["10000"] = {
         ["670828913"] = 343.8076885597965,
      }
   },

   ["model_8"] = {
      ["10000"] = {
         ["670828913"] = 329.129417006217807,
      }
   },

   ["model_9"] = {
      ["10000"] = {
         ["670828913"] = 129.10766206161105,
      }
   },

   ["model_ninkovic"] = {
      ["10000"] = {
         ["670828913"] = 431.936920465970331,
      }
   },

   ["model_triaxial"] = {
      ["10000"] = {
         ["670828913"] = 698.3865985898103,
      }
   },

   ["model_newhist1"] = {
      ["10000"] = {
         ["670828913"] = 2726.180023873662,
      }
   },

   ["model_newhist2"] = {
      ["10000"] = {
         ["670828913"] = 1468.471832532997723,
      }
   },

   ["model_newhist3"] = {
      ["10000"] = {
         ["670828913"] = 2515.5228997740865,
      }
   },

   ["model_LMC"] = {
      ["10000"] = {
         ["670828913"] = 1714.7154077794053,
      }
   },

   ["model_bar"] = {
      ["10000"] = {
         ["670828913"] = 548.139288339098698,
      }
   },

   ["model_LMC_bar"] = {
      ["10000"] = {
         ["670828913"] = 2211.599321169714585,
      }
   }
}


function resultCloseEnough(a, b)
   return math.abs(a - b) < 1.0e-10
end

errFmtStr = [[
Result differs from expected:
   Expected = %20.15f  Actual = %20.15f  |Difference| = %20.15f
]]

function runCheckTest(testName, histogram, seed, nbody, ...)
   local fileResults, bodyResults
   local ret, result

   if not generatingResults then
      -- Check if the result exists first so we don't waste time on a useless test
      fileResults = assert(refResults[testName], "Didn't find result for test file")
      bodyResults = assert(fileResults[nbody], "Didn't find result with matching bodies")
      refResult = assert(bodyResults[seed], "Didn't find result with matching seed")
   end

   --eprintf("CHECKTEST - Before runFullTest\n")

   ret = runFullTest{
      nbodyBin  = nbodyBinary,
      testDir   = testDir,
      testName  = testName,
      histogram = histogram,
      seed      = seed,
      cached    = false,
      extraArgs = { nbody }
   }

   --eprintf(ret.."\n")
   --eprintf("CHECKTEST - Before findLikelihood\n")

   result = findLikelihood(ret, false)

   --eprintf("CHECKTEST - Before write(ret)\n")

   io.stdout:write(ret)

   if generatingResults then
      io.stderr:write(string.format("Test result: %d, %d, %s: %20.15f\n", nbody, seed, testName, result))
      return false
   end

   if result == nil then
      return true
   end

   --eprintf("CHECKTEST - Before notClose\n")

   local notClose = not resultCloseEnough(refResult, result)
   if notClose then
      io.stderr:write(string.format(errFmtStr, refResult, result, math.abs(result - refResult)))
   end

   return notClose
end

-- return true if passed
function testProbabilistic(resultFile, testName, histogram, nbody, iterations)
   local testTable, histTable, answer
   local resultTable = persisence.load(resultFile)
   assert(resultTable, "Failed to open result file " .. resultFile)

   testTable = assert(resultTable[testName], "Did not find result for test " .. testName)
   histTable = assert(testTable[nbody], "Did not find result for nbody " .. tostring(nbody))
   answer = assert(histTable[nbody], "Did not find result for histogram " .. histogram)

   local minAccepted = answer.mean - 3.0 * answer.stddev
   local maxAccepted = answer.mean + 3.0 * answer.stddev

   local result = 0.0
   local z = (result - answer.mean) / answer.stddev


   return true
end



function getResultName(testName)
   return string.format("%s__results.lua", testName)
end

if runCheckTest(testName, histogramName, testSeed, testBodies) then
   os.exit(1)
end
