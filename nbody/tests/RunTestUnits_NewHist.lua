
require "NBodyTesting"
require "persistence"

local arg = {...}

assert(#arg == 11, "Test driver expected 11 arguments got " .. #arg)

local nbodyBinary = arg[1]
local testDir = arg[2]
local testName = arg[3]
local histogramName = arg[4]
local testBodies = arg[5]
local b = arg[6]
local r = arg[7]
local vx = arg[8]
local vy = arg[9]
local vz = arg[10]

local nbodyFlags = getExtraNBodyFlags()
eprintf("NBODY_FLAGS = %s\n", nbodyFlags)

math.randomseed(os.time())

-- Pick one of the random seeds used in generating these tests
local testSeeds = { "670828913", "886885833", "715144259", "430281807", "543966758" }
local testSeed = testSeeds[math.random(1, #testSeeds)]
--local testSeed = testSeeds[1]


refResults = {
   ["newhist_model1"] = {
      ["100"] = {
         ["670828913"] = 1921.399019920596629,
         ["886885833"] = 2507.259567413587774,
         ["715144259"] = 2277.483069844485726,
         ["430281807"] = 1921.399019920596629,
         ["543966758"] = 1880.892858295289443
      },

      ["1024"] = {
         ["670828913"] = 4745.453858370087801,
         ["886885833"] = 5024.409732087724478,
         ["715144259"] = 5068.041460645472398,
         ["430281807"] = 5205.379022531156807,
         ["543966758"] = 5205.379022531156807
      },

      ["10000"] = {
         ["670828913"] = 12988.972447398160512,
         ["886885833"] = 13810.623768933184692,
         ["715144259"] = 13233.105481851351215,
         ["430281807"] = 13211.225057524350632,
         ["543966758"] = 13196.675140474046202
      }
   },

   ["newhist_model2"] = {
      ["100"] = {
         ["670828913"] = 1609.590270961241231,
         ["886885833"] = 1651.020703246448420,
         ["715144259"] = 1637.497574410253719,
         ["430281807"] = 1691.965879082630636,
         ["543966758"] = 1673.613780795002867
      },

      ["1024"] = {
         ["670828913"] = 2227.660833501358411,
         ["886885833"] = 2315.765891129145984,
         ["715144259"] = 2218.976540229048169,
         ["430281807"] = 2221.856526162277078,
         ["543966758"] = 2322.599989602272672
      },

      ["10000"] = {
         ["670828913"] = 5487.594694945784795,
         ["886885833"] = 5448.047330770828012,
         ["715144259"] = 5644.313735281130903,
         ["430281807"] = 5619.263388670427958,
         ["543966758"] = 5506.644390509201003
      }
   },

   ["newhist_model3"] = {
      ["100"] = {
         ["670828913"] = 20831.672944609566912,
         ["886885833"] = 4158.342887342659196,
         ["715144259"] = 20831.672944609566912,
         ["430281807"] = 20831.672944609566912,
         ["543966758"] = 20831.672944609566912
      },

      ["1024"] = {
         ["670828913"] = 19683.357772942861629,
         ["886885833"] = 19683.357772942861629,
         ["715144259"] = 20831.672944609566912,
         ["430281807"] = 20831.672944609566912,
         ["543966758"] = 19683.357772942861629
      },

      ["10000"] = {
         ["670828913"] = 20609.774305598719366,
         ["886885833"] = 20609.774305598719366,
         ["715144259"] = 20542.412593231350911,
         ["430281807"] = 20559.213177104018541,
         ["543966758"] = 20508.890734055767098
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


