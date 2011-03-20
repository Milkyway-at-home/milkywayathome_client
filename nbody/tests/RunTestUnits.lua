
arg = {...}

assert(#arg == 4, "Test driver expected 4 arguments got " .. #arg)

nbodyBinary  = arg[1]
testName = arg[2]
testSeed = arg[3]
testBodies = arg[4]


--testDir = "tests"
testDir = "."

refResults = {
   ["orphan_model_1"] = {
      ["100"] = {
         ["670828913"] = -13833.731309987899294,
         ["886885833"] = 4,
         ["715144259"] = 4,
         ["430281807"] = 4,
         ["543966758"] = 4
      },

      ["1024"] = {
         ["670828913"] = 4,
         ["886885833"] = 4,
         ["715144259"] = 4,
         ["430281807"] = 4,
         ["543966758"] = 4
      }
   },

   ["orphan_model_2"] = {

   },

   ["orphan_model_3"] = {

   },

   ["orphan_model_4"] = {

   },

   ["orphan_model_5"] = {

   },

   ["orphan_model_6"] = {

   },

   ["orphan_model_7"] = {

   }
}

function getTestDir(baseName)
   return testDir .. "/" .. baseName .. "/"
end

function getTestFileName(baseName)
   return getTestDir(baseName) .. baseName .. ".lua"
end

function getTestHistogramName(baseName)
   return getTestDir(baseName) .. "histogram"
end

function os.readProcess(bin, ...)
   local args, cmd
   args = table.concat({...}, " ")
   -- Redirect stderr to stdout, since popen only gets stdout
   cmd = table.concat({ bin, args, "2>&1" }, " ")
   local f = assert(io.popen(cmd, "r"))
   local s = assert(f:read('*a'))
   f:close()
   return s
end


function runFullTest(testName, seed, ...)
   return os.readProcess(nbodyBinary,
                         "-i",
                         "--checkpoint-interval=-1",
                         "-t",
                         "-f", getTestFileName(testName),
                         "-h", getTestHistogramName(testName),
                         "--seed", seed,
                         table.concat({...}, " ")
                      )
end

-- Find the likelihood from the output of the process
function findLikelihood(str)
   local m = str:match("<search_likelihood>([+-]?%d+[.]?%d+e?[+-]?%d+)</search_likelihood>")
   local lineSep = string.rep("-", 80) .. "\n"
   if m == nil then
      io.stderr:write("Didn't match likelihood in output\nOffending output:\n" .. lineSep)
      io.stderr:write(str .. lineSep)
      return nil
   else
      return tonumber(m)
   end
end

function resultCloseEnough(a, b)
   return math.abs(a - b) < 1.0e-12
end

errFmtStr = [[
Result differs from expected:
   Expected = %20.15f  Actual = %20.15f  |Difference| = %20.15f
]]

function runCheckTest(testName, seed, nbody, ...)
   local fileResults = assert(refResults[testName], "Didn't find result for test file")
   local bodyResults = assert(fileResults[nbody], "Didn't find result with matching bodies")
   local refResult = assert(bodyResults[seed], "Didn't find result with matching seed")

   local ret = runFullTest(testName, seed, nbody, ...)
   local result = findLikelihood(ret)

   if result == nil then
      return false
   end

   local close = resultCloseEnough(refResult, result)
   if not close then
      io.stderr:write(string.format(errFmtStr, refResult, result, math.abs(result - refResult)))
   else
      print("close enough")
   end

   return close
end

if runCheckTest(testName, testSeed, testBodies) then
   os.exit(0)
else
   os.exit(1)
end



