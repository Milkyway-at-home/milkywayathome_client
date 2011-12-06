--
-- Copyright (C) 2011  Matthew Arsenault
--
-- This file is part of Milkway@Home.
--
-- Milkyway@Home is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
--
-- Milkyway@Home is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
--
--


require "persistence"
require "curry"

-- Utility for helping get errors from printf/eprintf only with error reported at actual use site
local function stringFormatCatchError(...)
   local args = { ... }
   -- trap the errors from string.format so we can tell where the error actually comes from
   local success, str = pcall(function() return string.format(unpack(args)) end)
   if success then
      return str
   else
      io.stderr:write(str .. "\n") -- str is the error from string.format
      error(str, 3) -- call level is not here, and then above printf/eprintf
   end
end

function printf(...)
   io.stdout:write(stringFormatCatchError(...))
end

function eprintf(...)
   io.stderr:write(stringFormatCatchError(...))
end



-- Given a function of an arbitrary number of arguments, and a number
-- of lists corresponding to that number of arguments, returns that
-- function applied to every combination of arguments from the lists.
function buildAllCombinations(f, ...)
   -- Apply each function in fs to each value in tbl
   local function mapApply(fs, tbl)
      return foldl(function(acc, f)
                      return foldl(function(acc, v)
                                      acc[#acc + 1] = f(v)
                                      return acc
                                   end,
                                   acc,
                                   tbl)
                   end,
                   { },
                   fs)
   end

   return foldl(mapApply, { curry(f, #{...}) }, {...})
end


function getKeyNames(t)
   local keys = { }
   local i = 1
   for k, _ in pairs(t) do
      keys[i] = k
      i = i + 1
   end
   return keys
end

function printTable(tbl)
   assert(type(tbl) == "table", "Expected table to print")
   io.stderr:write("------------------------------\n" .. tostring(tbl) .. "\n")
   local function printTableMain(t, level)
      for k, v in pairs(t) do
         local prefix = ""
         for i = 0, level do
            prefix = prefix .. "  "
         end

         io.stderr:write(prefix .. tostring(k) .. " = " .. tostring(v) .. "\n")
         if type(v) == "table" then
            printTableMain(v, level + 1)
            io.stderr:write("\n")
         end
      end
   end
   printTableMain(tbl, 0)
   io.stderr:write("------------------------------\n")
end

function deepcopy(object)
   local lookup_table = {}
   local function _copy(object)
      if type(object) ~= "table" then
         return object
      elseif lookup_table[object] then
         return lookup_table[object]
      end
      local new_table = {}
      lookup_table[object] = new_table
      for index, value in pairs(object) do
         new_table[_copy(index)] = _copy(value)
      end
      return setmetatable(new_table, getmetatable(object))
   end
   return _copy(object)
end

function mergeTables(m1, ...)
   local i = #m1
   for _, m in ipairs({...}) do
      for _, v in ipairs(m) do
         i = i + 1
         m1[i] = v
      end
   end
   return m1
end

function runTest(t)
   local status, fatal = "NBODY_SUCCESS", false  -- 0 steps is OK

   local ctx, pot, st = getTestNBodyState(t)
   for i = 1, t.nSteps do
      status = st:step(ctx, pot)
      fatal = statusIsFatal(status)
      if fatal then
         break
      end
   end

   return st:hashSortBodies(), status, fatal
end


function runTestSet(tests, resultTable)
   resultTable.hashtable = { }

   local i = 1

   for _, t in ipairs(tests) do
      print("Running test ", i, 100 * i / #tests)
      printResult(t)
      local resultHash, status, failed = runTest(t)
      local testHash = hashNBodyTest(t)
      resultTable.hashtable[testHash] = t
      resultTable.hashtable[testHash].result = resultHash
      resultTable.hashtable[testHash].status = status
      resultTable.hashtable[testHash].failed = failed
      i = i + 1

   end

   return resultTable
end

function generateResultsToFile(tests, resultTable, file)
   resultTable = runTestSet(tests, resultTable)
   persistence.store(file, resultTable)
end

function loadResultsFromFile(file)
   return assert(persistence.load(file), "Failed to load results from file " .. file)
end

function findTestResult(result, resultTable)
   local key, foundResult
   assert(type(resultTable) == "table", "Expected result table to be a table")
   assert(type(result) == "table", "Expected result to look up to be a table")

   key = hashNBodyTest(result)
   foundResult = resultTable.hashtable[key]

   if foundResult == nil then
      io.stderr:write("Result for test not found in result table\n")
   end

   return foundResult
end


function printResult(t)
   local fmt =
[[
{
   theta       = %.15f,
   treeRSize   = %.15f,
   seed        = %u,
   nbody       = %u,
   allowIncest = %s,
   nSteps      = %u,
   criterion   = "%s",
   potential   = "%s",
   doublePrec  = %s,
   model       = "%s",
   useQuad     = %s
]]

   local str
   str = string.format(fmt,
                       t.theta,
                       t.treeRSize,
                       t.seed,
                       t.nbody,
                       tostring(t.allowIncest),
                       t.nSteps,
                       t.criterion,
                       t.potential,
                       tostring(t.doublePrec),
                       t.model,
                       tostring(t.useQuad))

   local resultFmt =
[[

   result      = "%s",
   status      = "%s",
   failed      = %s
]]

   if (t.result ~= nil) then
      str = str .. string.format(resultFmt, t.result, tostring(t.status), tostring(t.failed))
   end

   str = str .. "}\n"

   print(str)
end

function checkTestResult(test, resultTable)
   local expected = findTestResult(test, resultTable)
   if expected == nil then
      error("Test result not in result table")
   end

   assert(expected.result ~= nil, "Hash missing from expected result")
   assert(expected.status ~= nil, "Status missing from expected result")
   assert(expected.failed ~= nil, "Failed status missing from expected result")

   local doesNotMatch = false
   local hash, status, failed = runTest(test)

   if hash ~= expected.result then
      io.stderr:write("Hash does not match expected:\n")
      doesNotMatch = true
   end

   if status ~= expected.status then
      io.stderr:write("Status does not match expected:\n")
      doesNotMatch = true
   end

   if failed ~= expected.failed then
      io.stderr:write("Failed status does not match expected:\n")
      doesNotMatch = true
   end

   if doesNotMatch then
      io.stderr:write("Failing test:\n")
      io.stderr:write(string.format("Got hash = %s, status = %s, failed = %s\nExpected:\n",
                                    hash, status, tostring(failed)))
      printResult(test)
   end

   return doesNotMatch
end

-- Marshal between a table with "x" "y" and "z" and the userdata
-- Vector type since persistence can't handle userdata
function vectorToTable(v)
   local vTable = { }
   vTable.x = v.x
   vTable.y = v.y
   vTable.z = v.z
   return vTable
end

function tableToVector(vTable)
   return Vector.create(vTable.x, vTable.y, vTable.z)
end


function os.readProcess(bin, ...)
   local args, cmd
   args = table.concat({...}, " ")

   -- Redirect stderr to stdout, since popen only gets stdout
   cmd = table.concat({ bin, args, "2>&1" }, " ")
   local f, rc = assert(io.popen(cmd, "r"))
   local s = assert(f:read('*a'))
   f:close()
   return s, rc
end

function getExtraNBodyFlags()
   nbodyFlags = os.getenv("NBODY_FLAGS")
   if nbodyFlags == nil then
      nbodyFlags = ""
   end
   return nbodyFlags
end


function calcMean(table)
   local total = 0.0
   assert(type(table) == "table")

   if #table == 0 then
      return 0.0
   end

   for k, v in ipairs(table) do
      total = total + v
   end

   return total / #table
end

function calcStddev(table, mean)
   local total = 0.0

   -- Bad things happen if this ends up nan
   if #table == 0 then
      return 0.0
   end


   if mean == nil then
      mean = calcMean(table)
   end

   for k, v in ipairs(table) do
      total = total + sqr(v - mean)
   end

   return sqrt(total / (#table - 1))
end

function calcStats(table)
   local avg = calcMean(table)
   return avg, calcStddev(table, avg)
end



function randomSeed()
   return math.random(0, 32767)
end

function findMin(table)
   assert(type(table) == "table" and #table > 0)
   local low = table[1]
   for _, i in ipairs(table) do
      if i < low then
         low = i
      end
   end

   return low
end


function findMinMax(table)
   assert(type(table) == "table" and #table > 0)

   local low = table[1]
   local high = table[1]

   for _, i in ipairs(table) do
      if i < low then
         low = i
      end

      if i > high then
         high = i
      end
   end

   return low, high
end

function runSimple(arg)
   return os.readProcess(arg.nbodyBin or "milkyway_nbody",
                         "--checkpoint-interval=-1",
                         "--debug-boinc",
                         "-t",
                         "-f", arg.input,
                         table.concat(arg.extraArgs, " ")
                      )
end

function runFullTest(arg)
   local testPath, histogramPath = arg.testName, arg.histogram
   if arg.testDir then
      testPath = string.format("%s/%s.lua", arg.testDir, arg.testName)
      histogramPath = string.format("%s/%s", arg.testDir, arg.histogram)
   else
      testPath = arg.testName .. ".lua"
   end

   mkdir(".checkpoints")

   local cpFlag = string.format("--no-clean-checkpoint -c %s/%s__%d",
                                ".checkpoints",
                                arg.testName,
                                arg.seed)
   if not arg.cached then
      cpFlag = cpFlag .. " --ignore-checkpoint"
   end

   return os.readProcess(arg.nbodyBin or "milkyway_nbody",
                         "--checkpoint-interval=-1", -- Disable checkpointing
                         "-g", -- Prevent stderr from getting consumed with BOINC
                         "-t",
                         "-f", testPath,
                         "-h", histogramPath,
                         "--seed", arg.seed,
                         cpFlag,
                         getExtraNBodyFlags(),
                         table.concat(arg.extraArgs, " ")
                      )

end

-- Find the likelihood from the output of the process
function findNumber(str, name)
   local pattern = string.format("<%s>(.+)</%s>", name, name)
   local m = str:match(pattern)
   local lineSep = string.rep("-", 80) .. "\n"
   if m == nil then
      eprintf("Didn't match '%s' in output\nOffending output:\n%s\n%s%s",
              name,
              lineSep,
              str,
              lineSep
           )
      return nil
   else
      return tonumber(m)
   end
end

function findLikelihood(str, emd)
   if emd then
      name = "emd"
   else
      name = "search_likelihood"
   end

   -- likelihood is negated
   return -findNumber(str, name)
end



