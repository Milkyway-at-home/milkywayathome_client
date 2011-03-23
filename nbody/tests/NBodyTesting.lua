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

