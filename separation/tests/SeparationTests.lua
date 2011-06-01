--
-- Copyright (C) 2011 Matthew Arsenault
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

require "ResultSets"

argv = {...}

binName = argv[1]
extraFlags = argv[2] -- Extra flags to pass, such as choosing which SSE path to use
testName = argv[3]

assert(binName, "Binary name not set")




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


-- Find numbers in between xml tags called tagName
function findResults(str, tagName)
   assert(str, "Expected string to search for results")
   assert(tagName, "Expected tagName")

   local _, innerTag = str:match("<" .. tagName .. "%s*(.-)>(.-)</" .. tagName .. ">")
   local i, results = 1, { }

   assert(innerTag ~= nil, "Expected to find tag " .. innerTag)

   for num in innerTag:gmatch("[+-]?%d+[.]?%d+e?[+-]?%d+%s+") do
      results[i] = tonumber(num)
      i = i + 1
   end

   return results
end

-- Find the expected results from the printed output
function findSeparationResults(str)
   local bgInt, stInt, bgLike, stLike, searchLike
   local results = { }

   bgInt = findResults(str, "background_integral")
   stInt = findResults(str, "stream_integral")
   bgLike = findResults(str, "background_likelihood")
   stLike = findResults(str, "stream_only_likelihood")
   searchLike = findResults(str, "search_likelihood")

   assert(#bgInt == 1, "Expected to find one background_integral")
   assert(#bgLike == 1, "Expected to find one background_only_likelihood")
   assert(#searchLike == 1, "Expected to find one search_likelihood")

   results.background_integral = bgInt[1]
   results.stream_integral = stInt
   results.background_likelihood = bgLike[1]
   results.stream_only_likelihood = stLike
   results.search_likelihood = searchLike[1]

   return results
end


function resultCloseEnough(a, b)
   return math.abs(a - b) < 1.0e-12
end

local errFmtStr = [[
Result '%s' differs from expected:
   Actual = %20.15f   Expected = %20.15f   |Difference| = %20.15f

]]

function checkSeparationResults(results, reference)
   local function compareField(name)
      local closeEnough

      assert(results[name] ~= nil, string.format("Field '%s' not set for results", name))
      assert(reference[name] ~= nil, string.format("Field '%s' not set for reference", name))

      local resultField, refField = results[name], reference[name]

      assert(type(resultField) == type(refField), "Expected results and reference to be the same type")

      if type(refField) == "table" then
         closeEnough = true
         local closeTmp

         assert(#resultField == #refField and #resultField > 0, "Expected reference result table to have same > 0 size")


         -- Report all of the items which don't match rather than stopping on the first
         for i = 1, #refField do
            closeTmp = resultCloseEnough(resultField[i], refField[i])

            if not closeTmp then
               io.stderr:write(string.format(errFmtStr,
                                             string.format("%s\[%d\]", name, i - 1),
                                             resultField[i],
                                             refField[i],
                                             math.abs(resultField[i] - refField[i])
                                          )
                            )
            end

            closeEnough = closeEnough and closeTmp
         end
      elseif type(refField) == "number" then
         closeEnough = resultCloseEnough(results[name], reference[name])
         if not closeEnough then
            io.stderr:write(
               string.format(errFmtStr,
                             name,
                             results[name],
                             reference[name],
                             math.abs(reference[name] - results[name]))
            )
         end
      else
         error("Result must be a number or table of numbers")
      end

      return closeEnough
   end

   -- and with correct second to avoid short circuiting
   local correct = true

   correct = compareField("background_integral") and correct
   correct = compareField("stream_integral") and correct
   correct = compareField("background_likelihood") and correct
   correct = compareField("stream_only_likelihood") and correct
   correct = compareField("search_likelihood") and correct

   return correct
end

function runTest(test, checkResults)
   assert(test, "No test found!")
   local output

   if test.parameters ~= nil then
      output = os.readProcess(binName,
                              extraFlags,
                              "-i", -- FIXME: Avoid stale checkpoints a better way?
                              "-g",
                              "-a", test.file,
                              "-s", test.stars,
                              "-np", #test.parameters,
                              "-p", table.concat(test.parameters, " ")
                           )
   else
      output = os.readProcess(binName, extraFlags, "-i", "-g", "-a", test.file, "-s", test.stars)
   end

   print(output)

   -- If we aren't checking the reference results, then just print the
   -- output of the process.
   if not checkResults then
      return false
   end

   local results = findSeparationResults(output)
   assert(results ~= nil, "Results missing from output?")

   local check = checkSeparationResults(results, test.results)

   return check
end

rc = 0
for name, test in pairs(testSet) do
   io.stdout:write(string.rep("-", 80) .. "\n")
   io.stdout:write(string.format("Beginning test %s:\n\n", name))

   if runTest(test, false) then
      rc = 1
   end
   io.stdout:write(string.format("Test %s completed:\n", name))
   io.stdout:write(string.rep("-", 80) .. "\n")
end

os.exit(rc)



