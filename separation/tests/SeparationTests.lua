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


argv = {...}

binName = argv[1]
testName = argv[2]

assert(binName, "Binary name not set")


param_set1 = { 0.571713, 12.312119, -3.305187, 148.010257, 22.453902, 0.42035,
               -0.468858, 0.760579, -1.361644, 177.884238, 23.882892, 1.210639,
               -1.611974, 8.534378, -1.361644, 177.884238, 10.882892, 1.210639,
               -1.611974, 8.534378
         }


param_set2 = { 0.571713, 12.312119, -3.305187, 148.010257, 22.453902, 0.42035, -0.468858,
               0.760520, -1.361644, 177.884238, 23.882892, 1.210639, -1.611974, 8.534378
         }



param_set3 = { 0.34217373320392042, 25.9517910846623, -2.1709414738826602, 38.272511356953906,
               30.225190442596112, 2.2149060013372885, 0.32316169064291655, 2.7740244716285285
            }

param_set4 = { 0.40587961154742185, 17.529961843393409, -1.8575145272144837, 29.360893891378243,
               31.228263575178566, -1.551741065334, 0.064096152599308373, 2.55428209912781
            }


param_set5 = { 0.73317163557524425, 14.657212876628332, -1.7054653473950408, 16.911711745343633,
               28.077212666463502, -1.2032908515814611, 3.5273606439247287, 2.2248214505875008
            }



testSet = {
   ["small_test11"] = {
      file       = "astronomy_parameters-11-small.txt",
      stars      = "stars-11.txt",
      parameters = param_set1,

      results = {
         background_integral = 0.000344404677767,
         stream_integral     = { 34.162598012270841, 667.258968797279977, 336.167834127368394 },

         background_likelihood  = -3.308234992702566,
         stream_only_likelihood = { -160.798151443575989, -4.185806766566943, -4.119821303303995 },

         search_likelihood = -3.125097883803232
      }
   },

   ["test12"] = {
      file       = "astronomy_parameters-12.txt",
      stars      = "stars-12.txt",
      parameters = param_set1,

      results = nil
   }

}

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

function runTest(testName, checkResults)
   local test = testSet[testName]
   assert(test, "Test " .. testName .. " not found!")

   io.stderr:write(string.rep("-", 80) .. "\n")
   io.stderr:write(string.format("Beginning test %s:\n", testName))

   local output = os.readProcess(binName,
                                 "-i", -- FIXME: Avoid stale checkpoints a better way?
                                 "-g",
                                 "-a", test.file,
                                 "-s", test.stars,
                                 "-np", #test.parameters,
                                 "-p", table.concat(test.parameters, " ")
                              )

   -- If we aren't checking the reference results, then just print the
   -- output of the process.
   if not checkResults then
      print(string.format("Test %s output:", testName))
      print(output)
      return false
   end

   local results = findSeparationResults(output)
   assert(results ~= nil, "Results missing from output?")

   local check = checkSeparationResults(results, test.results)

   io.stderr:write(string.format("Test %s completed:\n", testName))
   io.stderr:write(string.rep("-", 80) .. "\n")

   return check
end


if not runTest("small_test11", true) then
   os.exit(1)
end



