--
-- Copyright (C) 2011  Matthew Arsenault
--
-- This file is part of Milkway@Home.
--
-- Milkyway@Home is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.

-- Milkyway@Home is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
--
--

require "NBodyTesting"
SP = require "SamplePotentials"
SM = require "SampleModels"
require "persistence"

local generatingResults = false



local function comboTests()
   local tests = SP.makeAllPotentials()
   print("Number combinations", #tests)

   local results = { }
   local prng = DSFMT.create()

   local i = 1
   local r
   for _, v in ipairs(tests) do
      r = prng:randomVector(50)
      results[i] = { ["potential"] = v, ["r"] = r, ["accel"] = v:acceleration(r) }
      i = i + 1
   end

   return results
end

local function randomTests()
   local nTests = 1000000
   local prng = DSFMT.create()
   local p
   tests = { }
   for i = 1, nTests do
      p = SP.randomPotential(prng)
      r = prng:randomVector(50)
      tests[i] = { ["potential"] = p, ["position"] = r, ["acceleration"] = p:acceleration(r) }
      printTable(tests[i])
      i = i + 1
   end

   return tests
end

-- FIXME: Storing these userdatas doesn't work
function generatePotentialComboTests()
   local results = comboTests()
   persistence.store("potential_test", results)
end

function checkPotentialTests()
   local results = persistence.load("potential_test")
   local r, p

   for _, v in ipairs(results) do
      if v.potential:acceleration(v.r) ~= v.acceleration then
         io.stderr:write("Error in acceleration:\nExpected:\n")
         printTable(v)
         io.stderr:write(string.format("Got acceleration = %20.15f\n", v.acceleration))
         os.exit(1)
      end
   end
end


if generatingResults then
   generatePotentialComboTests()
else
   checkPotentialTests()
end


