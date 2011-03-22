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

local generatingResults = true

-- returns (ctx, potential, st)
function getTestNBodyState(t)
   local ctx, pot, model, bodies
   pot = SP.samplePotentials[t.potential]
   model = SM.sampleModels[t.model]
   bodies, eps2, dt = model(t.nbody, t.seed)

   ctx = NBodyCtx.create{
      timestep    = dt,
      timeEvolve  = 42.0,     -- Irrelevant, tests aren't run by the C stuff but avoid the safety check
      theta       = t.theta,
      eps2        = eps2,
      treeRSize   = t.treeRSize,
      criterion   = t.criterion,
      useQuad     = t.useQuad,
      allowIncest = t.allowIncest,
      quietErrors = true
   }
   return ctx, pot, NBodyState.create(ctx, pot, bodies)
end

local resultTable = {
   potentials   = SP.samplePotentialNames,
   models       = SM.sampleModelNames,
   nbody        = { 100, 1024 },
   nSteps       = { 1, 3, 8 },
   seeds        = { 1234567890, 609746760, 1000198000 },
   thetas       = { 1.0, 0.9, 0.5, 0.3 },
   treeRSizes   = { 8.0, 4.0, 2.0, 1.0 },
   criterion    = { "SW93", "NewCriterion", "BH86", "Exact" },
   useQuads     = { true, false },
   allowIncests = { true }  -- Might as well allow it for the tests.
}

-- Get list of all tests
local function generateFullTestSet()
   return buildAllCombinations(
      function(potential, model, nbody, nSteps, seed, theta, rsize, crit, useQuad, allowIncest)
         local c = { }
         c.doublePrec  = true
         c.potential   = potential
         c.model       = model
         c.nbody       = nbody
         c.nSteps      = nSteps
         c.seed        = seed
         c.theta       = theta
         c.treeRSize   = rsize
         c.criterion   = crit
         c.useQuad     = useQuad
         c.allowIncest = allowIncest
         return c
      end,
      resultTable.potentials,
      resultTable.models,
      resultTable.nbody,
      resultTable.nSteps,
      resultTable.seeds,
      resultTable.thetas,
      resultTable.treeRSizes,
      resultTable.criterion,
      resultTable.useQuads,
      resultTable.allowIncests)
end

if generatingResults then
   local fullTests = generateFullTestSet()
   print("Running ", #fullTests)
   generateResultsToFile(fullTests, resultTable, "context_test_results")
else
   local set = generateFullTestSet()
   local refTable = loadResultsFromFile("context_test_results")

   local i, n = 1, #set
   for _, t in ipairs(set) do
      print(100 * i / n)
      checkTestResult(t, refTable)
      i = i + 1
   end
end




