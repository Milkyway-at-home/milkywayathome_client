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

require "nbody_testing"
require "sample_potentials"
require "sample_models"

testModels = sampleModels  -- FIXME
testPotentials = samplePotentials
-- returns (ctx, st)
function getTestNBodyState(t)
   local ctx, pot, model, bodies
   pot = testPotentials[t.potential]
   model = testModels[t.model]
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

resultTable = { potentials   = samplePotentialNames,
                models       = sampleModelNames,
                nbody        = { 100, 200 },
                nSteps       = { 0, 1, 3, 10 },
                seeds        = { 1234567890, 609746760, 1000198000 },
                thetas       = { 1.0, 0.9, 0.7, 0.5, 0.3 },
                treeRSizes   = { 4.0, 8.0, 2.0, 1.0 },
                criterion    = { "sw93", "NewCriterion", "BH86", "Exact" },
                useQuads     = { true, false },
                allowIncests = { true, false }
             }


tests = buildAllCombinations(
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



smallerList = { }
j = 1
for i = 1, 38000, 1000 do
   smallerList[j] = tests[i]
   j = j + 1
end

--resultTable = generateTestResults(tests, resultTable)
--printTable(resultTable)


generateResultsToFile(smallerList, { }, "small_results")
smallLoad = loadResultsFromFile("small_results")
printTable(smallLoad)



brokenTest = {
   theta = 0.7,
   treeRSize = 4,
   allowIncest = true,
   useQuad = true,
   criterion = "NewCriterion",
   model = "modelB",
   nSteps = 1,
   nbody = 100,
   potential = "potentialA",
   seed = 609746760
}

sampleLookup = {
   theta = 1,
   treeRSize = 2,
   seed = 1234567890,
   nbody = 100,
   model = "modelB",
   allowIncest = true,
   nSteps = 1,
   criterion = "BH86",
   potential = "potentialA",
   doublePrec = true,
   useQuad = true
}

printTable(findTestResult(sampleLookup, smallLoad))

printResult(sampleLookup)

print("-------------")
printResult(findTestResult(sampleLookup, smallLoad))


