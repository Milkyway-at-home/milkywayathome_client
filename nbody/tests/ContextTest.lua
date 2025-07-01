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

local generatingResults = false

-- returns (ctx, st)
function getTestNBodyState(t)
   local ctx, potential, model, bodies, st
   pot = SP.samplePotentials[t.potential]
   model = SM.sampleModels[t.model]
   bodies, eps2, dt = model(t.nbody, t.seed)
   local prng = DSFMT.create(t.seed)

   ctx = NBodyCtx.create{
      timestep    = dt,
      timeEvolve  = 42.0,     -- Irrelevant, tests aren't run by the C stuff but avoid the safety check
      theta       = t.theta,
      eps2        = eps2,
      b           = 53.5,
      r           = 28.6,
      vx          = -156,
      vy          = 79,
      vz          = 107,
      treeRSize   = t.treeRSize,
      criterion   = t.criterion,
      useQuad     = t.useQuad,
      BestLikeStart = 0.95,
      BetaSigma     = 2.5,
      VelSigma      = 2.5,
      DistSigma     = 2.5,
      PMSigma       = 2.5,
      BetaCorrect   = 1.111,
      VelCorrect    = 1.111,
      DistCorrect   = 1.111,
      PMCorrect     = 1.111,
      IterMax       = 6,
      allowIncest = t.allowIncest,
      quietErrors = true,
      LMC         = t.LMC,
      LMCfunction = 1,
      LMCmass     = 449865.888,
      LMCscale    = 15.0,
      LMCscale2   = 16.6,
      LMCDynaFric = t.LMCDynaFric
   }
   --Add potential to context
   ctx:addPotential(pot)

   if (ctx.LMC) then
      st = NBodyState.createRandomLMC(ctx, 51, prng, bodies)
   else
       st = NBodyState.create(ctx, bodies)
   end

   return ctx, st
end

local resultTable = {
   potentials   = SP.samplePotentialNames,
   models       = SM.sampleModelNames,
   nbody        = { 100, 1024 },
   nSteps       = { 1, 4 },
   seeds        = { 1234567890, 609746760 },
   thetas       = { 1.0, 0.9, 0.5, 0.3 },
   treeRSizes   = { 8.0, 4.0, 2.0, 1.0 },
   criterion    = { "SW93", "TreeCode", "BH86", "Exact" },
   useQuads     = { true, false },
   allowIncests = { true },  -- Might as well allow it for the tests.
   LMCs         = { true, false },
   LMCDynaFrics = { true }   -- Setting to true to save time
}

-- Get list of all tests
local function generateFullTestSet()
   return buildAllCombinations(
      function(potential, model, nbody, nSteps, seed, theta, rsize, crit, useQuad, allowIncest, LMC, LMCDynaFric)
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
         c.LMC         = LMC
         c.LMCDynaFric = LMCDynaFric
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
      resultTable.allowIncests,
      resultTable.LMCs,
      resultTable.LMCDynaFrics)
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




