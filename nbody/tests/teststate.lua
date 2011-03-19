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

function makePotentialA()
   local disk, halo, spherical
   disk = Disk.miyamotoNagai{
      mass        = 4.45865888e5,
      scaleLength = 6.5,
      scaleHeight = 0.26
   }

   halo = Halo.logarithmic{
      vhalo       = 73,
      scaleLength = 12.0,
      flattenZ    = 1.0
   }

   spherical = Spherical.spherical{
      mass  = 1.52954402e5,
      scale = 0.7
   }

   return Potential.create{
      disk      = disk,
      halo      = halo,
      spherical = spherical
   }
end

function makePotentialB()
   local potA = deepcopy(makePotentialA())
   -- Replace disk with exponential disk
   potA.disk = Disk.exponential{
      scaleLength = 7,
      mass        = 5.0e5
   }
   return potA
end

testPotentials = { potentialA = makePotentialA(), potentialB = makePotentialB() }

-- Each model returns (model, eps2, timestep)
testModels = {
   modelA =
      function(nbody, seed)
         local mod, eps2, dt
         local r0, mass = 0.2, 16
         mod = predefinedModels.plummer{
            nbody = nbody,
            prng = DSFMT.create(seed),
            position = Vector.create(-22.0415, -3.35444, 19.9539),
            velocity = Vector.create(118.444, 168.874, -67.6378),
            mass = mass,
            scaleRadius = r0
         }

         return mod, calculateEps2(nbody, r0), calculateTimestep(mass, r0)
      end,

   modelB =
      function(nbody, seed)
         local eps2, pos, vel, prng, m1, m2
         local smallR0, bigR0 = 0.2, 0.5
         local smallMass, bigMass = 12, 12
         pos = Vector.create(-22.0415, -3.35444, 19.9539)
         vel = Vector.create(118.444, 168.874, -67.6378)
         prng = DSFMT.create(seed)

         m1 = predefinedModels.plummer{
            nbody = nbody / 5,
            prng = prng,
            position = pos,
            velocity = vel,
            mass = smallMass,
            scaleRadius = smallR0,
            ignore = true
         }

         m2 = predefinedModels.plummer{
            nbody = 4 * nbody / 5,
            prng = prng,
            position = pos,
            velocity = vel,
            mass = bigMass,
            scaleRadius = bigR0,
            ignore = true
         }

         eps2 = calculateEps2(nbody, smallR0)
         dt   = calculateTimestep(smallMass + bigMass, smallR0)
         return mergeModels(m1, m2), eps2, dt
      end
}

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

function runTest(t)
   local st, ctx, pot
   local err = false
   ctx, pot, st = getTestNBodyState(t)

   for i = 1, t.nSteps do
      if st:step(ctx, pot) then
         err = true
         break
      end
   end

   return st:hashBodies(), err
end

resultTable = { potentials   = getKeyNames(testPotentials),
                models       = getKeyNames(testModels),
                nbody        = { 100, 200 },
                nSteps       = { 0, 1, 3, 10, 100 },
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

resultTable = generateTestResults(tests, resultTable)
printTable(resultTable)


--generateResultsToFile(smallerList, { }, "small_results")

--smallLoad = loadResultsFromFile("small_results")

--printTable(smallLoad)



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

