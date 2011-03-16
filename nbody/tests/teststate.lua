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

--------------------------------------------------------------------------------

function getKeyAnswerPair(rawCtx)
   local hash = rawCtxHash(rawCtx)
   rawCtx.result = "I am a sample result"
   return { hash, rawCtx }
end

-- contextCombinations = buildAllCombinations(
--    function (theta, rsize, crit, useQuad, allowIncest)
--       local c = { }
--       c.theta = theta
--       c.treeRSize = rsize
--       c.criterion = crit
--       c.useQuad = useQuad
--       c.allowIncest = allowIncest
--       return c
--    end,
--    resultTable.thetas,
--    resultTable.treeRSizes,
--    resultTable.allCriterion,
--    resultTable.useQuads,
--    resultTable.allowIncests)






function makePotentialA()
      local disk = Disk.miyamotoNagai{ mass = 4.45865888e5,
                                    scaleLength = 6.5,
                                    scaleHeight = 0.26
                                  }

   local halo = Halo.logarithmic{ vhalo = 73,
                                  scaleLength = 12.0,
                                  flattenZ = 1.0
                                }

   local spherical = Spherical.spherical{ mass = 1.52954402e5,
                                          scale = 0.7
                                        }

   local pot = Potential.create{ disk = disk,
                                 halo = halo,
                                 spherical = spherical
                               }
   return pot
end

function makePotentialB()
   local potA = makePotentialA()
   -- Replace disk with exponential disk
   potA.disk = Disk.exponential{ scaleLength = 7,
                                 mass = 5.0e5
                              }
   return potA
end

testPotentials = { potentialA = makePotentialA(), potentialB = makePotentialB() }

function mergeModels(m1, ...)
   local i = #m1
   for _, m in ipairs({...}) do
      for _, v in ipairs(m) do
         i = i + 1
         m1[i] = v
      end
   end
   return m1
end


-- Each model returns (model, eps2, timestep)
testModels = {
   modelA =
      function(nbody, seed)
         local mod, eps2, dt
         mod = predefinedModels.plummer(nbody, {
                                           prng = DSFMT.create(seed),
                                           position = Vector.create(-22.0415, -3.35444, 19.9539),
                                           velocity = Vector.create(118.444, 168.874, -67.6378),
                                           mass = 16,
                                           scaleRadius = 0.2
                                        })
         eps2 = 0.01
         dt = 0.01
         return mod, eps2, dt
      end,

   modelB =
      function(nbody, seed)
         local eps2, pos, vel, prng, m1, m2
         pos = Vector.create(-22.0415, -3.35444, 19.9539)
         vel = Vector.create(118.444, 168.874, -67.6378)
         prng = DSFMT.create(seed)

         m1 = predefinedModels.plummer(nbody / 5, {
                                          prng = prng,
                                          position = pos,
                                          velocity = vel,
                                          mass = 12,
                                          scaleRadius = 0.2,
                                          ignore = true
                                       })

         m2 = predefinedModels.plummer(4 * nbody / 5, {
                                          prng = prng,
                                          position = pos,
                                          velocity = vel,
                                          mass = 12,
                                          scaleRadius = 0.5,
                                          ignore = true
                                       })

         eps2 = 0.01
         dt = 0.01
         return mergeModels(m1, m2), eps2, dt
      end
}

-- returns (ctx, st)
function getTestNBodyState(t)
   local ctx, pot, model, bodies
   pot = testPotentials[t.potential]
   model = testModels[t.model]
   bodies, eps2, dt = model(t.nbody, t.seed)

   ctx = NBodyCtx.create{ timestep = dt,
                       -- Irrelevant, tests aren't run by the C stuff but avoid the safety check
                          timeEvolve = 42.0,
                          theta = t.theta,
                          eps2 = eps2,
                          treeRSize = t.treeRSize,
                          criterion = t.criterion,
                          useQuad = t.useQuad,
                          allowIncest = t.allowIncest,
                       }
   return ctx, NBodyState.create(ctx, pot, bodies)
end

function runTest(t)
   local st, ctx, i
   ctx, st = getTestNBodyState(t)
   print("initial hash", st:hashBodies())

   for i = 1, t.nSteps do
      if st:step(ctx) then
         print("Error running step")
         break
      end
      print("hash step", i, st:hashBodies())
   end

   print("Final hash", st:hashBodies())
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


resultTable.hashtable = { }
for _, v in ipairs(tests) do
   resultTable.hashtable[hashNBodyTest(v)] = v
end

--table.foreach(hashTestTable, function(_, a) table.foreach(a, print) end)
--table.foreach(hashTestTable, function(k, _) print(k) end)


-- --persistence.store("result_table", resultTable)
-- loaded = persistence.load("result_table")
-- assert(loaded ~= nil, "Failed to load result_table")

-- print("arstarstarstarstrast")
-- table.foreach(tests[45], print)

-- lookupKey = hashNBodyTest(tests[45])
-- print("Lookup", lookupKey)
-- print("looked up", loaded.hashtable[lookupKey])

-- table.foreach(loaded.hashtable[lookupKey], print)

print("There are lots", #tests)
runTest(tests[9003])


