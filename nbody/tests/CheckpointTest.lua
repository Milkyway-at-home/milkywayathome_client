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

require "NBodyTesting"
SM = require "SampleModels"
SP = require "SamplePotentials"

function randomNBodyCtx(prng)
   if prng == nil then
      prng = DSFMT.create()
   end

   return NBodyCtx.create{
      timestep    = prng:random(1.0e-5, 1.0e-4),
      timeEvolve  = prng:random(0, 10),
      theta       = prng:random(0, 1),
      eps2        = prng:random(1.0e-9, 1.0e-3),
      treeRSize   = prng:randomListItem({ 4, 8, 2, 16 }),
      criterion   = prng:randomListItem({"NewCriterion", "SW93", "BH86", "Exact"}),
      useQuad     = prng:randomBool(),
      allowIncest = true,
      quietErrors = true
   }
end

function runNSteps(st, n, ctx)
   for i = 1, n do
      st:step(ctx)
   end
   return ctx, st
end

function runInterruptedSteps(st, totalSteps, ctx, prng)
   local tmpDir = os.getenv("TMP") or ""
   local checkpoint = tmpDir .. os.tmpname()

   for i = 1, totalSteps do
      st:step(ctx)
      if prng:randomBool() then
         local tmp = tmpDir .. os.tmpname()
         st:writeCheckpoint(ctx, checkpoint, tmp)
         ctx, st = NBodyState.readCheckpoint(checkpoint)
         os.remove(checkpoint)
      end
   end

   return ctx, st
end


local nTests = 5

for i = 1, nTests do
   local testSteps, st, stCopy
   local ctx, m
   local prng = DSFMT.create()

   m = SM.randomPlummer(prng, 500)
   ctx = randomNBodyCtx(prng)
   ctx:addPotential(SP.randomPotential(prng))


   st = NBodyState.create(ctx, m)
   stClone = st:clone()

   testSteps = floor(prng:random(0, 51))

   ctx, st = runNSteps(st, testSteps, ctx)
   ctxClone, stClone = runInterruptedSteps(stClone, testSteps, ctx, prng)

   assert(ctx == ctxClone,
          string.format("Checkpointed context does not match:\nctx 1 = %s\n ctx 2 = %s\n",
                        tostring(ctx),
                        tostring(ctxClone))
       )

   assert(st == stClone,
          string.format("Checkpointed state does not match:\nstate 1 = %s\n state 2 = %s\n",
                        tostring(st),
                        tostring(stClone))
       )
end


