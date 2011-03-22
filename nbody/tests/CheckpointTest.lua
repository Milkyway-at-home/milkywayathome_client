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

function runNSteps(st, n, ctx, pot)
   for i = 1, n do
      st:step(ctx, pot)
   end

   return ctx, st
end

function runInterruptedSteps(st, totalSteps, ctx, pot, prng)
   local checkpoint = os.tmpname()
   for i = 1, totalSteps do
      st:step(ctx, pot)
      if prng:randomBool() then
         st:writeCheckpoint(ctx, checkpoint, os.tmpname())
         ctx, st = NBodyState.readCheckpoint(checkpoint)
      end
   end

   os.remove(checkpoint)

   return ctx, st
end

local nTests = 2

for i = 1, nTests do
   local testSteps, st, stCopy
   local ctx, pot, m
   local prng = DSFMT.create()

   ctx = randomNBodyCtx(prng)
   pot = SP.randomPotential(prng)
   ctx.potential = pot  -- Kind of dumb hack so I don't have to fix separate potential in everything else
   m = SM.randomPlummer(prng, 500)

   st = NBodyState.create(ctx, pot, m)

   stClone = st:clone()
   testSteps = round(prng:random(0, 50))

   ctx, st = runNSteps(st, testSteps, ctx, pot)
   ctxClone, stClone = runInterruptedSteps(stClone, testSteps, ctx, pot, prng)

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


