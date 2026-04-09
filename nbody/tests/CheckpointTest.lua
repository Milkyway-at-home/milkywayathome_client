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

function erf(x)       --Pulled from https://hewgill.com/picomath/lua/erf.lua.html
    -- constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    -- Save the sign of x
    sign = 1
    if x < 0 then
        sign = -1
    end
    x = math.abs(x)

    -- A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)

    return sign*y
end

dec = 9.0   -- -- number of decimals to round to (default: 9.0)
function round(num, places)
    local mult = 10.0^(places)
    return floor(num * mult + 0.5) / mult
  end
rscale_l            = {round( 2.9,     dec)}  -- Baryonic Radius (kpc)
light_r_ratio       = {round( 0.2,     dec)}  -- Baryonic Radius / (Baryonic Radius + Dark Matter Radius)
mass_l              = {round( 2429.198,dec)}  -- Baryonic Mass (Structure Mass Units)
light_mass_ratio    = {round( 0.0830,  dec)}  -- Baryonic Mass / (Baryonic Mass + Dark Matter Mass)
orbit_parameter_l   = {round( 302.801, dec)}  -- Galactocentric l
orbit_parameter_b   = {round( -44.328, dec)}  -- Galactocentric b
orbit_parameter_r   = {round( 62.4,    dec)}  -- Galactocentric r
orbit_parameter_vx  = {round( 21.99,   dec)}  -- Galactocentric vx
orbit_parameter_vy  = {round( -201.36, dec)}  -- Galactocentric vy
orbit_parameter_vz  = {round( 171.25,  dec)}  -- Galactocentric vz

function randomNBodyCtx(prng)
   if prng == nil then
      prng = DSFMT.create()
   end
   sigma = prng:random(1.5,3.0)
   correct = math.sqrt(2*3.1415926535)/(math.sqrt(2*3.1415926535)*erf(sigma/math.sqrt(2)) - 2*sigma*math.exp(-sigma*sigma/2))
   return NBodyCtx.create{
      dwarfn = 1,
      timestep      = prng:random(1.0e-5, 1.0e-4),
      timeEvolve    = prng:random(0, 10),
      theta         = prng:random(0, 1),
      eps2          = prng:random(1.0e-9, 1.0e-3),
      b           = orbit_parameter_b,
      r           = orbit_parameter_r,
      vx          = orbit_parameter_vx,
      vy          = orbit_parameter_vy,
      vz          = orbit_parameter_vz,
      -- b             = prng:random(40.0,60.0),
      -- r             = prng:random(10.0,30.0),
      -- vx            = prng:random(-200.0,200.0),
      -- vy            = prng:random(-200.0,200.0),
      -- vz            = prng:random(-200.0,200.0),
      treeRSize     = prng:randomListItem({ 4, 8, 2, 16 }),
      criterion     = prng:randomListItem({"TreeCode", "SW93", "BH86", "Exact"}),
      useQuad       = prng:randomBool(),
      BestLikeStart = prng:random(0.85,0.99),
      BetaSigma     = sigma,
      VelSigma      = sigma,
      DistSigma     = sigma,
      PMSigma       = sigma,
      MomentumSigma = sigma,
      BetaCorrect   = correct,
      VelCorrect    = correct,
      DistCorrect   = correct,
      PMCorrect     = correct,
      MomentumCorrect = correct,
      IterMax       = prng:randomListItem({ 2, 3, 4, 5 }),
      allowIncest   = true,
      quietErrors   = true,
      LMC           = prng:randomBool(),
      LMCfunction   = 1,
      LMCmass       = prng:random(1.0e5,1.0e6),
      LMCscale      = prng:random(1.0,20.0),
      LMCscale2     = prng:random(1.0,20.0),
      LMCDynaFric   = prng:randomBool(),
      coulomb_log   = 15
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
   local checks = 0

   for i = 1, totalSteps do
      st:step(ctx)
      if prng:randomBool() then
         checks = checks + 1
         local tmp = tmpDir .. os.tmpname()
         st:writeCheckpoint(ctx, checkpoint, tmp)
         ctx, st = NBodyState.readCheckpoint(checkpoint)
         os.remove(checkpoint)
      end
   end

   return ctx, st, checks
end


local nTests = 20
local prng = DSFMT.create(20251003)

for i = 1, nTests do
   local testSteps, st, stCopy
   local ctx, m
   local lmcpos, lmcvel

   m = SM.randomPlummer(prng, 500)
   ctx = randomNBodyCtx(prng)
   ctx:addPotential(SP.randomPotential(prng))

   testSteps = floor(prng:random(0, 51))

   if (ctx.LMC) then
      st = NBodyState.createRandomLMC(ctx, testSteps + 1, prng, m)
   else
      st = NBodyState.create(ctx, m)
   end

   stClone = st:clone()

   assert(st == stClone,
          string.format("Failed to properly clone state:\nstate 1 = %s\nstate 2 = %s\n",
                        tostring(st),
                        tostring(stClone))
       )

   ctx, st = runNSteps(st, testSteps, ctx)
   ctxClone, stClone, nCheckpoints = runInterruptedSteps(stClone, testSteps, ctx, prng)

   assert(ctx == ctxClone,
          string.format("Checkpointed context does not match after %d checkpoints:\nctx 1 = %s\nctx 2 = %s\n",
                        nCheckpoints,
                        tostring(ctx),
                        tostring(ctxClone))
       )

   assert(st == stClone,
          string.format("Checkpointed state does not match after %d checkpoints:\nstate 1 = %s\nstate 2 = %s\n%s",
                        nCheckpoints,
                        tostring(st),
                        tostring(stClone),
                        tostring(ctx))
       )
end


