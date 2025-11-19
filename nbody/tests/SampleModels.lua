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

local SampleModels = { }


-- Each model returns (model, eps2, timestep)
SampleModels.sampleModels = {
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
         return mergeTables(m1, m2), eps2, dt
      end
}

SampleModels.sampleModelNames = getKeyNames(SampleModels.sampleModels)

SampleModels.randomPlummer = function(prng, maxN)
   if prng == nil then
      prng = DSFMT.create()
   end

   if maxN == nil then
      maxN = 10000
   end

   return predefinedModels.plummer{
      nbody       = prng:random(10, maxN),
      prng        = prng,
      position    = prng:randomVector(10),
      velocity    = prng:randomVector(100),
      mass        = prng:random(0.1, 50),
      scaleRadius = prng:random(0.01, 10),
      ignore      = prng:randomBool()
   }
end


return SampleModels


