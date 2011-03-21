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

function buildAllDisks()
   local miyamotoNagaiDisks, exponentialDisks
   miyamotoNagaiDisks = buildAllCombinations(
      function(m, a, b)
         return Disk.miyamotoNagai{ mass = m, scaleLength = a, scaleHeight = b }
      end,
      { 4.5e5, 5.0e6, 3.0e4 },
      { 6.5, 9.0, 3.0 },
      { 0.26, 0.5, 0.1 }
   )

   exponentialDisks = buildAllCombinations(
      function(m, b)
         return Disk.exponential{ mass = m, scaleLength = b }
      end,
      { 2.24933e5, 3.0e6, 1.5e4 },
      { 4, 7, 3 }
   )

   return mergeTables(miyamotoNagaiDisks, exponentialDisks)
end

function buildAllSphericals()
   return buildAllCombinations(
      function(m, s)
         return Spherical.spherical{ mass = m, scale = s }
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 },
      { 0.7, 1.1, 0.4 }
   )
end


function buildAllHalos()
   local logHalos, nfwHalos, triaxialHalos

   logHalos = buildAllCombinations(
      function(vhalo, scaleLength, flattenZ)
         return Halo.logarithmic{ vhalo = vhalo, scaleLength = scaleLength, flattenZ = flattenZ }
      end,
      { 73, 100, 50 },
      { 12, 20, 6 },
      { 1, 1.4, 0.7 }
   )

   nfwHalos = buildAllCombinations(
      function(vhalo, scaleLength)
         return Halo.nfw{ vhalo = vhalo, scaleLength = scaleLength }
      end,
      { 155, 200, 120 },
      { 22.25, 30, 14 }
   )

   triaxialHalos = buildAllCombinations(
      function(vhalo, scaleLength, flattenZ, flattenX, flattenY, triaxAngle)
         return Halo.triaxial{ vhalo = vhalo, scaleLength = scaleLength,
                               flattenZ = flattenZ, flattenX = flattenX, flattenY = flattenY,
                               triaxAngle = triaxAngle
                            }
      end,
      { 116, 150, 87 },
      { 16.3, 20, 12 },
      { 1.43, 1.7, 1.2 },
      { 1.26, 1.5, 1.1 },
      { 1.0, 1.2, 0.7 },
      { 96, 120, 66 }
   )

   return mergeTables(logHalos, nfwHalos, traixialHalos)
end

function makeAllPotentials()
   return buildAllCombinations(Potential.create, allSphericals, allDisks, allHalos)
end


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

samplePotentials = { potentialA = makePotentialA(), potentialB = makePotentialB() }

samplePotentialNames = getKeyNames(samplePotentials)

function randomHalo(prng)
   if prng == nil then
      prng = DSFMT.create()
   end

   local type = round(prng:random(0, 2))

   if type == 0 then
      return Halo.logarithmic{
         vhalo       = prng:random(),
         scaleLength = prng:random(),
         flattenZ    = prng:random()
      }
   elseif type == 1 then
      return Halo.nfw{
         vhalo       = prng:random(),
         scaleLength = prng:random()
      }
   elseif type == 2 then
      return Halo.triaxial{
         vhalo       = prng:random(),
         scaleLength = prng:random(),
         flattenX    = prng:random(),
         flattenY    = prng:random(),
         flattenZ    = prng:random(),
         triaxAngle  = prng:random()
      }
   else
      assert(false)
   end
end

function randomDisk(prng)
   if prng == nil then
      prng = DSFMT.create()
   end

   local type = round(prng:random(0, 1))

   if type == 0 then
      return Disk.miyamotoNagai{
         mass        = prng:random(),
         scaleLength = prng:random(),
         scaleHeight = prng:random()
      }
   elseif type == 1 then
      return Disk.exponential{
         mass        = prng:random(),
         scaleLength = prng:random()
      }
   else
      assert(false)
   end
end

function randomSpherical(prng)
   if prng == nil then
      prng = DSFMT.create()
   end

   return Spherical.spherical{
      mass  = prng:random(),
      scale = prng:random()
   }
end

function randomPotential(prng)
   return Potential.create{
      spherical = randomSpherical(prng),
      disk      = randomDisk(prng),
      halo      = randomHalo(prng)
   }
end



