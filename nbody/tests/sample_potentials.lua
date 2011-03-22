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

   return mergeTables(logHalos, nfwHalos, triaxialHalos)
end

function makeAllPotentials()
   return buildAllCombinations(Potential.create, allSphericals, allDisks, allHalos)
end

-- Artificial tests, should have at least every combination of
-- spherical, disk, halo component types.
samplePotentials = {
   potentialA = Potential.create{
      spherical = Spherical.spherical{
         mass  = 1.5e5,
         scale = 0.7
      },

      disk = Disk.miyamotoNagai{
         mass        = 4.5e5,
         scaleLength = 6.0,
         scaleHeight = 0.3
      },

      halo = Halo.logarithmic{
         vhalo       = 80,
         scaleLength = 15,
         flattenZ    = 1.0
      }
   },

   potentialB = Potential.create{
      spherical = Spherical.spherical{
         mass  = 2.0e5,
         scale = 0.7
      },

      disk = Disk.exponential{
         mass        = 9.0e5,
         scaleLength = 6
      },

      halo = Halo.logarithmic{
         vhalo       = 60,
         scaleLength = 10.0,
         flattenZ    = 1.1
      }
   },

   potentialC = Potential.create{
      spherical = Spherical.spherical{
         mass  = 2.0e6,
         scale = 1.0
      },

      disk = Disk.exponential{
         scaleLength = 9,
         mass        = 9.0e5
      },

      halo = Halo.nfw{
         vhalo       = 150,
         scaleLength = 25
      }
   },

   potentialD = Potential.create{
      spherical = Spherical.spherical{
         mass  = 2.0e6,
         scale = 1.0
      },

      disk = Disk.miyamotoNagai{
         mass        = 2e6,
         scaleLength = 7.0,
         scaleHeight = 0.2
      },

      halo = Halo.triaxial{
         vhalo       = 120,
         scaleLength = 18,
         flattenZ = 1.45,
         flattenX = 1.3,
         flattenY = 1.0,
         triaxAngle = 96
      }
   },

   potentialE = Potential.create{
      spherical = Spherical.spherical{
         mass  = 2e5,
         scale = 0.4
      },

      disk = Disk.miyamotoNagai{
         mass        = 7e5,
         scaleLength = 9.0,
         scaleHeight = 0.5
      },

      halo = Halo.nfw{
         vhalo       = 120,
         scaleLength = 18,
      }
   },

   potentialF = Potential.create{
      spherical = Spherical.spherical{
         mass  = 2.0e6,
         scale = 1.0
      },

      disk = Disk.exponential{
         mass        = 3.0e7,
         scaleLength = 8,
      },

      halo = Halo.triaxial{
         vhalo       = 140,
         scaleLength = 25,
         flattenZ = 1.3,
         flattenX = 1.4,
         flattenY = 1.1,
         triaxAngle = 112
      }
   }
}

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



