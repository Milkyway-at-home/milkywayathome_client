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

SamplePotentials = { }

SamplePotentials.buildAllDisks = function()
   local noDisks, miyamotoNagaiDisks, doubleExponentialDisks, sech2ExponentialDisks, freemanDisks
   noDisks = buildAllCombinations(
      function(m)
         return Disk.none{ mass = m }
      end,
      { 4.5e5, 5.0e6, 3.0e4 }
   )

   miyamotoNagaiDisks = buildAllCombinations(
      function(m, a, b)
         return Disk.miyamotoNagai{ mass = m, scaleLength = a, scaleHeight = b }
      end,
      { 4.5e5, 5.0e6, 3.0e4 },
      { 6.5, 9.0, 3.0 },
      { 0.26, 0.5, 0.1 }
   )

   doubleExponentialDisks = buildAllCombinations(
      function(m, a, b)
         return Disk.doubleExponential{ mass = m, scaleLength = a, scaleHeight = b }
      end,
      { 2.24933e5, 3.0e6, 1.5e4 },
      { 3.6, 2.0, 6.5 },
      { 0.26, 0.5, 0.1}
   )

   sech2ExponentialDisks = buildAllCombinations(
      function(m, a, b)
         return Disk.sech2Exponential{ mass = m, scaleLength = a, scaleHeight = b }
      end,
      { 2.24933e5, 3.0e6, 1.5e4 },
      { 3.6, 2.0, 6.5 },
      { 0.26, 0.5, 0.1}
   )

   freemanDisks = buildAllCombinations(
      function(m, a)
         return Disk.freeman{ mass = m, scaleLength = a }
      end,
      { 2.24933e5, 3.0e6, 1.5e4 },
      { 3.6, 2.0, 6.5 }
   )

   return mergeTables(noDisks, miyamotoNagaiDisks, doubleExponentialDisks, sech2ExponentialDisks, freemanDisks)
end

SamplePotentials.buildAllDisk2s = function()
   local noDisk2s, miyamotoNagaiDisk2s, doubleExponentialDisk2s, sech2ExponentialDisk2s, freemanDisk2s
   noDisk2s = buildAllCombinations(
      function(m)
         return Disk.none{ mass = m }
      end,
      { 4.5e5, 5.0e6, 3.0e4 }
   )

   miyamotoNagaiDisk2s = buildAllCombinations(
      function(m, a, b)
         return Disk.miyamotoNagai{ mass = m, scaleLength = a, scaleHeight = b }
      end,
      { 4.5e5, 5.0e6, 3.0e4 },
      { 6.5, 9.0, 3.0 },
      { 0.26, 0.5, 0.1 }
   )

   doubleExponentialDisk2s = buildAllCombinations(
      function(m, a, b)
         return Disk.doubleExponential{ mass = m, scaleLength = a, scaleHeight = b }
      end,
      { 2.24933e5, 3.0e6, 1.5e4 },
      { 3.6, 2.0, 6.5 },
      { 0.26, 0.5, 0.1}
   )

   sech2ExponentialDisk2s = buildAllCombinations(
      function(m, a, b)
         return Disk.sech2Exponential{ mass = m, scaleLength = a, scaleHeight = b }
      end,
      { 2.24933e5, 3.0e6, 1.5e4 },
      { 3.6, 2.0, 6.5 },
      { 0.26, 0.5, 0.1}
   )

   freemanDisk2s = buildAllCombinations(
      function(m, a)
         return Disk.freeman{ mass = m, scaleLength = a }
      end,
      { 2.24933e5, 3.0e6, 1.5e4 },
      { 3.6, 2.0, 6.5 }
   )

   return mergeTables(noDisk2s, miyamotoNagaiDisk2s, doubleExponentialDisk2s, sech2ExponentialDisk2s, freemanDisk2s)
end

SamplePotentials.buildAllSphericals = function()
   local noSphericals, herquistSphericals, plummerSphericals
   noSphericals = buildAllCombinations(
      function(m)
         return Spherical.none{ mass = m }
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 }
   )

   herquistSphericals = buildAllCombinations(
      function(m, s)
         return Spherical.herquist{ mass = m, scale = s }
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 },
      { 0.7, 1.1, 0.4 }
   )

   plummerSphericals = buildAllCombinations(
      function(m, s)
         return Spherical.plummer{ mass = m, scale = s }
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 },
      { 0.7, 1.1, 0.4 }
   )

   return mergeTables(noSphericals, herquistSphericals, plummerSphericals)
end


SamplePotentials.buildAllHalos = function()
   local noHalos, logHalos, nfwHalos, triaxialHalos, ASHalos, WEHalos, nfwmHalos, plummerHalos, hernquistHalos, ninkovicHalos

   noHalos = buildAllCombinations(
      function(m)
         return Halo.none{ mass = m }
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 }
   )

   logHalos = buildAllCombinations(
      function(v, L, fZ)
         return Halo.logarithmic{ vhalo = v, scaleLength = L, flattenZ = fZ }
      end,
      { 73, 100, 50 },
      { 12, 20, 6 },
      { 1, 1.4, 0.71 }
   )

   nfwHalos = buildAllCombinations(
      function(v, L)
         return Halo.nfw{ vhalo = v, scaleLength = L }
      end,
      { 155, 200, 120 },
      { 22.25, 30, 14 }
   )

   triaxialHalos = buildAllCombinations(
      function(v, L, fZ, fX, fY, angle)
         return Halo.triaxial{ vhalo = v, scaleLength = L,
                               flattenZ = fZ, flattenX = fX, flattenY = fY,
                               triaxAngle = angle
                            }
      end,
      { 116, 150, 87 },
      { 16.3, 20, 12 },
      { 1.43, 1.7, 1.2 },
      { 1.26, 1.5, 1.1 },
      { 1.0, 1.2, 0.7 },
      { 96, 120, 66 }
   )

   ASHalos = buildAllCombinations(
      function(m, L, l, g)
         return Halo.allenSantillan{ mass = m, scaleLength = L, lambda = l, gamma = g}
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 },
      { 16.3, 20, 12 },
      { 100.0, 200.0, 50.0},
      { 2.0, 3.0, 4.0 }
   )

   WEHalos = buildAllCombinations(
      function(m, L)
         return Halo.wilkinsonEvans{ mass = m, scaleLength = L }
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 },
      { 16.3, 20, 12 }
   )

   nfwmHalos = buildAllCombinations(
      function(m, L)
         return Halo.nfwmass{ mass = m, scaleLength = L }
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 },
      { 16.3, 20, 12 }
   )

   plummerHalos = buildAllCombinations(
      function(m, L)
         return Halo.plummer{ mass = m, scaleLength = L }
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 },
      { 16.3, 20, 12 }
   )

   hernquistHalos = buildAllCombinations(
      function(m, L)
         return Halo.hernquist{ mass = m, scaleLength = L }
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 },
      { 16.3, 20, 12 }
   )

   ninkovicHalos = buildAllCombinations(
      function(p, L, l)
         return Halo.ninkovic{ rho0 = p, scaleLength = L, lambda = l }
      end,
      { 1.52954402e5, 1.7e6, 1.1e4 },
      { 16.3, 20, 12 },
      { 70.0, 96.5, 150.0}
   )

   return mergeTables(noHalos, logHalos, nfwHalos, triaxialHalos, ASHalos, WEHalos, nfwmHalos, plummerHalos, hernquistHalos, ninkovicHalos)
end

SamplePotentials.makeAllPotentials = function()
  return buildAllCombinations(Potential.create,
                              SamplePotentials.buildAllSphericals(),
                              SamplePotentials.buildAllDisks(),
                              SamplePotentials.buildAllDisk2s(),
                              SamplePotentials.buildAllHalos())

end

-- Artificial tests, should have at least every combination of
-- spherical, disk, halo component typs. (750 COMBINATIONS!?!)

function Spherical_piece(i)
   if (i==0) then
      return Spherical.none{ mass = 3.0e5 }
   elseif (i==1) then
      return Spherical.hernquist{ mass = 1.5e5, scale = 0.8 }
   elseif (i==2) then
      return Spherical.plummer{ mass = 2.0e5, scale = 0.6 }
   else
      assert(false)
   end
end

function Disk_piece(j)
   if (j==0) then
      return Disk.none{ mass = 3.0e5 }
   elseif (j==1) then
      return Disk.miyamotoNagai{ mass = 4.5e5, scaleLength = 6.0, scaleHeight = 0.3 }
   elseif (j==2) then
      return Disk.doubleExponential{ mass = 4.0e5, scaleLength = 4.5, scaleHeight = 0.3 }
   elseif (j==3) then
      return Disk.sech2Exponential{ mass = 4.5e5, scaleLength = 6.0, scaleHeight = 0.3 }
   elseif (j==4) then
      return Disk.freeman{ mass = 4.5e5, scaleLength = 6.0 }
   else
      assert(false)
   end
end

function Disk2_piece(k)
   if (k==0) then
      return Disk.none{ mass = 3.0e5 }
   elseif (k==1) then
      return Disk.miyamotoNagai{ mass = 3.0e5, scaleLength = 3.0, scaleHeight = 0.3 }
   elseif (k==2) then
      return Disk.doubleExponential{ mass = 3.0e5, scaleLength = 1.5, scaleHeight = 0.3 }
   elseif (k==3) then
      return Disk.sech2Exponential{ mass = 3.0e5, scaleLength = 3.0, scaleHeight = 0.3 }
   elseif (k==4) then
      return Disk.freeman{ mass = 3.0e5, scaleLength = 2.5 }
   else
      assert(false)
   end
end

function Halo_piece(l)
   if (l==0) then
      return Halo.none{ mass = 3.0e6 }
   elseif (l==1) then
      return Halo.logarithmic{ vhalo = 80, scaleLength = 15, flattenZ = 1.0 }
   elseif (l==2) then
      return Halo.nfw{ vhalo = 90, scaleLength = 12 }
   elseif (l==3) then
      return Halo.triaxial{ vhalo = 120, scaleLength = 18, flattenX = 1.3, flattenY = 1.0, flattenZ = 1.45, triaxAngle = 96 }
   elseif (l==4) then
      return Halo.allenSantillan{ mass = 3.0e6, scaleLength = 20, lambda = 200, gamma = 2.0}
   elseif (l==5) then
      return Halo.wilkinsonEvans{ mass = 3.0e6, scaleLength = 15 }
   elseif (l==6) then
      return Halo.nfwmass{ mass = 3.0e6, scaleLength = 15 }
   elseif (l==7) then
      return Halo.plummer{ mass = 3.0e6, scaleLength = 15 }
   elseif (l==8) then
      return Halo.hernquist{ mass = 3.0e6, scaleLength = 15 }
   elseif (l==9) then
      return Halo.ninkovic{ rho0 = 26.0, scaleLength = 15, lambda = 96 }
   else
      assert(false)
   end
end

--Following for-loop generates all possible potential combinations
sphere_types = 3
disk_types = 5
halo_types = 10
all_combo = false        --Keep this set to false if you DON'T want to run over all possible potential combinations

SamplePotentials.samplePotentials = { }

if (all_combo) then
   pot_count = 1
   for n=0,sphere_types*disk_types*disk_types*halo_types+1 do
      i = floor(n/(disk_types*disk_types*halo_types))
      j = floor((n-(disk_types*disk_types*halo_types)*i)/(disk_types*halo_types))
      k = floor((n-(disk_types*disk_types*halo_types)*i-(disk_types*halo_types)*j)/(halo_types))
      l = n-(disk_types*disk_types*halo_types)*i-(disk_types*halo_types)*j-halo_types*k

      SamplePotentials.samplePotentials["potential"..tostring(pot_count)] = Potential.create{
         spherical = Spherical_piece(i),
         disk      = Disk_piece(j),
         disk2     = Disk2_piece(k),
         halo      = Halo_piece(l)
      }

      pot_count = pot_count + 1
      if (pot_count > (sphere_types*disk_types*disk_types*halo_types)) then
         break
      end
   end
else

   disk1_count = 1
   for n=0, disk_types+1 do
      SamplePotentials.samplePotentials["potential_prim_disk"..tostring(disk1_count)] = Potential.create{
         spherical = Spherical_piece(0),
         disk      = Disk_piece(n),
         disk2     = Disk2_piece(0),
         halo      = Halo_piece(0)
      }

      disk1_count = disk1_count + 1
      if (disk1_count > disk_types) then
         break
      end
   end

   disk2_count = 1
   for n=0, disk_types+1 do
      SamplePotentials.samplePotentials["potential_sec_disk"..tostring(disk2_count)] = Potential.create{
         spherical = Spherical_piece(0),
         disk      = Disk_piece(0),
         disk2     = Disk2_piece(n),
         halo      = Halo_piece(0)
      }

      disk2_count = disk2_count + 1
      if (disk2_count > disk_types) then
         break
      end
   end

   sphere_count = 1
   for n=0, sphere_types+1 do
      SamplePotentials.samplePotentials["potential_sphere"..tostring(sphere_count)] = Potential.create{
         spherical = Spherical_piece(n),
         disk      = Disk_piece(0),
         disk2     = Disk2_piece(0),
         halo      = Halo_piece(0)
      }

      sphere_count = sphere_count + 1
      if (sphere_count > sphere_types) then
         break
      end
   end

   halo_count = 1
   for n=0, halo_types+1 do
      SamplePotentials.samplePotentials["potential_halo"..tostring(halo_count)] = Potential.create{
         spherical = Spherical_piece(0),
         disk      = Disk_piece(0),
         disk2     = Disk2_piece(0),
         halo      = Halo_piece(n)
      }

      halo_count = halo_count + 1
      if (halo_count > halo_types) then
         break
      end
   end

end

SamplePotentials.samplePotentialNames = getKeyNames(SamplePotentials.samplePotentials)

SamplePotentials.randomHalo = function(prng)
   if (prng == nil) then
      prng = DSFMT.create()
   end

   local typ = floor(prng:random(0, 10))

   if (typ == 0) then
      return Halo.none{
         mass        = prng:random(1.0e4, 5.0e7)
      }
   elseif (typ == 1) then
      return Halo.logarithmic{
         vhalo       = prng:random(10, 200),
         scaleLength = prng:random(3, 100),
         flattenZ    = prng:random(0.708, 3)
      }
   elseif (typ == 2) then
      return Halo.nfw{
         vhalo       = prng:random(1, 300),
         scaleLength = prng:random(0.1, 40)
      }
   elseif (typ == 3) then
      return Halo.triaxial{
         vhalo       = prng:random(1, 200),
         scaleLength = prng:random(0.1, 30),
         flattenX    = prng:random(0, 4),
         flattenY    = prng:random(0, 3),
         flattenZ    = prng:random(0, 4),
         triaxAngle  = prng:random(0, 180)
      }
   elseif (typ == 4) then
      return Halo.allenSantillan{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(0.1, 30),
         lambda      = prng:random(20.0, 300.0),
         gamma       = prng:random(2.0, 5.0)
      }
   elseif (typ == 5) then
      return Halo.wilkinsonEvans{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(0.1, 30)
      }
   elseif (typ == 6) then
      return Halo.nfwmass{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(0.1, 30)
      }
   elseif (typ == 7) then
      return Halo.plummer{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(0.1, 30)
      }
   elseif (typ == 8) then
      return Halo.hernquist{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(0.1, 30)
      }
   elseif (typ == 9) then
      return Halo.ninkovic{
         rho0        = prng:random(5.0, 50.0),
         scaleLength = prng:random(0.1, 30),
         lambda      = prng:random(30.0, 150.0)
      }
   else
      assert(false)
   end
end

SamplePotentials.randomDisk = function(prng)
   if (prng == nil) then
      prng = DSFMT.create()
   end

   local typ = floor(prng:random(0, 5))

   if (typ == 0) then
      return Disk.none{
         mass        = prng:random(1.0e4, 5.0e7)
      }
   elseif (typ == 1) then
      return Disk.miyamotoNagai{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(1, 10),
         scaleHeight = prng:random(0.1, 8)
      }
   elseif (typ == 2) then
      return Disk.doubleExponential{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(1, 10),
         scaleHeight = prng:random(0.1, 8)
      }
   elseif (typ == 3) then
      return Disk.sech2Exponential{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(1, 10),
         scaleHeight = prng:random(0.1, 8)
      }
   elseif (typ == 4) then
      return Disk.freeman{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(1, 10)
      }
   else
      assert(false)
   end
end

SamplePotentials.randomDisk2 = function(prng)
   if (prng == nil) then
      prng = DSFMT.create()
   end

   local typ = floor(prng:random(0, 5))

   if (typ == 0) then
      return Disk.none{
         mass        = prng:random(1.0e4, 5.0e7)
      }
   elseif (typ == 1) then
      return Disk.miyamotoNagai{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(1, 10),
         scaleHeight = prng:random(0.1, 8)
      }
   elseif (typ == 2) then
      return Disk.doubleExponential{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(1, 10),
         scaleHeight = prng:random(0.1, 8)
      }
   elseif (typ == 3) then
      return Disk.sech2Exponential{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(1, 10),
         scaleHeight = prng:random(0.1, 8)
      }
   elseif (typ == 4) then
      return Disk.freeman{
         mass        = prng:random(1.0e4, 5.0e7),
         scaleLength = prng:random(1, 10)
      }
   else
      assert(false)
   end
end

SamplePotentials.randomSpherical = function(prng)
   if (prng == nil) then
      prng = DSFMT.create()
   end

   local typ = floor(prng:random(0, 3))

   if (typ == 0) then
      return Spherical.none{
         mass        = prng:random(1.0e4, 5.0e7)
      }
   elseif (typ == 1) then
      return Spherical.hernquist{
         mass  = prng:random(1.0e4, 9e7),
         scale = prng:random(0, 4)
      }
   elseif (typ == 2) then
      return Spherical.plummer{
         mass  = prng:random(1.0e4, 9e7),
         scale = prng:random(0, 4)
      }
   else
      assert(false)
   end
end

SamplePotentials.randomPotential = function(prng)
   return Potential.create{
      spherical = SamplePotentials.randomSpherical(prng),
      disk      = SamplePotentials.randomDisk(prng),
      disk2     = SamplePotentials.randomDisk2(prng),
      halo      = SamplePotentials.randomHalo(prng)
   }
end


return SamplePotentials


