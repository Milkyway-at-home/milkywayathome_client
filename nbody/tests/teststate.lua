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


allCriterion = { "sw93", "NewCriterion", "BH86", "Exact" }
useQuads = { true, false }
thetas = { 1.0, 0.9, 0.7, 0.5, 0.3 }
allowIncests = useQuads
treeRsizes = { 4.0, 8.0, 2.0, 1.0 }

function createNBodyCtxRaw(theta, rsize, crit, useQuad, allowIncest)
   assert(theta ~= nil, "theta")
   assert(rsize ~= nil, "rsize")
   assert(crit ~= nil, "crit")
   assert(useQuad ~= nil, "quad")
   assert(allowIncest ~= nil, "incest")
   local rawCtx = { }
   rawCtx.theta = theta
   rawCtx.treeRSize = rsize
   rawCtx.criterion = crit
   rawCtx.useQuad = useQuad
   rawCtx.allowIncest = allowIncest
   return rawCtx
end

function whee(a)
   print("--------------------")
   table.foreach(a, print)
   print("--------------------")
   return nil
end

everything = buildAllCombinations(createNBodyCtxRaw,
                                  thetas,
                                  treeRsizes,
                                  allCriterion,
                                  useQuads,
                                  allowIncests)

 print("WHEEEEE", #everything)
 map(whee, everything)

function testF(a, b)
   local x = { }
   x.omg = b
   x.wow = a
   return x
end

function testF3(a, b ,c)
   local x = { }
   x.bag = a
   x.wow = c
   x.omg = b

   return x
end

testE = buildAllCombinations(testF3, { "a", "b" }, { "c", "d" }, { "e", "f", "g" })
-- a c, a d, b c, b d

print("arstarstarst", #testE, testE)
map(whee, testE)

print("Derp1", rawCtxHash(everything[1]))
print("Derp2", rawCtxHash(everything[1]))

function getKeyAnswerPair(rawCtx)
   local hash = rawCtxHash(rawCtx)
   rawCtx.result = "I am a sample result"
   return { hash, rawCtx }
end


