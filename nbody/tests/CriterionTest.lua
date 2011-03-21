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

require "nbody_testing"
require "sample_models"
require "sample_potentials"


local function createTestCritCtx(crit)
   local ctx = NBodyCtx.createTestCtx()
   ctx.criterion = crit
   ctx.theta = 1.0
   return ctx
end

local newCritCtx = createTestCritCtx("NewCriterion")
local sw93Ctx = createTestCritCtx("SW93")
local bh86Ctx = createTestCritCtx("BH86")
local exactCtx = createTestCritCtx("Exact")

print("whee = ", findRCrit(newCritCtx, Vector.create(3, 4, 5), 5, Vector.create(5, 6, 2), 2.5))






