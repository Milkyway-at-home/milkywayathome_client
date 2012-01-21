
require "NBodyTesting"
require "persistence"

args = { ... }

if #args ~= 3 then
   eprintf("usage: output input1 input2)")
   os.exit(1)
end



-- merge results tables into t1
-- avoids clobbering the "samples" and recalculates the statisics.
function mergeResultTables(t1, t2)
   for k,v in pairs(t2) do
      if type(v) == "table" then
         if type(t1[k] or false) == "table" then
            if k == "samples" then
               local t1Samples, t2Samples = t1["samples"], t2["samples"]
               local i = #t1Samples
               for kk, vv in ipairs(t2Samples) do
                  i = i + 1
                  t1Samples[i] = vv
               end
            else
               mergeResultTables(t1[k] or {}, t2[k] or {})
            end
         else
            t1[k] = v
         end
      else
         t1[k] = v
      end

      if k == "mean" then
         t1["mean"] = calcMean(t1["samples"])
      elseif k == "stddev" then
         t1["stddev"] = calcStddev(t1["samples"])
      elseif k == "min" then
         local min, _ = findMinMax(t1["samples"])
         t1["min"] = min
      elseif k == "max" then
         local _, max = findMinMax(t1["samples"])
         t1["max"] = max
      end
   end
   return t1
end

local t1 = persistence.load(args[2])
local t2 = persistence.load(args[3])

--printTable(t1)

local merged = mergeResultTables(t1, t2)
persistence.store(args[1], merged)




