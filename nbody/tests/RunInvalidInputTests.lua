
require "NBodyTesting"

local args = { ... }

assert(args and #args > 1, "There must be at least 2 arguments")

local verbose = true

local nbodyBin = args[1]

local nbodyArgs = {
   nbodyBin = nbodyBin,
   extraArgs = { "--verify-file" }
}

local failCount = 0

-- try running file verify on each given file argument
for i = 2, #args do
   print("Trying ", args[i])
   nbodyArgs.input = args[i]

   local output, exitCode = runSimple(nbodyArgs)

   if exitCode == 0 then
      eprintf("Input '%s' should have failed but did not\n", args[i])
      failCount = failCount + 1
   else
      printf("Inpt '%s' failed as expected\n", args[i])
   end

   if verbose then
      io.stdout:write(output)
   end
end

eprintf("%d/%d tests failed.\n", failCount, #args - 1)

if failCount ~= 0 then
   os.exit(1)
end

