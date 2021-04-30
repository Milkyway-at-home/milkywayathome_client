
require "NBodyTesting"

args = {...}

nbodyBin = assert(args[1], "Missing binary name")
inputTest = "ArgumentTestInput.lua"

local simple = os.readProcess(nbodyBin,
                              "--debug-boinc",
                              "--verify-file",
                              "--input-file", inputTest,
                              "-np 5 -p hello -4 -3.14 1.23 13434")
if simple:find("File is OK") == nil then
   eprintf("simple failed: '%s'\n", simple)
   os.exit(1)
else
   printf("simple passed\n")
end



local startNthread = os.readProcess(nbodyBin,
                                    "--debug-boinc",
                                    "--nthreads 5",
                                    "--verify-file",
                                    "--input-file", inputTest,
                                    "-np 5 -p hello -4 -3.14 1.23 13434")
if startNthread:find("File is OK") == nil then
   eprintf("start nthread failed:\n")
   error(startNthread)
else
   printf("start nthread passed\n")
end


local endNthread = os.readProcess(nbodyBin,
                                  "--debug-boinc",
                                  "--verify-file",
                                  "--input-file", inputTest,
                                  "-np 5 -p hello -4 -3.14 1.23 13434",
                                  "--nthreads 5")
if endNthread:find("File is OK") == nil then
   eprintf("end nthread failed:\n")
   error(endNthread)
else
   printf("end nthread passed\n")
end


local multiEnd = os.readProcess(nbodyBin,
                                "--debug-boinc",
                                "--verify-file",
                                "--input-file", inputTest,
                                "-np 5 -p hello -4 -3.14 1.23 13434",
                                "--nthreads 5 --progress")
if multiEnd:find("File is OK") == nil then
   eprintf("multiend failed:\n")
   error(multiEnd)
else
   printf("multiend passed\n")
end




