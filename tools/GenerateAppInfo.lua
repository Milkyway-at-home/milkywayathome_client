
argv = { ... }

separation_bin = argv[1]
nbody_bin = argv[2]
nbody_graphics_bin = argv[3]

print(separation_bin, nbody_bin, nbody_graphics_bin)
print("arsotenareotneiarstnarntarnst")

-- <ncpu>4</ncpu>
-- <avg_ncpus>4</avg_ncpus>
-- <max_ncpus>4</max_ncpus>

if separation_bin == "NOTFOUND" then
   separation_bin = nil
end

if nbody_bin == "NOTFOUND" then
   nbody_bin = nil
end

if nbody_graphics_bin == "NOTFOUND" then
   nbody_graphics_bin = nil
end

if separation_bin == nil and nbody_bin == nil then
   io.stderr:write("Nothing to install\n")
   os.exit(1)
end

base_template =
[[    <app>
    <name>%s</name>
    <user_friendly_name>%s</user_friendly_name>
      </app>]]

executable_template =
[[    <file_info>
    <name>%s</name>
    <executable/>
    </file_info>]]


-- Order of magnitude estimate
function flopsEstimateFromPlanClass(pc)
   if pc == "cuda_opencl" then
      return 1.0e11
   elseif pc == "ati14" then
      return 1.0e11
   elseif pc == "mt" then
      return 100.0e9
   else
      return 100.0e9
   end
end

function getFlops(flops)
   return string.format("    <flops>%f</flops>", flops)
end

beginAppVersion = "    <app_version>"
endAppVersion = "    </app_version>"

function getExeFileInfo(name)
   if name == nil then
      return nil
   else
      return string.format(executable_template, name)
   end
end

function getAppName(name)
   return string.format("    <app_name>%s</app_name>", name)
end

function getPlanClass(class)
   if class == nil then
      return nil
   else
      return string.format("    <plan_class>%s</plan_class>", class)
   end
end

function findVersionFromName(name)
   return name:match("_(%d.%d+)_")
end

function findPlanClassFromName(name)
   assert(name, "Name not found")
   return name:match("__(.+)")
end

function getVersionNum(version)
   if version == nil then
      return nil
   else
      return string.format("    <version_num>%d</version_num>", 100 * version)
   end
end

fileRefMainProgramTemplate = [[
    <file_ref>
      <file_name>%s</file_name>
      <main_program/>
    </file_ref>]]

fileRefGraphicsTemplate = [[
    <file_ref>
      <file_name>%s</file_name>
      <open_name>graphics_app</open_name>
    </file_ref>]]

function getFileRefMainProgram(file)
   assert(file, "Main file must exist")
   return string.format(fileRefMainProgramTemplate, file)
end

function getFileRefGraphics(file)
   if file == nil then
      return nil
   else
      return string.format(fileRefGraphicsTemplate, file)
   end
end

function buildAppInfo(name, ufName, binName, graphicsName)

   local pc = getPlanClass(findPlanClassFromName(binName))


   local parts = {
      string.format(base_template,
                    name,
                    ufName),
      getExeFileInfo(binName),
      getExeFileInfo(graphicsName),
      beginAppVersion,
        getAppName(name),
        getVersionNum(findVersionFromName(binName)),
        pc,
        getFlops(flopsEstimateFromPlanClass(pc)),
        getFileRefMainProgram(binName),
        getFileRefGraphics(graphicsName),
      endAppVersion
     }


   local realParts = { }
   for k, v in pairs(parts) do
      if v ~= nil then
         realParts[#realParts + 1] = v
      end
   end

   return table.concat(realParts, "\n")
end

function buildNbodyXML()
   return buildAppInfo("milkyway_nbody",
                       "MilkyWay@Home N-body",
                       nbody_bin,
                       nbody_graphics_bin
                    )
end

function buildSeparationXML()
   return buildAppInfo("milkyway_separation",
                       "MilkyWay@Home Separation",
                       separation_bin,
                       nil
                    )
end



 header = [[
 <?xml version="1.0"?>
 <app_info>]]

 endHeader = [[
 </app_info>
 ]]

if nbody_bin ~= nil then
   nbody_xml = buildNbodyXML()
   assert(nbody_xml)
end

if separation_bin ~= nil then
   separation_xml = buildSeparationXML()
   assert(separation_xml)
end


app_info_parts = { header, nbody_xml, "\n", separation_xml, endHeader }
app_info_xml = { }
for k, v in pairs(app_info_parts) do
   if v ~= nil then
      app_info_xml[#app_info_xml + 1] = v
   end
end


app_info = io.open("app_info.xml", "w")
assert(app_info, "Failed to open app_info.xml for writing")
app_info:write(table.concat(app_info_xml, "\n"))
app_info:close()


