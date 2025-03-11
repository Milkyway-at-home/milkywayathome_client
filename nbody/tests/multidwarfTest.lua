-- 
-- >> MultidwarfTest <<
-- checks that dwarfs generate correctly with the 「multidwarfs」 update under various generation parameters
-- many settings are kept as the base values from for_developers.lua
--



--require "NBodyTesting"

SP = require "SamplePotentials"

dec = 9.0   -- digits to round to for calculations

print("0")

-- stolen from for_developers.lua
function round(num, places)
    local mult = 10.0^(places)
    return floor(num * mult + 0.5) / mult
  end
  
print("1")

if prng == nil then
    prng = DSFMT.create()
end

print("2")

local ntests = 20
modelComponents = 2     --
SunGCDist = 8.0
sigma = 2.5
sigmaIter = 6
correct = 1.111

ndwarfs = 4             --
totalBodies = 1000      --
dt = 0.001      --

LMC_presence = false
LMC_mass     = 449865.888
LMC_rscale   = 15.0
LMC_dynafric = true

print("3")

-- maybe there's a way to initialize all these arrays more efficiently? idrk
-- currently uses real values for SMC/Sgr./Fornax/Sculptor Dwarfs
orbit_b = {round(-44.33, dec), round(-14.17, dec), round(-65.65, dec), round(-83.16, dec)}
orbit_r = {round(62.4, dec), round(25, dec), round(143, dec), round(88.91, dec)}
orbit_vx = {round(21.99, dec), round(223.97, dec), round(-27.04, dec), round(-22.11, dec)}
orbit_vy = {round(-201.36, dec), round(-5.34, dec), round(-172.14, dec), round(197.28, dec)}
orbit_vz = {round(171.25, dec), round(185.78, dec), round(101.21, dec), round(-102.1, dec)}

orbit_l = {round(302.8, dec), round(5.57, dec), round(237.1, dec), round(287.54, dec)}
dwarf_mass_l = {round(2429.198, dec), round(107.041, dec), round(80.159, dec), round(9.384, dec)}
dwarf_lm_ratio = {round(0.08307, dec), round(0.05949, dec), round(0.14255, dec), round(0.06729, dec)}
dwarf_rscale_l = {round(2.9, dec), round(1.53, dec), round(1.425, dec), round(0.725, dec)}
dwarf_lr_ratio = {round(0.2, dec), round(0.2, dec), round(0.2, dec), round(0.2, dec)}

dwarf_mass_d = {round(26813.5851948958, dec), round(1692.2698085392, dec), round(482.1629922834, dec), round(130.0720855996, dec)}
dwarf_rscale_d = {round(14.5, dec), round(7.65, dec), round(7.125, dec), round(3.625, dec)}
-- for i = 1, ndwarfs do
--     orbit_b[i] = round(prng:random(-90.0,90.0), dec)
--     orbit_r[i] = round(prng:random(10.0,60.0), dec)
--     orbit_vx[i] = round(prng:random(-200.0,200.0), dec)
--     orbit_vy[i] = round(prng:random(-200.0,200.0), dec)
--     orbit_vz[i] = round(prng:random(-200.0,200.0), dec)

--     orbit_l[i] = round(prng:random(0.0, 360.0), dec)
--     dwarf_mass_l[i] = round((prng:random(0.95, 5.0))^3, dec)    -- actual range 0.857-125
--     dwarf_lm_ratio[i] = round((prng:random(0.075, 0.55))^3, dec)  -- actual range 0.00042-0.166
--     dwarf_rscale_l[i] = round(prng:random(0.1, 2.0), dec)
--     dwarf_lr_ratio[i] = round(0.2, dec)

--     dwarf_mass_d[i] =  (dwarf_mass_l[i]/dwarf_lm_ratio[i]) - dwarf_mass_l[i]
--     dwarf_rscale_d[i] = (dwarf_rscale_l[i]/dwarf_lr_ratio[i])
-- end

print("4")

-- ripped out of for_developers.lua
function soft()
    if (ModelComponents == 1) then --plugs in two-comp. analog for single-comp. run so i don't have to edit the eps2 function
        sp = calculateEps2(totalBodies, dwarf_scale_l[1], dwarf_rscale_d[1], dwarf_mass_l[1]/2, dwarf_mass_d[1]/2, 0)
    else
        sp = calculateEps2(totalBodies, dwarf_rscale_l[1], dwarf_rscale_d[1], dwarf_mass_l[1], dwarf_mass_d[1], 0)
    end
    return sp
end

print("5")

ctx = NBodyCtx.create{
    dwarfn        = ndwarfs,
    timestep      = round(dt, 3),
    timeEvolve    = round(dt, 3),
    timeBack      = round(dt, 3),
    theta         = 1.0,
    eps2          = soft(),
    b             = orbit_b,
    r             = orbit_r,
    vx            = orbit_vx,
    vy            = orbit_vy,
    vz            = orbit_vz,
    sunGCDist     = SunGCDist,
    criterion     = "TreeCode",
    useQuad       = true,
    useBestLike   = false,
    BestLikeStart = prng:random(0.85,0.99),
    BetaSigma     = sigma,
    VelSigma      = sigma,
    DistSigma     = sigma,
    BetaCorrect   = correct,
    VelCorrect    = correct,
    DistCorrect   = correct,
    IterMax       = sigmaIter,
    allowIncest   = true,
    quietErrors   = true,
    LMC           = LMC_presence,
    LMCmass       = LMC_mass,
    LMCscale      = LMC_rscale,
    LMCDynaFric   = LMC_dynafric
}

print("6")

potential = Potential.create{
    spherical = Spherical.hernquist{ mass  = 1.52954402e5, scale = 0.7 },
    disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
    disk2     = Disk.none{ mass = 3.0e5 },
    halo      = Halo.logarithmic{ vhalo = 74.61, scaleLength = 12.0, flattenZ = 1.0 }
}

print("7")

function makeDwarfs(ctx, potential)
    local firstModel = {}
    local finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity = {}, {}

    if (LMC_body) then
        local potential = potential
        local position = lbrToCartesianTable(ctx, Vector.creates(orbit_l, orbit_b, orbit_r))
        local velocity = Vector.creates(orbit_vx, orbit_vy, orbit_vz)
        local LMCposition = Vector.create(-1.1, -41.1, -27.9)
        local LMCvelocity = Vector.create(-57, -226, 221)
        local LMCmass = LMC_Mass
        local LMCscale = LMC_scaleRadius
        local LMCDynaFric = LMC_DynamicalFriction and 1 or 0
        local coulomb_log = CoulombLogarithm
        local ftime = evolveTime
        local tstop = revOrbTime
        local dt = ctx.timestep / 10.0      
        local masses    = dwarfMass  
        local rscales = rscale_t

        finalPosition, finalVelocity, LMCfinalPosition, LMCfinalVelocity = reverseOrbitS_LMC(potential, position, velocity, LMCposition, LMCvelocity, LMCmass, LMCscale, LMCDynaFric, coulomb_log, ftime, tstop, dt, masses, rscales)      
    else
        local potential = potential
        local position  = lbrToCartesianTable(ctx, Vector.creates(orbit_l, orbit_b, orbit_r))
        local velocity  = Vector.creates(orbit_vx, orbit_vy, orbit_vz)
        local tstop     = revOrbTime
        local dt        = ctx.timestep / 10.0 
        local masses    = dwarfMass

        finalPosition, finalVelocity = reverseOrbitS(potential, position, velocity, tstop, dt, masses)
    end
    
    if(ModelComponents == 2) then 
        for i = 1, ndwarfs do
            local Model = predefinedModels.mixeddwarf{
                nbody       = totalBodies,
                prng        = prng,
                position    = finalPosition[i],
                velocity    = finalVelocity[i],
                comp1       = Dwarf.plummer{mass = dwarf_mass_l[i], scaleLength = dwarf_rscale_l[i]}, -- Dwarf Options: plummer, nfw, general_hernquist
                comp2       = Dwarf.plummer{mass = dwarf_mass_d[i], scaleLength = dwarf_rscale_d[i]}, -- Dwarf Options: plummer, nfw, general_hernquist
                ignore      = true
                }
            for _, row in ipairs(Model) do
                table.insert(firstModel, row)
            end
            print(string.format("Dwarf %d bodies generation finished", i))
        end
        
    elseif(ModelComponents == 1) then
        for i = 1, n do
            local Model = predefinedModels.plummer{
                nbody       = totalBodies,
                prng        = prng,
                position    = finalPosition[i],
                velocity    = finalVelocity[i],
                mass        = mass_l[i],
                scaleRadius = rscale_l[i],
                ignore      = true
                }
            for _, row in ipairs(Model) do
                table.insert(firstModel, row)
                print(row)
            end
            print(string.format("Dwarf %d bodies generation finished", i))
        end
    end
    return firstModel
end

print("8")

st = makeDwarfs(ctx, potential)
st:step(ctx)

print("Bodies generated successfully.")