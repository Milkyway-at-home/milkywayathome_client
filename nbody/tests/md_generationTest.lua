-- 
-- >> MultidwarfTest <<
-- checks that dwarfs generate correctly with the 「multidwarfs」 update under various generation parameters
-- many settings are kept as the base values from for_developers.lua
--



require "NBodyTesting"

dec = 9.0   -- digits to round to for calculations
dt = 0.001
-- stolen from for_developers.lua
function round(num, places)
    local mult = 10.0^(places)
    return floor(num * mult + 0.5) / mult
  end

function cross(x1, y1, z1, x2, y2, z2)
    local x3 = y1*z2-y2*z1
    local y3 = x2*z1-x1*z2
    local z3 = x1*y2-x2*y1
    return {x=x3, y=y3, z=z3}
end

function abs(x)
    if(x >= 0) then
        return x
    else
        return -x
    end
end
argSeed = 34086709
prng = DSFMT.create(argSeed)

modelComponents = 1     --
SunGCDist = 8.0
sigma = 2.5
sigmaIter = 6
correct = 1.111

ndwarfs = 4             --
nbodies = 4000000      --
dt = 0.001      --

LMC_presence = false
LMC_mass     = 449865.888
LMC_rscale   = 15.0
LMC_dynafric = true

----------------------
-- DWARF PARAMETERS --
----------------------

-- maybe there's a way to initialize all these arrays more efficiently? idrk
-- currently uses real values for SMC/Sagittarius/Fornax/Sculptor Dwarfs
orbit_b = {round(-44.33, dec), round(-14.17, dec), round(-65.65, dec), round(-83.16, dec)}
orbit_r = {round(62.4, dec), round(25, dec), round(143, dec), round(88.91, dec)}
orbit_vx = {round(21.99, dec), round(223.97, dec), round(-27.04, dec), round(-22.11, dec)}
orbit_vy = {round(-201.36, dec), round(-5.34, dec), round(-172.14, dec), round(197.28, dec)}
orbit_vz = {round(171.25, dec), round(185.78, dec), round(101.21, dec), round(-102.1, dec)}

orbit_l = {round(302.8, dec), round(5.57, dec), round(237.1, dec), round(287.54, dec)}
-- vv for 1-comp mode vv --
dwarf_mass_l = {round(29241.283, dec), round(1799.464, dec), round(562.332, dec), round(139.458, dec)}
dwarf_rscale_l = {round(2.9, dec), round(1.53, dec), round(1.425, dec), round(0.725, dec)}
-- vv for 2-comp mode vv --
-- dwarf_mass_l = {round(2429.198, dec), round(107.041, dec), round(80.159, dec), round(9.384, dec)}
-- dwarf_lm_ratio = {round(0.08307, dec), round(0.05949, dec), round(0.14255, dec), round(0.06729, dec)}
-- dwarf_rscale_l = {round(2.9, dec), round(1.53, dec), round(1.425, dec), round(0.725, dec)}
-- dwarf_lr_ratio = {round(0.2, dec), round(0.2, dec), round(0.2, dec), round(0.2, dec)}

dwarf_mass_d = {}
dwarf_rscale_d = {}
if(modelComponents == 2) then
    for d = 1, ndwarfs do
        dwarf_mass_d[d] = dwarf_mass_l[d]/dwarf_lm_ratio[d]-dwarf_mass_l[d]
        dwarf_rscale_d[d] = dwarf_rscale_l[d]/dwarf_lr_ratio[d]-dwarf_rscale_l[d]
    end
end
-- dwarf_mass_d = {round(26813.5851948958, dec), round(1692.2698085392, dec), round(482.1629922834, dec), round(130.0720855996, dec)}
-- dwarf_rscale_d = {round(14.5, dec), round(7.65, dec), round(7.125, dec), round(3.625, dec)}

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

----------------------
-- DWARF GENERATION --
----------------------

-- ripped out of for_developers.lua
function soft()
    if (modelComponents == 1) then --plugs in two-comp. analog for single-comp. run so i don't have to edit the eps2 function
        sp = calculateEps2(nbodies, dwarf_rscale_l[1], dwarf_rscale_l[1], dwarf_mass_l[1]/2, dwarf_mass_l[1]/2, 0)
    else
        sp = calculateEps2(nbodies, dwarf_rscale_l[1], dwarf_rscale_d[1], dwarf_mass_l[1], dwarf_mass_d[1], 0)
    end
    return sp
end


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


potential = Potential.create{
    spherical = Spherical.hernquist{ mass  = 1.52954402e5, scale = 0.7 },
    disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
    disk2     = Disk.none{ mass = 3.0e5 },
    halo      = Halo.logarithmic{ vhalo = 74.61, scaleLength = 12.0, flattenZ = 1.0 }
}

printf(string.format("Simulating %d dwarfs with %d bodies each...\n", ndwarfs, nbodies))

function makeDwarfs(ctx, potential)
    local firstModel = {}

    if (LMC_presence) then
        -- local potential = potential
        -- local position = lbrToCartesianTable(ctx, Vector.creates(orbit_l, orbit_b, orbit_r))
        -- local velocity = Vector.creates(orbit_vx, orbit_vy, orbit_vz)
        -- local LMCposition = Vector.create(-1.1, -41.1, -27.9)
        -- local LMCvelocity = Vector.create(-57, -226, 221)
        -- local LMCmass = LMC_Mass
        -- local LMCscale = LMC_scaleRadius
        -- local LMCDynaFric = LMC_DynamicalFriction and 1 or 0
        -- local coulomb_log = CoulombLogarithm
        -- local ftime = evolveTime
        -- local tstop = revOrbTime
        -- local dt = ctx.timestep / 10.0      
        -- local masses    = dwarfMass  
        -- local rscales = rscale_t
        assert(1 == 0, "work in progress - please set LMC_presence to false!")
    else
        position  = lbrToCartesianTable(ctx, Vector.creates(orbit_l, orbit_b, orbit_r))
        velocity  = Vector.creates(orbit_vx, orbit_vy, orbit_vz)
    end
    
    if(modelComponents == 2) then 
        for i = 1, ndwarfs do
            local Model = predefinedModels.mixeddwarf{
                nbody       = nbodies,
                prng        = prng,
                position    = position[i],
                velocity    = velocity[i],
                comp1       = Dwarf.plummer{mass = dwarf_mass_l[i], scaleLength = dwarf_rscale_l[i]}, -- Dwarf Options: plummer, nfw, general_hernquist
                comp2       = Dwarf.plummer{mass = dwarf_mass_d[i], scaleLength = dwarf_rscale_d[i]}, -- Dwarf Options: plummer, nfw, general_hernquist
                ignore      = true
                }
            for _, row in ipairs(Model) do
                table.insert(firstModel, row)
            end
            print(string.format("Dwarf %d bodies generation finished", i))
        end
        
    elseif(modelComponents == 1) then
        for i = 1, ndwarfs do
            local Model = predefinedModels.plummer{
                nbody       = nbodies,
                prng        = prng,
                position    = position[i],
                velocity    = velocity[i],
                mass        = dwarf_mass_l[i],
                scaleRadius = dwarf_rscale_l[i],
                ignore      = true
                }
            for _, row in ipairs(Model) do
                table.insert(firstModel, row)
            end
            --print(string.format("Dwarf %d bodies generation finished", i))
        end
    end
    return firstModel
end

st = makeDwarfs(ctx, potential)
-------------------------
-- OUTPUT VERIFICATION --
-------------------------

mass_tot = {}
pos_avg = {}
rscale_med = {}
L_tot = {}

for d = 1, ndwarfs do
    local firstParticleIndex = ((d-1)*nbodies)+1
    local finalParticleIndex = d*nbodies

    local m = 0
    local part_pos = {x={}, y={}, z={}}
    local avg_pos = {x=0, y=0, z=0}
    local p = {x={}, y={}, z={}}
    local tot_L = 0
    for i = 1, nbodies do
        particleIndex = i + (d-1)*nbodies
        local particle = st[particleIndex]
        m = m + particle.mass
        avg_pos.x = avg_pos.x + particle.position.x
        avg_pos.y = avg_pos.y + particle.position.y
        avg_pos.z = avg_pos.z + particle.position.z
        p.x[i] = particle.velocity.x*particle.mass
        p.y[i] = particle.velocity.y*particle.mass
        p.z[i] = particle.velocity.z*particle.mass
        part_pos.x[i] = particle.position.x
        part_pos.y[i] = particle.position.y
        part_pos.z[i] = particle.position.z
        if(i==nbodies) then
            avg_pos.x = avg_pos.x/nbodies
            avg_pos.y = avg_pos.y/nbodies
            avg_pos.z = avg_pos.z/nbodies
        end
    end
    mass_tot[d] = m
    pos_avg[d] = avg_pos

    --print(part_pos.x[1], part_pos.y[1], part_pos.z[1])

    local part_r = {}
    for i = 1, nbodies do
        part_r[i] = ((part_pos.x[i]-avg_pos.x)^2+(part_pos.y[i]-avg_pos.y)^2+(part_pos.z[i]-avg_pos.z)^2)^0.5
        --local L_vec = {x=1,y=1,z=1}
        local L_vec = cross(part_pos.x[i], part_pos.y[i], part_pos.z[i], p.x[i], p.y[i], p.z[i])
        tot_L = tot_L + ((L_vec.x)^2+(L_vec.y)^2+(L_vec.z)^2)^0.5
    end

    L_tot[d] = tot_L
    table.sort(part_r)

    if(nbodies %2 == 0) then
        rscale_med[d] = (part_r[nbodies/2]+part_r[nbodies/2+1])/2
    else
        rscale_med[d] = part_r[floor(nbodies/2)+1]
    end
end

-- print(mass_tot[1],mass_tot[2],mass_tot[3],mass_tot[4])
-- print(pos_avg[1].x, pos_avg[1].y,pos_avg[1].z, "\n", pos_avg[2].x, pos_avg[2].y,pos_avg[2].z,
--  "\n", pos_avg[3].x, pos_avg[3].y,pos_avg[3].z, "\n", pos_avg[4].x, pos_avg[4].y,pos_avg[4].z)
-- print(L_tot[1], L_tot[2], L_tot[3], L_tot[4])
-- print(rscale_med[1], rscale_med[2], rscale_med[3], rscale_med[4])

-- expected values --
mass_exp = {29241.283, 1799.464, 562.332, 139.458}
pos_exp = {
    x = {15.9462, 15.9872, -40.4811, -5.1567},
    y = {-37.5198, 2.3527, -49.5040, -10.0965},
    z = {-43.6455, -6.1609, -130.1761, -88.2643}
}
rscale_exp = {2.9, 1.53, 1.425, 0.725}
halfmass_rad_factor = 1.304766  -- conversion factor between half-mass* radius and scale radius (look it up)
                                -- *median (effectively)
L_exp = {463056766, 7938148.46, 16312544.4, 2585580.04}

dmass_threshold     = 1.0   --%
dpos_threshold      = 0.5   --kpc
drscale_threshold   = 1.0   --%
dL_threshold        = 1.0   --%

errstr = ""
for d=1, ndwarfs do
    print(string.format("Dwarf %d:", d)) --Δ
    delta_mass = 100*(mass_tot[d]/mass_exp[d]-1)    -- in percent
    print(string.format("Δmass:\t\t%3f%%", delta_mass))
    delta_x = pos_avg[d].x - pos_exp.x[d]           -- in kpc
    delta_y = pos_avg[d].y - pos_exp.y[d]           -- in kpc
    delta_z = pos_avg[d].z - pos_exp.z[d]           -- in kpc
    print(string.format("(Δx, Δy, Δz):\t(%f, %f, %f)", delta_x, delta_y, delta_z))
    delta_rscale = 100*(rscale_med[d]/halfmass_rad_factor/rscale_exp[d]-1)  -- in percent
    print(string.format("Δa (scale radius):\t%3f%% ", delta_rscale))
    print(string.format("L:\t\t\t%3f", L_tot[d]))
    delta_L = 100*(L_tot[d]/L_exp[d]-1)             -- in percent
    print(string.format("ΔL:\t\t\t%3f%%", delta_L))

    local dwarf_errstr = ""
    if(abs(delta_mass) > dmass_threshold) then
        dwarf_errstr = dwarf_errstr .. string.format("|Δm| = %f%% > %f%%\n", abs(delta_mass), dmass_threshold)
    end
    if(abs(delta_x) > dpos_threshold) then
        dwarf_errstr = dwarf_errstr .. string.format("|Δx| = %fkpc > %fkpc\n", abs(delta_x), dpos_threshold)
    end
    if(abs(delta_y) > dpos_threshold) then
        dwarf_errstr = dwarf_errstr .. string.format("|Δy| = %fkpc > %fkpc\n", abs(delta_y), dpos_threshold)
    end
    if(abs(delta_z) > dpos_threshold) then
        dwarf_errstr = dwarf_errstr .. string.format("|Δz| = %fkpc > %fkpc\n", abs(delta_z), dpos_threshold)
    end
    if(abs(delta_rscale) > drscale_threshold) then
        dwarf_errstr = dwarf_errstr .. string.format("|Δa| = %f%% > %f%%\n", abs(delta_rscale), drscale_threshold)
    end
    if(abs(delta_L) > dL_threshold) then
        dwarf_errstr = dwarf_errstr .. string.format("|ΔL| = %f%% > %f%%\n", abs(delta_L), dL_threshold)
    end
    if (dwarf_errstr~="") then
        errstr = errstr..string.format("Dwarf %d:\n", d)..dwarf_errstr
    end
end
assert(errstr == "", "\nParameter deviations over threshold:\n"..errstr)