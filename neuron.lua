--[[!
-- \file neuron.lua
-- 
-- \brief example script how to use NEURON API within membrane potential mapping plugin
--
--]]

-- script name 
scriptname = debug.getinfo(1).short_src

-- load ug script util
ug_load_script("../scripts/ug_util.lua")

-- check if plugins are loaded
AssertPluginsLoaded( {"MembranePotentialMapping"} )
AssertPluginsLoaded( {"PlasmaMembrane"} )

-- hoc setup
hoc_stim = util.GetParam("-hoc_stim", "/Users/stephangrein/git/ug4/build/stim.hoc")
hoc_geom = util.GetParam("-hoc_geom", "/Users/stephangrein/git/ug4/build/geom.hoc")
hoc_dt = util.GetParamNumber("-hoc_dt", 0.1)
hoc_tstop = util.GetParamNumber("-hoc_tstop", 1.0)
hoc_finitialize = util.GetParamNumber("-hoc_finitialize", -75.0)

-- create neuron instance
HocSetup = Transformator()

-- load geometry and stimulation protocol hoc files
HocSetup:load_geom(hoc_geom)
HocSetup:load_stim(hoc_stim)

-- setup hoc interpreter (tstart, tstop, dt, initial membrane potential)
HocSetup:setup_hoc(0, hoc_tstop, hoc_dt, hoc_finitialize)
-- print summarzing table of setup
HocSetup:print_setup(true)

-- purge the hoc interpreter and load another geometry/stimulation protocol
HocSetup:purge()

-- create mpm instances for current and next timestep
before = os.clock()
mpm_t_current = MembranePotentialMapper(HocSetup)
after = os.clock()

-- measure time for cosntructing of the tree with no points above
print()
print("Elapsed time to build current tree instance: " .. math.abs(after-before) .. " [s]")
print()
before = os.clock()
HocSetup:extract_vms(1, 3) -- advances 3 timesteps and extract the membrane potentials at the arriving timestep
after = os.clock()

-- measure elapsed time (most basic) during extract vms by neuron api
print()
print("Elapses time during extracting timesteps: " ..  math.abs(after-before) .. " [s]")
print()



before = os.clock()
mpm_t_next = MembranePotentialMapper(HocSetup)
after = os.clock()

-- measure time for cosntructing of the tree with no points above
print()
print("Elapsed time to build next tree instance: " .. math.abs(after-before) .. " [s]")
print()

before = os.clock()
print("Potential for x:=(1,2,3): " .. mpm_t_current:get_potential(1, 2, 3) .. " [mV]")
print("Potential for x:=(1,2,3): " .. mpm_t_next:get_potential(1, 2, 3) .. " [mV]")
after = os.clock()

-- measure elapsed time (most basic) during mapping of one point in two vm2ug instances
print()
print("Elapses time during mapping of one point in two vm2ug instances: " ..  math.abs(after-before) .. " [s]")
print()


HocSetup:print_setup(true)

t_current = 0
-- some neumann boundary
function NeumannBoundary(x, y, z, t, si) 
   if (t == t_current) then
      return mpm_t_current:get_potential(x, y, z)
   elseif (t == t_current + dt) then
      return mpm_t_next:get_potetntial(x, y, z)
   elseif (t == t_current + 2.0 * dt) then
      t_current = t_current + dt
      mpm_t_current = mpm_t_next;
      HocSetup:extract_vms(1)
      mpm_t_next:build_tree()
      return mpm_t_next:get_potential(x, y, z)
   end   
end

-- define vars for ug algebra
dim = 3
system = 1
algebra = "CPU"

-- init ug with an algebra type (CPU for systems of PDEe: 1)
InitUG(dim, AlgebraType(algebra, system));

