#
# Simple example scene (hello world)
# Simulation of a buoyant smoke density plume (with noise texture as smoke source)
#

#import pdb; pdb.set_trace()

from manta import *
import utils

# solver params
res = 64
gs  = vec3(res, int(1.5*res), res)
s   = FluidSolver(name='main', gridSize = gs)

# prepare grids
flags    = s.create(FlagGrid) # The flag grid stores cell type (Fluid/Obstacle/Air).
vel      = s.create(MACGrid)
velTmp   = s.create(VecGrid)
density  = s.create(RealGrid)
pressure = s.create(RealGrid)

# noise field, tweak a bit for smoke source
noise = s.create(NoiseField, loadFromFile=True)
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

# a cylinder shape is created, which will be used later as an inflow for smoke.
# center = (32, 6.4, 32), radius = 9, z = (0, 1.3, 0)
source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))

# initDomain creates an empty box with solid boundaries. 
# As most plugins expect the outmost cells of a simulation 
# to be obstacles, this should always be used.
flags.initDomain() 
# fillGrid then marks all inner cells as fluid.
flags.fillGrid()

if (GUI):
  gui = Gui()
  gui.show()
  
#main loop
for t in range(250):

  if(t == 0):
    residue = utils.InitSim(flags, vel, velTmp, noise, density, pressure,
                            utils.bWidth, utils.cgAccuracy,
                            utils.precondition, utils.cgMaxIterFac)
  #mantaMsg('\nFrame %i' % (s.frame))
  # Use the cylinder as an inflow for source the first 100 timesteps
  if t<100:
    densityInflow(flags=flags, density=density, noise=noise, 
                  shape=source, scale=1, sigma=0.5)
  # optionally, enforce inflow velocity
  #  The cylinder shape created above is projected to grid, 
  # and all cells within are assigned a smoke density of 1
  source.applyToGrid(grid=vel, value=vec3(0.1,0,0))

  # Then, the density and velocity grid are advected using 
  # second-order semi-Lagrangian advection, i.e. MacCormack advection.
  advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
  advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2, strength=1.0)
  
  # the boundary conditions at the obstacle are re-set
  setWallBcs(flags=flags, vel=vel)    
  # buoyant forces are added to the velocity grid
  addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-4,0), flags=flags)
  
  # the pressure projection is applied
  solvePressure( flags=flags, vel=vel, pressure=pressure )
  # The step function now tells the solver that one iteration is complete, 
  # which is important for the GUI and timing functions.
  s.step()

