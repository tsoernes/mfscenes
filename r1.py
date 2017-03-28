# Training script for creating 3D data 
# Usage:
#    manta must be called from it's directory (build) or else there will 
#    be import errors
#
#    manta ../scenes/r1.py
#
# "I would recommend to model the inflows with obstacle flags, and a custom fixed velocity 
# into the domain. 
# Then you can define custom outflow regions in the flag grid (FlagOutflow).
# The open boundaries currently use that, you can check out "plume_2d" as an example. 
# Esp. the openBounds option will be important in the velocity advection."
# 
# So..
# How can inflows be modeled with obstacle flags?
#   Maybe by creating objects (Cylinder) like in plume_2d?
#
# Set open bound (FlagEmpty??) behind the classroom front (and under the floor)
# (Remove obstacle flags from under the seats to model inflows there)
# Set fixed velocity (maybe it's necessary to have two velocity fields --
#   either to have different velocity between front and floor inflows
#   or to ensure that the inflow air exits the inflows at the correct (normal) angle
#   and not at 45 degrees
# velInflow = vec3(0.9, 0, 0)
# vel.setConst(velInflow)
# before (solvePressure, setWallBcs) and after:
# setInflowBcs(vel=vel,dir='xX',value=velInflow)
#
# useful examples: 
#   plume_2d
#   waveletturbulence 
#   freesurface (for outflow boundaries)

import argparse
import gc
from manta import *
import os, shutil, math, sys, random
from Emitter import *
from voxel_utils import VoxelUtils
import utils
import binvox_rw
import numpy as np

ap = argparse.ArgumentParser()

# Some arguments the user might want to set.
ap.add_argument("--dim", type=int, default=3)
ap.add_argument("--numFrames", type=int, default=2000)
#ap.add_argument("--frameStride", type=int, default=4)
# The CFL condition states that ∆t (time step) should
# be small enough so that the maximum motion in the velocity field
# is less than the width of a grid cell.
#ap.add_argument("--timeStep", type=float, default=0.1)
ap.add_argument("--timeStep", type=float, default=1)
ap.add_argument("--addNoise", type=bool, default=True)
ap.add_argument("--addVortConf", type=bool, default=True)
ap.add_argument("--voxelPath", type=str, default="../../voxelizer/r1voxels/")
ap.add_argument("--voxelNameRegex", type=str, default=".*_(64|128)\.binvox")
ap.add_argument("--voxelLayoutBoxSize", type=int, default=64)
ap.add_argument("--datasetNamePostfix", type=str, default="")
ap.add_argument("--plumeScale", type=float, default=0.25)

verbose = False
addVortConf = True

# TODO: Consider scaling model up to 128. In comparison to 'simpleplume' the
# model seems somewhat coarse.  

# The minecraft model has some extra blocks and dead space on the edges that we do not want
x = [12, 52] 
y = [0, 63]
z = [17, 47]
resX = x[1] - x[0]
resY = y[1] - y[0]
resZ = z[1] - z[0]

# PS: Adding dead zones in corners should not be necessarry if 
#   not setting random initial vel or noise

args = ap.parse_args()
print("\nUsing arguments:")
for k, v in vars(args).items():
  print("  %s: %s" % (k, v))
print("\n")

res = 64
gridSize = vec3(resX, resY, resZ)
layoutBoxSize = [res, res, res]
inputDims = [res, res, res]

"""
First, a solver object is created. Solver objects are the parent object for
grids, particle systems and most other simulation objects. It requires
gridSize as a parameter, for which we use the custom vec3 datatype.
Most functions expecting vec3 will also accept a python tuple or sequence
of 3 numbers. 
"""
solver = FluidSolver(name="main", gridSize=gridSize, dim=args.dim)
solver.timestep = args.timeStep

datasetName = "output_current_3d_r1"

datasetName = datasetName + "_model"
datasetName = datasetName + args.datasetNamePostfix
print("Outputting dataset '%s'" % (datasetName))

modelList = []

print("using " + args.voxelPath + " for voxel data with pattern "
    + args.voxelNameRegex)
modelListTest = VoxelUtils.create_voxel_file_list(
    args.voxelPath, args.voxelNameRegex)

# Next, the solver object is used to create grids. In the same way,
#   any other object can be created with the solver as a parent.
flags = solver.create(FlagGrid) # The flag grid stores cell type (Fluid/Obstacle/Air).
vel = solver.create(MACGrid)
# What's velTmp used for?
#velTmp = solver.create(VecGrid)  # Internally a Grid<Vec3>
density = solver.create(RealGrid)
pressure = solver.create(RealGrid)

timings = Timings()

if not verbose:
  setDebugLevel(-1)  # Disable debug printing altogether.
else:
  setDebugLevel(10)  # Print like crazy!

gui = Gui()
gui.show(True)


# Add inflows 
source = Cylinder( parent = solver, center=vec3(20,1,20), radius=6, z=vec3(0,
                                                                            0.5, 0))
# As most plugins expect the outmost cells of a simulation to be obstacles,
# this should always be used. 
# flags.initDomain(inflow="xX", phiWalls=phiWalls, boundaryWidth=0)
flags.initDomain() # creates an empty box with solid boundaries
flags.fillGrid() #  marks all inner cells as fluid

# Where is plume positioned?
# plumeRad = 0.15  # Should match rad in fluid_net_3d_sim.
# plumeScale = args.plumeScale
# plumeUp = True
# plumeFace = 'y'
# setPlumeBound(flags, density, vel, utils.bWidth, plumeRad, plumeScale, plumeUp,
#              plumeFace)

# Add Model Geometry
geom = binvox_rw.Voxels(np.zeros(inputDims), inputDims, [0, 0, 0],
  [1, 1, 1], "xyz", 'cur_geom')
VoxelUtils.create_grid_layout_stat(modelListTest, layoutBoxSize, geom, args.dim)

# Add the Minecraft model as solid blocks
for i in range(x[0], x[1]+1):
  for j in range(y[0], y[1]+1):
    for k in range(z[0], z[1]+1):
      if geom.data[i, j, k]:
        flags.setObstacle(i-x[0], j-y[0], k-z[0])


gc.collect()

# Do we really want or need random density? prolly not
#utils.CreateRandomDensity(density)

# Random emitters.
emitters = []
numEmitters = 0

emitterAmps = []
globalEmitterAmp = random.uniform(0.1, 1)
for e in range(0, numEmitters):
  # NOTE: to test out all these emitter properties you can use
  # scenes/EmitterTest to see a visualization.

  # eRad: controls the size of the emitter.
  eRad = random.randint(1, 3)

  # eVel: speed of movement around the domain.
  eVel = 10 ** random.uniform(0.5, 0.8)

  # eAmp: the amplitude of the force applied (in the same direction as vel).
  eAmp = 10 ** random.uniform(-0.3, 0.3) * globalEmitterAmp

  # eCurvature and eCurveTScale define the amount of movement in the
  # particle's path through the simulation.
  eCurvature = 10 ** random.uniform(-1, 1)
  eCurveTScale = 10 ** random.uniform(-1, 1)

  # Create the emitter.
  emitters.append(ForceEmitter(eRad, eVel, eAmp, args.dim == 3, resX, resY,
                               resZ, eCurvature, eCurveTScale, utils.bWidth))
  emitterAmps.append(eAmp)

if emitters:
  print('  Num Emitters: ' + str(len(emitters)))
  print('  Global emitter amp: ' + str(globalEmitterAmp))
  print('  Emitter min amp: ' + str(min(emitterAmps)))
  print('  Emitter max amp: ' + str(max(emitterAmps)))

# Noise field used to initialize velocity fields.
noise = None
noise = solver.create(NoiseField, loadFromFile=True)
noise.posScale = vec3(90)
#noise.posScale = vec3(45)
noise.clamp = False
#noise.clampNeg = 0
#noise.clampPos = 1
#noise.valOffset = 0.75
#noise.timeAnim = 0.2

# TODO: Add outflows


directory = "../../data/datasets/%s/%06d" % \
    (datasetName, 0)
if not os.path.exists(directory):
    os.makedirs(directory)

residue = 0
for t in range(args.numFrames):

  if (t + 1 != args.numFrames):
    sys.stdout.write('  Simulating %d of %d\r' % (t + 1, args.numFrames))
    sys.stdout.flush()
  else:
    print('  Simulating %d of %d' % (t + 1, args.numFrames))  # keep \n char


  #if(t == 0):
  if(False):
    residue = utils.InitSim(flags, vel, velTmp, noise, density, pressure,
                            utils.bWidth, utils.cgAccuracy,
                            utils.precondition, utils.cgMaxIterFac)


  densityInflow(flags=flags, density=density, noise=noise, shape=source,
                scale=0, sigma=1.0)
  source.applyToGrid(grid=density, value=1.0)
  advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2,
                     boundaryWidth=utils.bWidth)
  advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2,
                     strength=1.0, boundaryWidth=utils.bWidth)
                     #openBounds=True, boundaryWidth=utils.bWidth)

  # set zero normal velocity boundary condition on walls /
  # the boundary conditions at the obstacle are re-set
  setWallBcs(flags=flags, vel=vel)
  #if emitters:
  #  for em in emitters:
  #    em.update(solver.timestep, solver.timestep * t)
  #    em.addVelocities(vel, flags, utils.bWidth)
  #  setWallBcs(flags=flags, vel=vel)
  #setPlumeBound(flags, density, vel, utils.bWidth, plumeRad, plumeScale,
  #              plumeUp, plumeFace)

  #if addVortConf:
  #  vStrength = solver.dx()
  #  vorticityConfinement(vel=vel, flags=flags, strength=vStrength)

  # if t % args.frameStride == 0:
  #   filename = "%06d_divergent.bin" % t
  #   fullFilename = directory + "/" + filename
  #   writeOutSim(fullFilename, vel, pressure, density, flags)

  addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-4,0), flags=flags)
  setWallBcs(flags=flags, vel=vel)
  residue = solvePressure(flags=flags, vel=vel, pressure=pressure)
  # residue = solvePressure(flags=flags, vel=vel, pressure=pressure,
  #                         cgMaxIterFac=utils.cgMaxIterFac,
  #                         cgAccuracy=utils.cgAccuracy,
  #                         precondition=utils.precondition)

  # if residue < utils.cgAccuracy * 10 or not math.isnan(residue):
  #   # This (pretty much?) always fires 
  #   setWallBcs(flags=flags, vel=vel)  # This will be "in" the model.
  setWallBcs(flags=flags, vel=vel)  # This will be "in" the model.

  # if t % args.frameStride == 0:
  #   filename = "%06d.bin" % t
  #   fullFilename = directory + "/" + filename
  #   writeOutSim(fullFilename, vel, pressure, density, flags)

  solver.step()

  if verbose:
    timings.display()
