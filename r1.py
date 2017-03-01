# Training script for creating 2D or 3D data with or without geometry.
# Usage:
#
#    manta ../scenes/_trainingData.py --help
#
#    manta ../scenes/_trainingData.py --dim 3

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
ap.add_argument("--numFrames", type=int, default=256)
ap.add_argument("--frameStride", type=int, default=4)
# The CFL condition states that ∆t should
# be small enough so that the maximum motion in the velocity field
# is less than the width of a grid cell.
ap.add_argument("--timeStep", type=float, default=0.1)
ap.add_argument("--addNoise", type=bool, default=True)
ap.add_argument("--addVortConf", type=bool, default=True)
ap.add_argument("--voxelPath", type=str, default="../../voxelizer/r1voxels/")
ap.add_argument("--voxelNameRegex", type=str, default=".*_(64|128)\.binvox")
ap.add_argument("--voxelLayoutBoxSize", type=int, default=64)
ap.add_argument("--datasetNamePostfix", type=str, default="")

verbose = False
addVortConf = True

x = [12, 52] 
y = [0, 63]
z = [17, 47]
resX = x[1] - x[0]
resY = y[1] - y[0]
resZ = z[1] - z[0]

# Adding dead zones in corners not necessarry if on not setting random initial vel or noise?

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
of 3 numbers. When creating an mantaflow object, you can always specify a
name parameter, which is used in debug output and in the GUI.
"""
solver = Solver(name="main", gridSize=gridSize, dim=args.dim)
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
# any other object can be created with the solver s as a parent.
flags = solver.create(FlagGrid) # The flag grid stores cell type (Fluid/Obstacle/Air).
vel = solver.create(MACGrid)
# What's velTmp used for?
velTmp = solver.create(VecGrid)  # Internally a Grid<Vec3>
pressure = solver.create(RealGrid)
density = solver.create(RealGrid)

timings = Timings()

if not verbose:
  setDebugLevel(-1)  # Disable debug printing altogether.
else:
  setDebugLevel(10)  # Print like crazy!

if (GUI):
  gui = Gui()
  gui.show(True)


vel.clear()
velTmp.clear()
pressure.clear()
density.clear()

# As most plugins expect the outmost cells of a simulation to be obstacles,
# this should always be used. 
flags.initDomain() # creates an empty box with solid boundaries
flags.fillGrid() #  marks all inner cells as fluid

# Add Model Geometry:
geom = binvox_rw.Voxels(np.zeros(inputDims), inputDims, [0, 0, 0],
  [1, 1, 1], "xyz", 'cur_geom')
VoxelUtils.create_grid_layout_stat(modelListTest, layoutBoxSize, geom, args.dim)
 
for i in range(x[0], x[1]+1):
  for j in range(y[0], y[1]+1):
    for k in range(z[0], z[1]+1):
      if geom.data[i, j, k]:
        flags.setObstacle(i-x[0], j-y[0], k-z[0])


# Noise field used to initialize velocity fields.
noise = None
gc.collect()
#if args.addNoise:
#  noise = utils.CreateNoiseField(solver)

utils.CreateRandomDensity(density)

# Random emitters.
emitters = []
numEmitters = 5

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

print('  Num Emitters: ' + str(len(emitters)))
print('  Global emitter amp: ' + str(globalEmitterAmp))
print('  Emitter min amp: ' + str(min(emitterAmps)))
print('  Emitter max amp: ' + str(max(emitterAmps)))

residue = 0
offset = 0  # You can use this to add additional frames.
for t in range(args.numFrames):
  directory = "../../data/datasets/%s/%06d" % \
      (datasetName, offset)

  if (t + 1 != args.numFrames):
    sys.stdout.write('  Simulating %d of %d\r' % (t + 1, args.numFrames))
    sys.stdout.flush()
  else:
    print('  Simulating %d of %d' % (t + 1, args.numFrames))  # keep \n char

  if not os.path.exists(directory):
    os.makedirs(directory)

  if(t == 0):
    residue = utils.InitSim(flags, vel, velTmp, noise, density, pressure,
                            utils.bWidth, utils.cgAccuracy,
                            utils.precondition, utils.cgMaxIterFac)
    if math.isnan(residue):
      print ("T=0: residue isNan")
      # Try again but with the preconditioner off.
      residue = utils.InitSim(flags, vel, velTmp, noise, density, pressure,
                              utils.bWidth, utils.cgAccuracy,
                              False, utils.cgMaxIterFac)
    if residue > utils.cgAccuracy * 10 or math.isnan(residue):
      print("WARNING: Residue (%f) has blown up before starting sim)" % \
          (residue))
      print("--> Starting a new simulation")
      break

  advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2,
                     boundaryWidth=utils.bWidth)
  advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2,
                     openBounds=True, boundaryWidth=utils.bWidth)

  # set zero normal velocity boundary condition on walls /
  # the boundary conditions at the obstacle are re-set
  setWallBcs(flags=flags, vel=vel)

  for em in emitters:
    em.update(solver.timestep, solver.timestep * t)
    em.addVelocities(vel, flags, utils.bWidth)

  setWallBcs(flags=flags, vel=vel)

  if addVortConf:
    vStrength = solver.dx()
    vorticityConfinement(vel=vel, flags=flags, strength=vStrength)

  if t % args.frameStride == 0:
    filename = "%06d_divergent.bin" % t
    fullFilename = directory + "/" + filename
    writeOutSim(fullFilename, vel, pressure, density, flags)

  setWallBcs(flags=flags, vel=vel)  # This will be "in" the model.
  residue = solvePressure(flags=flags, vel=vel, pressure=pressure,
                          cgMaxIterFac=utils.cgMaxIterFac,
                          cgAccuracy=utils.cgAccuracy,
                          precondition=utils.precondition)

  if math.isnan(residue):
    print ("Residue isNan")
    # try again but with the preconditioner off.
    residue = solvePressure(flags=flags, vel=vel, pressure=pressure,
                            cgMaxIterFac=utils.cgMaxIterFac,
                            cgAccuracy=utils.cgAccuracy,
                            precondition=False)

  if residue < utils.cgAccuracy * 10 or not math.isnan(residue):
    setWallBcs(flags=flags, vel=vel)  # This will be "in" the model.
  else:
    # If we hit maxIter, than residue will be higher than our
    # specified accuracy.  This is OK, but we shouldn't let it grow too
    # much.
    #
    # If it does grow (it happens 1 in every ~10k frames), then the best
    # thing to do is just start a new simulation rather than crashing out
    # completely. This ends up being faster than letting the solver
    # arbitrarily increase the number of iterations.
    #
    # Admittedly this is a hack. These high residue frames could be (and
    # probably are) highly correlated with ConvNet failure frames and so are
    # likely examples that we should be trying to incorporate. Unfortunately,
    # they're rare, and they result in significantly longer simulation
    # times when the solver bumps up against the max iter, that it's just
    # not worth it to include them.
    print("WARNING: Residue (%f) has blown up on frame %d" % (residue, t))
    print("--> Starting a new simulation")

    # Remove the last recorded divergent frame.
    if t % args.frameStride == 0:
      filename = "%06d_divergent.bin" % t
      fullFilename = directory + "/" + filename
      os.remove(fullFilename)

    # Break out of frame loop (and start a new sim).
    break

  if t % args.frameStride == 0:
    filename = "%06d.bin" % t
    fullFilename = directory + "/" + filename
    writeOutSim(fullFilename, vel, pressure, density, flags)

  solver.step()

  if verbose:
    timings.display()
