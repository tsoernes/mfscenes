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
# before (solvePressure, setWallBcs) and after:
# setInflowBcs(vel=vel,dir='xX',value=velInflow)
#
# useful examples: 
#   karman (constant inflow velocity) (what does flags.initDomain(inflow="xX") do??)
#   plume_2d (outflow boundaries; applyToGrid density)
#   waveletTurbulence 
#   freesurface (for outflow boundaries)
# 
# flags.setData(x, y, z, FlagOutflow)
#
# Kapasitet R1: 478 personer
# Figure out velocity units
#     -Speed is irrelevant of timestep
#     -How does speed scale? vel=0.2 -> 5 blocks/time_unit; vel=0.1 -> 1 block/time_unit
#
#
# Residue < 0.01 is OK


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
ap.add_argument("--timeStep", type=float, default=0.5)
ap.add_argument("--addNoise", type=bool, default=True)
ap.add_argument("--addVortConf", type=bool, default=True)
ap.add_argument("--voxelPath", type=str, default="../../voxelizer/r1voxels/")
ap.add_argument("--voxelNameRegex", type=str, default=".*_(64|128)\.binvox")
ap.add_argument("--voxelLayoutBoxSize", type=int, default=64)
ap.add_argument("--datasetNamePostfix", type=str, default="")

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

args = ap.parse_args()
print("\nUsing arguments:")
for k, v in vars(args).items():
  print("  %s: %s" % (k, v))
print("\n")

res = 64
bWidth = 1
gridSize = vec3(resX, resY, resZ)
resV = [res, res, res]

"""
First, a solver object is created. Solver objects are the parent object for
grids, particle systems and most other simulation objects. It requires
gridSize as a parameter, for which we use the custom vec3 datatype.
Most functions expecting vec3 will also accept a python tuple or sequence
of 3 numbers. 
"""
solver = FluidSolver(name="main", gridSize=gridSize, dim=args.dim)
solver.timestep = args.timeStep

modelListTest = VoxelUtils.create_voxel_file_list(
    args.voxelPath, args.voxelNameRegex)

# Next, the solver object is used to create grids. In the same way,
#   any other object can be created with the solver as a parent.
flags = solver.create(FlagGrid) # The flag grid stores cell type (Fluid/Obstacle/Air).
vel = solver.create(MACGrid)
density = solver.create(RealGrid)
pressure = solver.create(RealGrid)

gui = Gui()
gui.show(True)

#density.setConst(1)



# As most plugins expect the outmost cells of a simulation to be obstacles,
# this should always be used. 
# flags.initDomain(inflow="xX", phiWalls=phiWalls, boundaryWidth=0)
flags.initDomain(boundaryWidth=bWidth) # creates an empty box with solid boundaries
flags.fillGrid() #  marks all inner cells as fluid


# Add Model Geometry
geom = binvox_rw.Voxels(np.zeros(resV), resV, [0, 0, 0],
    [1, 1, 1], "xyz", 'cur_geom')
VoxelUtils.create_grid_layout_stat(modelListTest, resV, geom, args.dim)

#VoxelUtils.uniform_scale(2.0, geom)

# Add the Minecraft model as solid blocks
for i in range(x[0], x[1]+1):
  for j in range(y[0], y[1]+1):
    for k in range(z[0], z[1]+1):
      if geom.data[i, j, k]:
        a = i-x[0]
        b = j-y[0]
        c = k-z[0]
        flags.setObstacle(a, b, c)

# Add inflows 
velInflow = vec3(0, 1, 0)
#source = Box( parent=solver, center=(20,2,3), size=(3,1,1))
source = Box( parent=solver, p0=(18,2,2), p1=(22,3,3))
# Set Inflow(8) + Fluid(1) flag to 'source'. Don't know if those flags are optimal.
source.applyToGrid(grid=flags, value=9)

def chair_inflow(p0, p1):
  print("Modelled char in at p1=%s, p2=%s" % (p0, p1))
  cbox = Box(parent=solver, p0=p0, p1=p1)
  cbox.applyToGrid(grid=flags, value=9)
  return cbox

cVelInflow = vec3(0, -1, 0)
yzs = [
    [9,2],
    [15,4],
    [20,7],
    [25,9],
    [28,11],
    [31,12],
    [36,15],
    [41,17],
    [49,20],
    [55,23]]

cIns = []
for x_lr in [[12,15], [24,27]]:
  for yz in yzs:
    cIns.append(chair_inflow((x_lr[0], yz[0], yz[1]), (x_lr[1], yz[0]+1, yz[1]+1)))

# Model outflows

# Roof
for x1 in range(6, resX, 7): # 40
  for y1 in range(3, resY, 8): # 63
    print ("Placed outflow at %d, %d" % (x1, y1))
    hole = Box( parent=solver, center=(x1,y1,29), size=(1,1,1))
    hole.applyToGrid(grid=flags, value=FlagOutflow|FlagEmpty)


gc.collect()


for t in range(args.numFrames):
  mantaMsg('  Simulating %d of %d\r' % (t + 1, args.numFrames))

  # Add particles with velocity to inflow
  source.applyToGrid(grid=density, value=1.0)
  source.applyToGrid(grid=vel, value=velInflow)
  for cIn in cIns:
    cIn.applyToGrid(grid=density, value=1.0)
    cIn.applyToGrid(grid=vel, value=cVelInflow)
  advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)
  advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2,
                     openBounds=True, boundaryWidth=bWidth)

  resetOutflow(flags=flags, real=density)
  # set zero normal velocity boundary condition on walls /
  # the boundary conditions at the obstacle are re-set
  setWallBcs(flags=flags, vel=vel)
  
  #if addVortConf:
  #  vStrength = solver.dx()
  #  vorticityConfinement(vel=vel, flags=flags, strength=vStrength)
  #  setWallBcs(flags=flags, vel=vel)
  addBuoyancy(density=density, vel=vel, gravity=vec3(0,0,-4e-3), flags=flags)
  solvePressure(flags=flags, vel=vel, pressure=pressure)
  setWallBcs(flags=flags, vel=vel)  

  solver.step()
