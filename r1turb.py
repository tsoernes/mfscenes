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
s = FluidSolver(name="main", gridSize=gridSize, dim=args.dim)
s.timestep = args.timeStep

modelListTest = VoxelUtils.create_voxel_file_list(
    args.voxelPath, args.voxelNameRegex)

# Next, the solver object is used to create grids. In the same way,
#   any other object can be created with the solver as a parent.
flags = s.create(FlagGrid) # The flag grid stores cell type (Fluid/Obstacle/Air).
vel = s.create(MACGrid)
pressure = s.create(RealGrid, show=False)
#density = s.create(RealGrid)

k = s.create(RealGrid)
eps = s.create(RealGrid)
prod = s.create(RealGrid)
nuT= s.create(RealGrid)
strain= s.create(RealGrid)
vc=s.create(MACGrid)
temp=s.create(RealGrid)

# noise field
noise = s.create(NoiseField)
noise.timeAnim = 0

# turbulence particles
turb = s.create(TurbulenceParticleSystem, noise=noise)

# As most plugins expect the outmost cells of a simulation to be obstacles,
# this should always be used. 
# flags.initDomain(inflow="xX", phiWalls=phiWalls, boundaryWidth=0)
flags.initDomain(boundaryWidth=bWidth) # creates an empty box with solid boundaries
flags.fillGrid() #  marks all inner cells as fluid


# Add Model Geometry
geom = binvox_rw.Voxels(np.zeros(resV), resV, [0, 0, 0],
    [1, 1, 1], "xyz", 'cur_geom')
VoxelUtils.create_grid_layout_stat(modelListTest, resV, geom, args.dim)

# Add the Minecraft model as solid blocks
for x1 in range(x[0], x[1]+1):
  for y1 in range(y[0], y[1]+1):
    for z1 in range(z[0], z[1]+1):
      if geom.data[x1, y1, z1]:
        a = x1-x[0]
        b = y1-y[0]
        c = z1-z[0]
        flags.setObstacle(a, b, c)

# Add front inflow
fVelInflow = vec3(0, 1, 0)
#frontIn = Box( parent=s, p0=(16,2,2), p1=(24,3,4))
frontIn = Box( parent=s, p0=(12,2,2), p1=(28,3,4))
# Set Inflow(8) + Fluid(1) flag to inflow shape. Don't know if those flags are optimal.
frontIn.applyToGrid(grid=flags, value=9)

# Add inflows beneath chairs
def addChairInflow(p0, p1):
  print("Modelled char in at p1=%s, p2=%s" % (p0, p1))
  cbox = Box(parent=s, p0=p0, p1=p1)
  cbox.applyToGrid(grid=flags, value=9)
  return cbox

cVelInflow = vec3(0, -0.5, 0.5)
# y-start, z-start coordinates for chair inflows
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
# x-start, x-end coordinates for left and right chair inflows
#x_lrs = [[12,15], [24,27]]
x_lrs = [[11,16], [23,28]]

chairIns = []
for x_lr in x_lrs:
  for yz in yzs:
    chairIns.append(addChairInflow((x_lr[0], yz[0], yz[1]), (x_lr[1], yz[0]+1, yz[1]+1)))


# Add roof outflows
for x1 in range(6, resX, 7): # 40
  for y1 in range(3, resY, 8): # 63
    print ("Placed outflow at %d, %d" % (x1, y1))
    hole = Box( parent=s, center=(x1,y1,29), size=(1,1,1))
    hole.applyToGrid(grid=flags, value=FlagOutflow|FlagEmpty)
#setOpenBound(flags, bWidth, 'Z')

sdfgrad = obstacleGradient(flags)
sdf = obstacleLevelset(flags)
bgr = s.create(Mesh)
sdf.createMesh(bgr)


# turbulence parameters
L0 = 0.01 # turbulent (eddie) lenght scale
intensity = 0.14 # the turbulence intensity of the u-component of velocity at the inlet which is taken as 0.14 in the absence of measured values
nu = 0.00001568 # Kinematic viscosity of air at 25c (nu_t is turbulent/eddy viscosity)
mult = 0.1
prodMult = 2.5
enableDiffuse = True


if (GUI):  
	gui = Gui()
	gui.setBackgroundMesh(bgr)
	gui.show()
	sliderL0 = gui.addControl(Slider, text='turbulent lengthscale', val=L0, min=0.001, max=0.5)
	sliderMult = gui.addControl(Slider, text='turbulent mult', val=mult, min=0, max=1)
	sliderProd = gui.addControl(Slider, text='production mult', val=prodMult, min=0.1, max=5)
	checkDiff = gui.addControl(Checkbox, text='enable RANS', val=enableDiffuse)

KEpsilonBcs(flags=flags,k=k,eps=eps,intensity=intensity,nu=nu,fillArea=True)


gc.collect()


for t in range(args.numFrames):
  mantaMsg('  Simulating %d of %d\r' % (t + 1, args.numFrames))
  if (GUI):
    mult = sliderMult.get()
    K0 = sliderL0.get()
    enableDiffuse = checkDiff.get()
    prodMult = sliderProd.get()

  turb.seed(frontIn,50) # particle generation rate
  for chairIn in chairIns:
    turb.seed(chairIn, 5)
  
  turb.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4)
  turb.synthesize(flags=flags, octaves=1, k=k, switchLength=5, L0=L0, scale=mult, inflowBias=fVelInflow)
  if chairIns:
    turb.synthesize(flags=flags, octaves=1, k=k, switchLength=5, L0=L0, scale=mult, inflowBias=cVelInflow)
  #turb.projectOutside(sdfgrad)
  turb.deleteInObstacle(flags)
  
  KEpsilonBcs(flags=flags,k=k,eps=eps,intensity=intensity,nu=nu,fillArea=False)
  advectSemiLagrange(flags=flags, vel=vel, grid=k, order=1)
  advectSemiLagrange(flags=flags, vel=vel, grid=eps, order=1)
  KEpsilonBcs(flags=flags,k=k,eps=eps,intensity=intensity,nu=nu,fillArea=False)
  KEpsilonComputeProduction(vel=vel, k=k, eps=eps, prod=prod, nuT=nuT, strain=strain, pscale=prodMult)
  KEpsilonSources(k=k, eps=eps, prod=prod)
  
  if enableDiffuse:
    KEpsilonGradientDiffusion(k=k, eps=eps, vel=vel, nuT=nuT, sigmaU=10.0);


  advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2,
                     openBounds=True, boundaryWidth=bWidth)
  setWallBcs(flags=flags, vel=vel)
  #setInflowBcs(vel=vel,dir='xXyYzZ', value=fVelInflow)
  frontIn.applyToGrid(grid=vel, value=fVelInflow)
  for chairIn in chairIns:
    chairIn.applyToGrid(grid=vel, value=cVelInflow)
  #addBuoyancy(density=density, vel=vel, gravity=vec3(0,0,-4e-3), flags=flags)
  solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5)
  setWallBcs(flags=flags, vel=vel)  
  frontIn.applyToGrid(grid=vel, value=fVelInflow)
  for chairIn in chairIns:
    chairIn.applyToGrid(grid=vel, value=cVelInflow)
  #setInflowBcs(vel=vel,dir='xXyYzZ', value=fVelInflow) # has problems with bwidth

  s.step() 

