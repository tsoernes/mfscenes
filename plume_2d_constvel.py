#
# Simple example scene for a 2D simulation
# Simulation of a buoyant smoke density plume with open boundaries at top & bottom
#
from manta import *

if (GUI):
	gui = Gui()
	gui.show( True ) 

# solver params
res = 64
gs = vec3(res,res,1)
s = Solver(name='main', gridSize = gs, dim=2)
s.timestep = 0.5
timings = Timings()

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)


# Is it necessary to add initial pressure/density????
#density.setConst(1)

bWidth=1
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()

# Model outflows
for i in range(1,63):
	flags.setData(i,63,0,FlagOutflow|FlagEmpty)
	flags.setData(i,62,0,FlagOutflow|FlagEmpty)

# setOpenBound and source.applyToGrid(grid=flags, FlagOutflow|FlagEmpty) 
#   should be equivalent, though perhaps not at the edges where an open cell
#   has an obstacle neighbor
#setOpenBound(flags, bWidth,'Y') 

# Model inflows
# How come, that when the inflow is turned off after X seconds, that the velocity is 
# higher in the stream vacuum than it was with the inflow turned on?
velInflow = vec3(0,1,0)
source = s.create( Box, center=(32, 1, 0.5), size=(10, 1, 0))
# Need to revise which flag is used here.. Inflow, empty..????
# it would probably affect setInflowBcs. Find out which flags 
# are used in the examples for inflows, especially for inflows 
# at the boundary
source.applyToGrid(grid=flags, value=9)
# probably not necessary since its done in loop: source.applyToGrid(grid=vel, value=velInflow)

# Add obstacle "walls" besides the inflows
# Inflows should be flush with the walls
#for i in range(1,22):
#	flags.setData(i,1,0, FlagObstacle) 
#for i in range(42,64):
#	flags.setData(42,1,0, FlagObstacle) 

def inOutBcs():
	# Ask on MantaFlow forum: 
	# What does 'dir' mean in this context? 
        # Does this work with inflows/outflows not at the grid boundaries? 
	setInflowBcs(vel=vel, dir="Y", value=velInflow)
 
steps = 2800
for t in range(steps):
	mantaMsg('\nFrame %i' % (s.frame))

	# The order of the steps below need to be verified
	#      Check examples and FluidForTheRestOfUS (smoke not liquid!)

	if t<int(steps/2):
		# Add particles to inflow
		source.applyToGrid(grid=density, value=1)
		# Set particle velocity to constant at inflow 
		source.applyToGrid(grid=vel, value=velInflow)


	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2) 
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, openBounds=True, boundaryWidth=bWidth)

	# Enforce a constant inflow/outflow at the grid boundaries
	#   This sets constant velocity (1) to outflows, which does not seem
	#   suitable for our model.  
	#inOutBcs()

	# Ensure empty flag in outflow cells, remove fluid particles and density
	resetOutflow(flags=flags,real=density) 

	# Set zero normal velocity boundary condition on walls (no-slip condition)
	# Does this also apply to obstacles that are not at the edge of the grid? 
	setWallBcs(flags=flags, vel=vel)    

	#inOutBcs()

	# Perform pressure projection of the velocity grid 
	solvePressure(flags=flags, vel=vel, pressure=pressure)
	
	#timings.display()    
	s.step()

