from abaqus import *
import sys
from odbAccess import *
import numpy as np
from abaqusConstants import *
from visualization import *

pts = np.loadtxt('points.txt') # 100 points among top elements, evenly distributed
pts = np.array(pts)
pts_len = len(pts)
pts = pts.astype(np.float32)
increments = np.linspace(9,69,7)
increments = increments.astype(np.int32)
odbPath = 'JobName.odb'

fhand1 = open('JobName1.txt','w') # lam_g1
fhand2 = open('JobName2.txt','w') # lam_1
fhand3 = open('JobName3.txt','w') # lam_g2
fhand4 = open('JobName4.txt','w') # lam_2
fhand5 = open('JobName5.txt','w') # x coord
fhand6 = open('JobName6.txt','w') # y coord
fhand7 = open('JobName7.txt','w') # z coord

o1=session.openOdb(name=odbPath)
odb = session.odbs[odbPath]
session.viewports['Viewport: 1'].setValues(displayedObject=o1)


for increment in increments:

	session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=increment)
	#session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(averageElementOutput=True, averagingThreshold=100)
	session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(averageElementOutput=False)
    	
	session.Path(name='Path-1', type=POINT_LIST, expression=(pts))
	pth = session.paths['Path-1']
	
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='SDV1', outputPosition=INTEGRATION_POINT)

	lam1g = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
	

	for ii in lam1g:	

		fhand1.write('%s '%ii[-1])
	fhand1.write('\n')

	del session.xyDataObjects['XYData-1']
	
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='SDV3', outputPosition=INTEGRATION_POINT)
	lam1 = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
	for ii in lam1:	

		fhand2.write('%s '%ii[-1])
	fhand2.write('\n')

	del session.xyDataObjects['XYData-1']

	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='SDV4', outputPosition=INTEGRATION_POINT)

	lam2g = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
	

	for ii in lam2g:	

		fhand3.write('%s '%ii[-1])
	fhand3.write('\n')

	del session.xyDataObjects['XYData-1']
	
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='SDV6', outputPosition=INTEGRATION_POINT)
	lam2 = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
	for ii in lam1:	

		fhand4.write('%s '%ii[-1])
	fhand4.write('\n')

	del session.xyDataObjects['XYData-1']
	
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='COORD', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'COORD1')) 
	xcoord = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
	for ii in xcoord:	

		fhand5.write('%s '%ii[-1])
	fhand5.write('\n')

	del session.xyDataObjects['XYData-1']
	
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='COORD', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'COORD2')) 
	ycoord = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)	
	for ii in ycoord:	

		fhand6.write('%s '%ii[-1])
	fhand6.write('\n')

	del session.xyDataObjects['XYData-1']

	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='COORD', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'COORD3')) 
	zcoord = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
	for ii in zcoord:	

		fhand7.write('%s '%ii[-1])
	fhand7.write('\n')

	del session.xyDataObjects['XYData-1']



fhand1.close()
fhand2.close()
fhand3.close()
fhand4.close()
fhand5.close()
fhand6.close()
fhand7.close()
session.odbs[odbPath].close()




