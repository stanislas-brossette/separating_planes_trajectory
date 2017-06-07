#!/usr/bin/python

import sys
import yaml
import numpy as np
from mayavi import mlab

from traits.api import HasTraits, Range, Instance, \
        on_trait_change
from traitsui.api import View, Item, Group

from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor, \
                MlabSceneModel

def lenWithNone(l):
    if l is None:
        return 0
    else:
        return len(l)


def plotBox(triangleMesh, box, col, alpha):
    c = box['position']
    s = box['size']
    xmin = float(c[0]) - float(s[0])/2
    xmax = float(c[0]) + float(s[0])/2
    ymin = float(c[1]) - float(s[1])/2
    ymax = float(c[1]) + float(s[1])/2
    zmin = float(c[2]) - float(s[2])/2
    zmax = float(c[2]) + float(s[2])/2

    x = [xmin, xmax, xmax, xmin, xmin, xmax, xmax, xmin] 
    y = [ymin, ymin, ymax, ymax, ymin, ymin, ymax, ymax]
    z = [zmin, zmin, zmin, zmin, zmax, zmax, zmax, zmax] 
    triangles = [
            (0, 1, 3), (1,2,3),
            (0, 1, 4), (1,4,5),
            (1, 2, 5), (2,5,6),
            (2, 3, 6), (3,6,7),
            (0, 3, 4), (3,4,7),
            (4, 5, 6), (4,6,7)
            ]

    if (triangleMesh == None):
        triangleMesh = mlab.triangular_mesh(x, y, z, triangles, color=col, opacity=alpha)
    else:
        triangleMesh.mlab_source.set(x=x,y=y,z=z)
    return triangleMesh

def plotPlane(triangleMesh, quiver, p, initBox, boxes, finalBox, obstacles, sizePlane, col, alpha):
    nVec = np.array(map(float,p['normal'])) #normal to the plane
    d = p['d'] #distance from 0 to the plane
    o = np.array(map(float,obstacles[p['obstacleBelow']]['position']))#center of the obstacle
    if 'boxAbove' in p:
        b = np.array(map(float,boxes[p['boxAbove']]['position'])) #center of the box
    elif 'box0Above' in p and 'box1Above' in p:
        if p['box0Above'] == -1 and p['box1Above'] != lenWithNone(boxes):
            b0 = np.array(map(float,initBox['position']))
            b1 = np.array(map(float,boxes[p['box1Above']]['position'])) #center of the box
        elif p['box0Above'] != -1 and p['box1Above'] == lenWithNone(boxes):
            b0 = np.array(map(float,boxes[p['box0Above']]['position'])) #center of the box
            b1 = np.array(map(float,finalBox['position']))
        elif p['box0Above'] != -1 and p['box1Above'] != lenWithNone(boxes):
            b0 = np.array(map(float,boxes[p['box0Above']]['position'])) #center of the box
            b1 = np.array(map(float,boxes[p['box1Above']]['position'])) #center of the box
        else:
            print("Something is wrong, there is no mobile box above this plane")
        b = (b0+b1)/2

    c = d*nVec #center of the plane

    # center = o + (b-o)*(np.dot(c-o,nVec)/np.dot(b-o,nVec))
    center = c

    if (nVec[1] == 0 and nVec[2] == 0):
        tVec= np.array([0, nVec[0], 0])
    else:
        tVec= np.array([0, nVec[2], -nVec[1]])

    tVec = tVec/np.linalg.norm(tVec)
    bVec = np.cross(nVec,tVec)
    p0 = center + (sizePlane/2)*tVec + (sizePlane/2)*bVec
    p1 = center + (sizePlane/2)*tVec - (sizePlane/2)*bVec
    p2 = center - (sizePlane/2)*tVec + (sizePlane/2)*bVec
    p3 = center - (sizePlane/2)*tVec - (sizePlane/2)*bVec
    x = [p0[0],p1[0],p2[0],p3[0]]
    y = [p0[1],p1[1],p2[1],p3[1]]
    z = [p0[2],p1[2],p2[2],p3[2]]
    triangles = [(0, 1, 2), (1,2,3)]
    if (triangleMesh == None):
        triangleMesh = mlab.triangular_mesh(x, y, z, triangles, color=col, opacity=alpha)
        quiver = mlab.quiver3d(center[0], center[1], center[2], nVec[0], nVec[1], nVec[2], color=col, opacity=alpha)
    else:
        triangleMesh.mlab_source.set(x=x,y=y,z=z)
        quiver.mlab_source.set(x=center[0], y=center[1], z=center[2], u=nVec[0], v=nVec[1], w=nVec[2])
    return triangleMesh, quiver

def plotFixedPlane(triangleMesh, quiver, p, sizePlane, col, alpha):
    nVec = np.array(p['normal']) #normal to the plane
    d = p['d'] #distance from 0 to the plane

    center = d*nVec #center of the plane

    if (nVec[1] == 0 and nVec[2] == 0):
        tVec= np.array([0, nVec[0], 0])
    else:
        tVec= np.array([0, nVec[2], -nVec[1]])

    tVec = tVec/np.linalg.norm(tVec)
    bVec = np.cross(nVec,tVec)
    p0 = center + (sizePlane/2)*tVec + (sizePlane/2)*bVec
    p1 = center + (sizePlane/2)*tVec - (sizePlane/2)*bVec
    p2 = center - (sizePlane/2)*tVec + (sizePlane/2)*bVec
    p3 = center - (sizePlane/2)*tVec - (sizePlane/2)*bVec
    x = [p0[0],p1[0],p2[0],p3[0]]
    y = [p0[1],p1[1],p2[1],p3[1]]
    z = [p0[2],p1[2],p2[2],p3[2]]
    triangles = [(0, 1, 2), (1,2,3)]
    # scene.mlab.quiver3d(center[0], center[1], center[2], nVec[0], nVec[1], nVec[2], color=col, opacity=alpha)
    if (triangleMesh == None):
        triangleMesh = mlab.triangular_mesh(x, y, z, triangles, color=col, opacity=alpha)
        quiver = mlab.quiver3d(center[0], center[1], center[2], nVec[0], nVec[1], nVec[2], color=col, opacity=alpha)
    else:
        triangleMesh.mlab_source.set(x=x,y=y,z=z)
        quiver.mlab_source.set(x=center[0], y=center[1], z=center[2], u=nVec[0], v=nVec[1], w=nVec[2])
    return triangleMesh, quiver

def colorIndex(i, N):
    val = (float(i)%float(N))/float(N)
    return (0, val, 1-val)

def colorObstacles(isVirtual):
    if isVirtual:
        return (1, 0.5, 0)
    else:
        return(1, 0, 0)


class MyModel(HasTraits):
    print('Python command: ', str(sys.argv))
    print 'reading file ', str(sys.argv[1])

    with open(str(sys.argv[1]), 'r') as stream:
        try:
            content = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    plotMobilePlanes = True;

    if len(sys.argv) > 2:
        if str(sys.argv[2]) in ['False', 'false', '0']:
            plotMobilePlanes = False

    nIter = content['nIter']
    iter = Range(0,nIter-1,1)

    initBox = content['InitialBox']
    finalBox = content['FinalBox']
    obstacles = content['Obstacles']
    fixedPlanes = content['FixedPlanes']
    mobileBoxes = []
    planes = []
    for i in range(0,nIter,1):
        mobileBoxes.append(content['MobileBoxes'+str(i)])
        planes.append(content['SeparatingPlanes'+str(i)])

    nFixedPlanes = 0
    if fixedPlanes is not None:
        nFixedPlanes = lenWithNone(fixedPlanes)
    nObstacles = 0
    if obstacles is not None:
        nObstacles = lenWithNone(obstacles)


    initBoxPlot = None
    finalBoxPlot = None
    fixedPlanesPlot = [None]*nFixedPlanes
    fixedPlanesQuiver = [None]*nFixedPlanes
    obstacleBoxesPlot = [None]*nObstacles

    mobileBoxesPlot = [None]*lenWithNone(mobileBoxes[0])
    mobilePlanesPlot = [None]*lenWithNone(planes[0])
    mobilePlanesQuiver = [None]*lenWithNone(planes[0])

    # When the scene is activated, or when the parameters are changed, we
    # update the plot.

    scene = Instance(MlabSceneModel, ())

    plot = Instance(PipelineBase)

    @on_trait_change('scene.activated')
    def create_plot(self):
        iterBoxes = self.mobileBoxes[self.iter]
        iterPlanes = self.planes[self.iter]
        if self.plot is None:
            # Initial and final
            self.initBoxPlot = plotBox(self.initBoxPlot, self.initBox, (0, 0, 1), 1)
            self.finalBoxPlot = plotBox(self.finalBoxPlot, self.finalBox, (0, 1, 0), 1)
            # Fixed Planes
            for i in range(0,self.nFixedPlanes,1):
                self.fixedPlanesPlot[i], self.fixedPlanesQuiver[i] = plotFixedPlane(self.fixedPlanesPlot[i], self.fixedPlanesQuiver[i], self.fixedPlanes[i], 2.0, (1, 1, 0), 1)
            # Obstacles
            for i in range(0,self.nObstacles,1):
                self.obstacleBoxesPlot[i] = plotBox(self.obstacleBoxesPlot[i], self.obstacles[i], colorObstacles(self.obstacles[i]["isVirtual"]), 0.6)
            # Mobile Boxes
            for i in range(0,lenWithNone(iterBoxes),1):
                self.mobileBoxesPlot[i] = plotBox(self.mobileBoxesPlot[i], iterBoxes[i], colorIndex(i, lenWithNone(iterBoxes)), 1)
            # Separating Planes
            if self.plotMobilePlanes:
                for i in range(0,lenWithNone(iterPlanes),1):
                    self.mobilePlanesPlot[i], self.mobilePlanesQuiver[i] = plotPlane(self.mobilePlanesPlot[i], self.mobilePlanesQuiver[i], iterPlanes[i], self.initBox, iterBoxes, self.finalBox, self.obstacles, 0.2, colorIndex(i, lenWithNone(iterBoxes)), 0.4)
            self.scene.mlab.view(azimuth=90, elevation=90, distance=4, focalpoint=[0,0,-0.4], roll=None, reset_roll=True, figure=None)

    @on_trait_change('iter')
    def update_plot(self):
        iterBoxes = self.mobileBoxes[self.iter]
        iterPlanes = self.planes[self.iter]
        for i in range(0,lenWithNone(iterBoxes),1):
            self.mobileBoxesPlot[i] = plotBox(self.mobileBoxesPlot[i], iterBoxes[i], colorIndex(i, lenWithNone(iterBoxes)), 1)
        if self.plotMobilePlanes:
            for i in range(0,lenWithNone(iterPlanes),1):
                self.mobilePlanesPlot[i], self.mobilePlanesQuiver[i] = plotPlane(self.mobilePlanesPlot[i], self.mobilePlanesQuiver[i], iterPlanes[i], self.initBox, iterBoxes, self.finalBox, self.obstacles, 0.2, colorIndex(i, lenWithNone(iterBoxes)), 0.4)


    # The layout of the dialog created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                Group('iter'),
                resizable=True,
                )

my_model = MyModel()
my_model.configure_traits()


