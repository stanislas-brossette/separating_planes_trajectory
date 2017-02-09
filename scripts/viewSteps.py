#!/usr/bin/python

import sys
import yaml
import numpy as np
from mayavi import mlab

def plotBox(box, col, alpha):
    c = box['position']
    s = box['size']
    xmin = c[0] - s[0]/2
    xmax = c[0] + s[0]/2
    ymin = c[1] - s[1]/2
    ymax = c[1] + s[1]/2
    zmin = c[2] - s[2]/2
    zmax = c[2] + s[2]/2

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

    return mlab.triangular_mesh(x, y, z, triangles, color=col, opacity=alpha)

def plotPlane(p, initBox, boxes, finalBox, obstacles, col, alpha, sizePlane):
    nVec = np.array(p['normal']) #normal to the plane
    d = p['d'] #distance from 0 to the plane
    o = np.array(obstacles[p['obstacleBelow']]['position'])#center of the obstacle
    if 'boxAbove' in p:
        b = np.array(boxes[p['boxAbove']]['position']) #center of the box
    elif 'box0Above' in p and 'box1Above' in p:
        if p['box0Above'] == -1 and p['box1Above'] != len(boxes):
            b0 = np.array(initBox['position'])
            b1 = np.array(boxes[p['box1Above']]['position']) #center of the box
        elif p['box0Above'] != -1 and p['box1Above'] == len(boxes):
            b0 = np.array(boxes[p['box0Above']]['position']) #center of the box
            b1 = np.array(finalBox['position'])
        elif p['box0Above'] != -1 and p['box1Above'] != len(boxes):
            b0 = np.array(boxes[p['box0Above']]['position']) #center of the box
            b1 = np.array(boxes[p['box1Above']]['position']) #center of the box
        else:
            print("Something is wrong, there is no mobile box above this plane")
        b = (b0+b1)/2

    c = d*nVec #center of the plane

    center = o + (b-o)*(np.dot(c-o,nVec)/np.dot(b-o,nVec))
    # center = c

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
    mlab.triangular_mesh(x, y, z, triangles, color=col, opacity=alpha)
    mlab.quiver3d(center[0], center[1], center[2], nVec[0], nVec[1], nVec[2], color=col, opacity=alpha)
    return
    
print 'reading file ', str(sys.argv[1])

with open(str(sys.argv[1]), 'r') as stream:
    try:
        content = yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)


mobileBoxes = content['MobileBoxes']
obstacles = content['Obstacles']
planes = content['SeparatingPlanes']

# Initial and final
if 'InitialBox' in content:
    initBox = content['InitialBox']
    plotBox(initBox, (1, 1, 1), 1)

if(hasattr(mobileBoxes, '__iter__')):
    for b in mobileBoxes:
        plotBox(b, (0, 0, 1), 1)
if 'FinalBox' in content:
    finalBox = content['FinalBox']
    plotBox(finalBox, (1, 1, 1), 1)

# Obstacles
if(hasattr(obstacles, '__iter__')):
    for o in obstacles:
        plotBox(o, (1, 0, 0), 1)

# Separating Planes
if(hasattr(planes, '__iter__')):
    for p in planes:
        plotPlane(p, initBox, mobileBoxes, finalBox, obstacles, (0, 1, 0), 0.4, 0.2)

mlab.show()
