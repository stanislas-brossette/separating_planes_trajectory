#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import arange, pi, cos, sin

from traits.api import HasTraits, Range, Instance, \
        on_trait_change
from traitsui.api import View, Item, Group

from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor, \
                MlabSceneModel

def box(posX, posY, posZ, size):
    xmin = float(posX) - float(size[0])/2
    xmax = float(posX) + float(size[0])/2
    ymin = float(posY) - float(size[1])/2
    ymax = float(posY) + float(size[1])/2
    zmin = float(posZ) - float(size[2])/2
    zmax = float(posZ) + float(size[2])/2

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
    return x, y, z, triangles

class MyModel(HasTraits):
    posX = Range(-10, 10, 5)
    posY = Range(-10, 10, 5)
    posZ = Range(-10, 10, 5)
    size = [1, 2, 3]

    scene = Instance(MlabSceneModel, ())

    plot = Instance(PipelineBase)


    # When the scene is activated, or when the parameters are changed, we
    # update the plot.
    @on_trait_change('posX,posY,posZ,scene.activated')
    def update_plot(self):
        x, y, z, t = box(self.posX, self.posY, self.posZ, self.size)
        if self.plot is None:
            self.plot = self.scene.mlab.triangular_mesh(x, y, z, t, color=(1,0,0))
        else:
            self.plot.mlab_source.set(x=x, y=y, z=z)


    # The layout of the dialog created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                Group(
                        '_', 'posX', 'posY', 'posZ'
                     ),
                resizable=True,
                )

my_model = MyModel()
my_model.configure_traits()


