"""
SPDX-License-Identifier: MIT
Portions derived from a pyqtgraph example (c) pyqtgraph developers
Modifications (c) 2025 Stefan Costea, LeCAD-PEG

The MIT license applies to this file. See the repository root for EUPL-1.2
covering the rest of the project.
"""

import sys
import openpmd_api as opmd
import numpy as np
import scipy
import time
from collections import OrderedDict

import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph import parametertree as ptree
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets


class HEXAPIC_PlotWidget(QtWidgets.QSplitter):

    def __init__(self, series, parent=None):
        QtWidgets.QSplitter.__init__(self, QtCore.Qt.Orientation.Horizontal)
        self.ctrlPanel = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        self.addWidget(self.ctrlPanel)
        self.fieldList = QtWidgets.QListWidget()
        self.fieldList.setSelectionMode(self.fieldList.SelectionMode.ExtendedSelection)
        # self.ptree = ptree.ParameterTree(showHeader=False)
        # self.params = ptree.Parameter.create(name='params', type='group', children=[self.filter, self.colorMap])
        # self.ptree.setParameters(self.params, showTop=False)
        self.it = None
        self.plot = gl.GLViewWidget()
        self.closeEvent = exit # monkey-patching to kill app when window closed
        ## Get system dimensions
        series.read_iterations()
        Nx = series.get_attribute("Nx")
        Ny = series.get_attribute("Ny")
        series.flush()
        self.Nx = Nx
        self.Ny = Ny
        ## Add grids to the view
        # gx = gl.GLGridItem()
        # #gx.scale(1,1,1)
        # gx.setSize(Nx,Ny,1)
        # # gx.translate(10, 0, 0)
        # #self.plot.addItem(gx)
        gsf = 10 # grid scale factor
        gy = gl.GLGridItem()
        gy.setSize(Nx/gsf,Ny/gsf,1)
        gy.scale(gsf,gsf,1)
        gy.rotate(90, 0, 1, 0)
        gy.translate(-Nx/2, 0, 0)
        self.plot.addItem(gy)
        gz = gl.GLGridItem()
        gz.setSize(Nx/gsf,Ny/gsf,1)
        gz.scale(gsf,gsf,1)
        gz.rotate(90, 1, 0, 0)
        gz.translate(0, -Ny/2, 0)
        self.plot.addItem(gz)
        ax = gl.GLAxisItem()
        ax.setSize(x=Nx, y=Ny, z=Nx)
        # self.plot.addItem(ax)
        ## Add axes limit labels
        lpos = [Nx/2*1.01, Ny/2*1.01, Nx/2*1.01]
        ltext = ["0", str(Nx-1), str(Ny-1), str(Nx-1)]
        txtpos = [(-lpos[0], -lpos[0], -0.1), (-lpos[0], lpos[0], -0.1), 
                  (lpos[1], -lpos[1], -0.1),  (-lpos[0], -lpos[1], lpos[2])] 
        self.ranges = []
        for i in range(4):
            txtitem = gl.GLTextItem(pos=txtpos[i], text=ltext[i])
            self.ranges.append(txtitem)
            self.plot.addItem(txtitem)
        ## Simple surface plot example
        self.p1 = gl.GLSurfacePlotItem(shader='shaded', color=(0.5, 0.5, 1, 1))
        self.p1.translate(-Nx/2, -Ny/2, 0)
        self.plot.addItem(self.p1)
        self.p2 = gl.GLLinePlotItem(color=(1.0, 0.2, 0.2, 1), width=2, antialias=True)
        self.p2.translate(-Nx/2, -Ny/2, 0)
        self.plot.addItem(self.p2)
        self.ctrlPanel.addWidget(self.fieldList)
        # self.ctrlPanel.addWidget(self.ptree)
        self.addWidget(self.plot)
        # fg = fn.mkColor(getConfigOption('foreground'))
        # fg.setAlpha(150)
        #self.fieldList.currentRowChanged.connect(self.updatePlot)
        self.initFields = True

    def setFields(self, fields, mouseOverField=None):
        """
        Set the list of field names/units to be processed.
        
        The format of *fields* is the same as used by 
        :meth:`~pyqtgraph.widgets.ColorMapWidget.ColorMapParameter.setFields`
        """
        self.fields = OrderedDict(fields)
        self.mouseOverField = mouseOverField
        self.fieldList.clear()
        for f,opts in fields:
            item = QtWidgets.QListWidgetItem(f)
            item.opts = opts
            item = self.fieldList.addItem(item)
        
    def updatePlot(self):
        
        nsp = series.get_attribute("nsp")
        # part = {"position" : {}, "velocity" : {}}
        # part_src = []
        # for key in part.keys():
        #     for i in range(nsp):
        #         part[key][str(i)]={}
        #         for j in ["x", "y", "z"]:
        #             part_src.append(self.it.particles[str(i)][key][j])
        #             part[key][str(i)][j] = part_src[-1].load_chunk()

        if self.initFields:
            self.items = list(self.it.meshes.items())
            self.items = [x[0] for x in self.items if x[0]!='time']
            fields = []
            for x in self.items:
                fields.append((x, {'units': '-'}))
            self.setFields(fields)
            self.initFields = False
            self.xgrid = [0]
            self.ygrid = [0]
            
        alldata = []
        for x in self.items:
            alldata.append(self.it.meshes[x][opmd.Mesh_Record_Component.SCALAR].load_chunk())

        self.it.close() # complete the flush

        idx = self.fieldList.currentRow()
        if idx < 0:
            idx = 0

        data=alldata[idx]

        if self.xgrid[-1]!=np.shape(data)[0]-1 or self.ygrid[-1]!=np.shape(data)[1]-1:
            #print(self.xgrid[-1], np.shape(data)[0]-1, self.ygrid[-1], np.shape(data)[1]-1)
            self.xgrid = np.asarray(list(range(np.shape(data)[0])))
            self.ygrid = np.asarray(list(range(np.shape(data)[1])))
            self.plot.removeItem(self.p1)
            self.p1 = gl.GLSurfacePlotItem(shader='shaded', color=(0.5, 0.5, 1, 1))
            self.p1.translate(-self.Nx/2, -self.Ny/2, 0)
            self.plot.addItem(self.p1)
            self.p1.setData(x=self.xgrid, y=self.ygrid)
            self.ranges[1].setData(text=str(self.xgrid[-1]))
            self.ranges[2].setData(text=str(self.ygrid[-1]))
            self.plot.setCameraPosition(distance=self.xgrid[-1]+self.ygrid[-1])

        nmax = np.nanmax(data)
        nmin = np.nanmin(data)
        val = np.max([np.abs(nmax),np.abs(nmin)])
        if (nmax-nmin)!=0: 
            data = data/val*self.Nx/2
            self.p1.setData(z=data)
            self.ranges[3].setData(text=str(nmax))

    def get_trajectory(self, pos, pidx):
        data = self.p2.pos
        ppos = np.asarray([[pos['x'][pidx], pos['y'][pidx], pos['z'][pidx]]])
        print(ppos)
        if data is None:
            data = ppos
        else:
            data = np.append(data, ppos, axis=0)
        return data

## Create a GL View widget to display data
app = pg.mkQApp("GLSurfacePlot Example")

timer = QtCore.QTimer()

filename = sys.argv[1]
series = opmd.Series(filename, opmd.Access_Type.read_linear)

HPlot = HEXAPIC_PlotWidget(series)

HPlot.resize(1000,800)
HPlot.show()

for it in series.read_iterations():
    HPlot.it = it
    HPlot.updatePlot()
    QtGui.QGuiApplication.processEvents()
    time.sleep(0.01) # limit to 100 FPS


if __name__ == '__main__':
    pg.exec()
