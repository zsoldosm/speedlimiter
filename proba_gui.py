#!/usr/bin/env python

"""
Waypoint and map editor
"""

import subprocess
import os
import pyqtgraph as pg
import pyqtgraph.Qt as qtgqt
import pyqtgraph.dockarea as darea
import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QSlider, QWidget, QFileDialog, QInputDialog
from PyQt5 import QtCore, QtGui
#import rospy, rospkg, roslaunch
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import math

class PlotHandler():
    def __init__(self, *args, **kwargs):
        pg.setConfigOptions(antialias=True)
        self.app = QtGui.QApplication([])
        self.filename = ""
        self.ulimit = 0.0
        self.dlimit = 0.0

    
    def initializePlot(self):
        self.win = QtGui.QMainWindow()
        self.win.setUpdatesEnabled(True)
        area = darea.DockArea()
        """red = (200, 66, 66); redB = pg.mkBrush(200, 66, 66, 200)
        blue = (6, 106, 166); blueB = pg.mkBrush(6, 106, 166, 200)
        green = (16, 200, 166); greenB = pg.mkBrush(16, 200, 166, 200)
        yellow = (244, 244, 160); yellowB = pg.mkBrush(244, 244, 160, 200)"""
        self.win.setWindowTitle("Velocity editor")
        self.win.resize(1000, 800)
        self.win.setCentralWidget(area)
        dock1 = darea.Dock("Params", size = (1,1))  # give this dock minimum possible size
        dock2 = darea.Dock("Diagrams to compare", size = (500,400)) # size is only a suggestion
        area.addDock(dock1, "top")
        area.addDock(dock2, "bottom", dock1)
        widg1 = pg.LayoutWidget()
        self.csv1Label = QtGui.QLabel("none"); self.csv1Label.setAlignment(pg.QtCore.Qt.AlignCenter)
        self.usliderLabel = QtGui.QLabel("Acceleration limit"); self.usliderLabel.setAlignment(pg.QtCore.Qt.AlignCenter)
        self.dsliderLabel = QtGui.QLabel("Deceleration limit"); self.dsliderLabel.setAlignment(pg.QtCore.Qt.AlignCenter)
        self.selectFileBtn = QtGui.QPushButton("Select file")
        self.selectFileBtn.setStyleSheet("font: 10pt; color: rgb(40, 40, 40)")
        self.saveFileBtn = QtGui.QPushButton("Save file")
        self.saveFileBtn.setStyleSheet("font: 10pt; color: rgb(40, 40, 40)")
        self.update_d = QtGui.QPushButton("Generate!")
        self.uslider = QSlider(QtCore.Qt.Horizontal)
        self.uslider.setRange(0, 100)
        self.uslider.setSingleStep(1)
        self.uslidervalue = QtGui.QLabel("- m/s^2"); self.uslidervalue.setAlignment(pg.QtCore.Qt.AlignCenter)
        self.usliderLabel.setStyleSheet("font: 10pt; color: rgb(40, 40, 40)")
        self.dslider = QSlider(QtCore.Qt.Horizontal)
        self.dslider.setRange(0, 100)
        self.dslider.setSingleStep(1)
        self.dslidervalue = QtGui.QLabel("- m/s^2"); self.dslidervalue.setAlignment(pg.QtCore.Qt.AlignCenter)
        self.dsliderLabel.setStyleSheet("font: 10pt; color: rgb(40, 40, 40)")
        widg1.setStyleSheet("background-color: rgb(255, 255, 255); color: rgb(40, 40, 40);")
        dock1.setStyleSheet("background-color: rgb(255, 255, 255);")
        dock2.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.csv1Label.setStyleSheet("font: 12pt; color: rgb(40, 40, 40)")
        self.uslidervalue.setStyleSheet("font: 10pt; color: rgb(40, 40, 40)")
        self.dslidervalue.setStyleSheet("font: 10pt; color: rgb(40, 40, 40)")
        widg1.addWidget(self.selectFileBtn, row=0, col=1)
        widg1.addWidget(self.saveFileBtn, row=0, col=3)
        widg1.addWidget(self.csv1Label, row=0, col=2)
        widg1.addWidget(self.usliderLabel, row=1, col=1)
        widg1.addWidget(self.update_d, row = 1, col = 4)
        widg1.addWidget(self.uslider, row=1, col=2)
        widg1.addWidget(self.uslidervalue, row=1, col=3)
        widg1.addWidget(self.dsliderLabel, row=2, col=1)
        widg1.addWidget(self.dslider, row=2, col=2)
        widg1.addWidget(self.dslidervalue, row=2, col=3)
        dock1.addWidget(widg1)
        self.state = None
        self.widg2 = MplCanvas(self, width=5, height=4, dpi=100)
        dock2.addWidget(self.widg2)
        self.widg2.setHidden(True)
        if ((self.filename != "") and ((self.ulimit + self.dlimit) > 0.0)):
            self.update_plot()

        self.selectFileBtn.clicked.connect(self.openFileNameDialog)
        self.saveFileBtn.clicked.connect(self.saveFileDialog)
        self.textSpeedArray = np.empty(1, dtype=object)
        self.uslider.valueChanged.connect(self.uValueHandler)
        self.dslider.valueChanged.connect(self.dValueHandler)
        self.update_d.clicked.connect(self.update_plot)
        self.win.show()

    def uValueHandler(self,value):   
        uscaledValue = float(value)/20
        self.uslidervalue.setText(str(uscaledValue) + " [m/s^2]")
        self.ulimit = uscaledValue

    def dValueHandler(self,value):   
        dscaledValue = float(value)/20
        self.dslidervalue.setText(str(dscaledValue) + " [m/s^2]")
        self.dlimit = dscaledValue

    def openFileNameDialog(self):
        filename, _filter = QFileDialog.getOpenFileName(None,"Open...", filter='CSV files (*.csv)')
        if filename:
            self.csv1Label.setText(os.path.basename(str(filename)))
            self.filename = str(filename)

    def saveFileDialog(self):
        fileName, _ = QFileDialog.getSaveFileName(None,"Save to...", filter='CSV Files (*.csv)')
        if fileName:
            save_fname = str(fileName)
            lim = np.asarray(self.limiter(), order='F')
            df1 = pd.DataFrame({
                'x': self.data()['x'],
                'y': self.data()['y'],
                'z': self.data()['z'],
                'yaw': self.data()['yaw'],
                'velocity': lim,
                'change_flag': self.data()['change_flag'] 
            })
            df1.to_csv(save_fname, index=False)

    
    def update_plot(self):
        self.widg2.setHidden(False)
        self.widg2.vt1.cla()
        self.widg2.vt2.cla()
        self.widg2.at1.cla()
        self.widg2.at2.cla()
        self.plot_diags()
        self.widg2.draw_idle()
        self.widg2.showNormal()


    def plot_diags(self):
        x1 = self.time()
        y1 = self.data()["velocity"]
        x2 = self.time()
        y2 = self.limiter()
        x3 = self.time()[:-1]
        y3 = self.del_accelerate()
        x4 = self.time()[:-1]
        y4 = self.del_accelerate_l()
        self.widg2.vt1.grid(axis='both')
        self.widg2.vt2.grid(axis='both')
        self.widg2.at1.grid(axis='both')
        self.widg2.at2.grid(axis='both')
        self.widg2.vt1.plot(x1, y1, 'x-b')
        self.widg2.vt1.set_xlabel('time [sec]')
        self.widg2.vt1.set_ylabel('Original [km/h]')
        self.widg2.vt2.plot(x2, y2, 'x-g')
        self.widg2.vt2.set_xlabel('time [sec]')
        self.widg2.vt2.set_ylabel('Limited [km/h]')
        self.widg2.at1.plot(x3, y3, 'ob')
        self.widg2.at1.set_xlabel('time [sec]')
        self.widg2.at1.set_ylabel('Original [m/s^2]')
        self.widg2.at2.plot(x4, y4, 'og')
        self.widg2.at2.set_xlabel('time [sec]')
        self.widg2.at2.set_ylabel('Limited [m/s^2]')

    def data(self):
        """CSV kiterjesztésű, adott formátumú fájlból sebesség - értékek kinyerése."""

        names = ["x", "y", "z", "yaw", "velocity", "change_flag"]
        df = pd.read_csv(self.filename, sep=",", names=names, skiprows=1)
        x = len(df) - 1
        df._set_value(0, "velocity", 0.0)
        df._set_value(x, "velocity", 0.0)
        return df

    def delta_speed(self):
        """Mintavételezési pontok közötti sebesség - különbségek kiszámítása."""
        ds = []
        for i in range(1, len(self.data())):
            ds.append((self.data()["velocity"][i])-(self.data()["velocity"][i-1]))
            ds[i-1] = round(ds[i-1], 4)
        return ds

    def delta_speed_l(self):
        """Mintavételezési pontok közötti limitált sebesség - különbségek kiszámítása."""
        dsl = []
        for i in range(1, len(self.limiter())):
            dsl.append((self.limiter()[i])-self.limiter()[i-1])
            dsl[i-1] = round(dsl[i-1], 4)
        return dsl

    def delta_dist(self):
        """X, Y koordinátákból mintavételezési pontok között megtett út kiszámítása Pitagorasz - tétellel."""
        dd = []
        dat = self.data()
        for j in range(1, len(dat)):
            s1 = abs(dat["x"][j] - (dat["x"][j-1]))
            s0 = abs(dat["y"][j] - (dat["y"][j-1]))
            del_s = round(math.sqrt(s1**2 + s0**2), 4)
            dd.append(del_s)
        return dd

    def distance(self):
        """Adott pontig megtett út kiszámítása delta_dist - ből."""
        d = []
        del_d = self.delta_dist()
        for item in del_d:
            d.append(item)
        for i in range(0, len(del_d)):
            for j in range(0,i):
                d[i] += del_d[j] 
            d[i] = round(d[i], 4)
            a = [0]
        return a + d

    def delta_time(self):
        """Mintavételezési pontok között eltelt idő kiszámítása."""
        v = (self.data()["velocity"])
        v.pop(0)
        s = self.delta_dist()
        t = []
        for i in range(1, len(v)):
            x = round((s[i-1] / v[i])*3.6, 4)   #km/h -> m/s
            t.append(x)
        cor = (s[len(s)-1] / (v[len(s)-1]))*3.6
        cor = round(cor, 4)
        a = [cor]
        return t + a

    def time(self):
        """Adott pontig eltelt idő kiszámítása delta_time - ból."""
        ti = []
        del_t = self.delta_time()
        for item in del_t:
            ti.append(item)
        for i in range(0, len(del_t)):
            for j in range(0,i):
                ti[i] += del_t[j] 
            ti[i] = round(ti[i], 4)
        a = [0]
        return a + ti

    def del_accelerate(self):
        """Mintavételezési pontok közötti gyorsulás - különbségek kiszámítása."""
        dv = self.delta_speed()
        dt = self.delta_time()
        acc = []
        for i in range(0,len(dv)):
            x = (dv[i]/3.6) / dt[i]
            x = round(x, 4)
            acc.append(x)
        return acc

    def del_accelerate_l(self):
        """Mintavételezési pontok közötti limitált gyorsulás - különbségek kiszámítása."""
        dv = self.delta_speed_l()
        dt = self.delta_time()
        accl = []
        for i in range(0,len(dv)):
            x = (dv[i]/3.6) / dt[i]
            x = round(x, 4)
            accl.append(x)
        return accl

    def limiter(self):
        """Egymás utáni pontok gyorsulás értékeinek összehasonlítása, limitálása pozitív és negatív irányból."""
        v = self.data()['velocity']
        dt = self.delta_time()
        v_lim = []
        for item in v:
            v_lim.append(item)
        for i in range(0, len(v)-1):
            if (((v[i+1] - v_lim[i])/dt[i]) >= (self.ulimit * 3.6)):
                v_lim[i+1] = round((v_lim[i] + (self.ulimit * 3.6 * dt[i])), 4)
            else: v_lim[i+1] = v[i+1]

        v_lim_2 = []
        for item in v_lim:
            v_lim_2.append(item)
        v_lim_2[-1] = 0

        for j in range(len(v)-1, 0, -1):
            if ((v_lim[j-1] - v_lim_2[j]) >= (self.dlimit * 3.6 * dt[j-1])):
                v_lim_2[j-1] = round((v_lim_2[j] + (self.dlimit * 3.6 * dt[j-1])), 4)
            else: v_lim_2[j-1] = v_lim[j-1]

        return v_lim_2


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        fig.suptitle('Original vs. Limited Velocity and Acceleration')
        self.vt1 = fig.add_subplot(2,2,1)
        self.vt2 = fig.add_subplot(2,2,3)
        self.at1 = fig.add_subplot(2,2,2)
        self.at2 = fig.add_subplot(2,2,4)
        super(MplCanvas, self).__init__(fig)

if __name__ == "__main__":
    print("GUI started")
    ph = PlotHandler()
    ph.initializePlot()

    if (sys.flags.interactive != 1) or not hasattr(QtCore, "PYQT_VERSION"):
        QApplication.instance().exec_()