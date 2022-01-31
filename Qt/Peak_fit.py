from multiprocessing.pool import IMapIterator
from re import M
import sys
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QWidget,
    QHBoxLayout, QMenu, QFileDialog
)
from PyQt6.QtGui import QPalette, QColor, QAction
from pyqtgraph import PlotWidget
import pyqtgraph as pg
#from PyQt6.QtCore import Qt


### Import non QT related stuff
sys.path.append('../src/diamond_analysis')
import numpy as np
import matplotlib.pyplot as plt
import sys 
from fit_peaks import *
from load_data import *

### import custom QT stuff
sys.path.append('resources')
from multiple_checkboxes import MultipleCheckbox

class Color(QWidget):

    def __init__(self, color):
        super(Color, self).__init__()
        self.setAutoFillBackground(True)

        palette = self.palette()
        palette.setColor(QPalette.ColorRole.Window, QColor(color))
        self.setPalette(palette)

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.initUI()
        

    def initUI(self):
        # Set the main layout as a Hbox
        hbox = QHBoxLayout()
        hbox.setSpacing(10)

        # Define the plot layout 
        widget_graph = pg.PlotWidget(background='w')
        
        hbox.addWidget(widget_graph,7)
        hbox.addWidget(QWidget(),2)

        widget = QWidget()
        widget.setLayout(hbox)
        self.setCentralWidget(widget)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('File')

        ### Set up the menu bar
        impMenu = QMenu('Load', self)
        impAct = QAction('Load data', self)
        impMenu.addAction(impAct)


        ### Setup the menu bar on-click actions
        impAct.triggered.connect(self.load_data)

        fileMenu.addMenu(impMenu)

        self.setGeometry(300, 300, 1050, 850)
        self.setWindowTitle("Peak fit")
        self.show()
    
    def load_data(self):
        path = QFileDialog.getExistingDirectory(self,"Choose Directory")
        if path != ('', ''):
            #self.dir_lineout = (path[0].rpartition('/'))[0] + '/'
            self.dir_lineout = path
            print(f'Lineout directory was set to: {self.dir_lineout}')
            self.files_drive, self.files_ambient, self.runs = list_runs(self.dir_lineout, pattern_run, pattern_bkg, pattern_ign)
            print(f'Found runs: {self.runs} in the lineout directory' )
            dlg = MultipleCheckbox(self.runs)
            dlg.setWindowTitle("Choose the runs you want to load.")
            dlg.exec()


app = QApplication(sys.argv)
w = MainWindow()
w.show()
app.exec()