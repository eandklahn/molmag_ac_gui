#third-party packages
import scipy.constants as sc

from PyQt5.QtWinExtras import QWinTaskbarButton
from PyQt5.QtGui import QIcon, QFont, QDoubleValidator
from PyQt5.QtWidgets import (QMainWindow, QWidget, QApplication, QPushButton, QLabel, QAction, QComboBox, QStackedWidget,
                             QDoubleSpinBox, QFormLayout, QCheckBox, QVBoxLayout, QMessageBox, QSplitter, QGridLayout,
                             QHBoxLayout, QFileDialog, QDialog, QLineEdit, QListWidget, QListWidgetItem, QTabWidget,
                             QScrollArea, QStatusBar)

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT as NavigationToolbar

#set constants
kB = sc.Boltzmann

class PlottingWindow(QWidget):

    def __init__(self, make_cax=False):
    
        super(PlottingWindow, self).__init__()
        
        self.layout = QVBoxLayout()
        
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.tools = NavigationToolbar(self.canvas, self)
        if make_cax:
            self.grid = plt.GridSpec(20,1)
            self.ax = self.fig.add_subplot(self.grid[:17,0])
            self.cax = self.fig.add_subplot(self.grid[-1,0])
        else:
            self.ax = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.95)
        
        self.layout.addWidget(self.canvas)
        
        self.tool_lo = QHBoxLayout()
        self.tool_lo.addWidget(self.tools)
        self.tool_lo.addStretch()
        
        self.reset_axes_btn = QPushButton('Reset axes')
        self.reset_axes_btn.clicked.connect(self.reset_axes)
        self.tool_lo.addWidget(self.reset_axes_btn)
        self.layout.addLayout(self.tool_lo)
        
        self.setLayout(self.layout)
    
    def clear_canvas(self):
        
        self.ax.clear()
        self.canvas.draw()
    
    def reset_axes(self):
       
       s = 0
       if len(self.ax.lines)<1: pass
       else:
           while True:
               start = self.ax.lines[s]
               if len(start.get_xdata())<1:
                   s += 1
               else:
                   break
           
           x = start.get_xdata()
           y = start.get_ydata()
           
           new_x = [x.min(), x.max()]
           new_y = [y.min(), y.max()]
           
           for i in range(s, len(self.ax.lines)):
               x = self.ax.lines[i].get_xdata()
               y = self.ax.lines[i].get_ydata()
               
               if len(x)>1 and len(y)>1:
                   if x.min()<new_x[0]: new_x[0] = x.min()
                   if x.max()>new_x[1]: new_x[1] = x.max()
                   if y.min()<new_y[0]: new_y[0] = y.min()
                   if y.max()>new_y[1]: new_y[1] = y.max()
           
           if new_x[0] == new_x[1]:
               new_x[0] -= 0.5
               new_x[1] += 0.5
           if new_y[0] == new_y[1]:
               new_y[0] -= 0.5
               new_y[1] += 0.5
               
           self.ax.set_xlim(new_x[0]-0.05*(new_x[1]-new_x[0]),new_x[1]+0.05*(new_x[1]-new_x[0]))
           self.ax.set_ylim(new_y[0]-0.05*(new_y[1]-new_y[0]),new_y[1]+0.05*(new_y[1]-new_y[0]))
           
           self.canvas.draw()
                            
class GuessDialog(QDialog):

    def __init__(self,
                 parent=None,
                 guess=None):
        super(GuessDialog, self).__init__()
        
        self.layout = QFormLayout()
        self.setWindowTitle('Guess parameters')
        self.validator = QDoubleValidator()
        
        self.values = []
        
        for pair in guess.items():
            lbl = QLabel(pair[0])
            text = QLineEdit()
            text.setValidator(self.validator)
            if pair[0]=='Ueff':
                text.setText(str(pair[1]/kB))
            else:
                text.setText(str(pair[1]))
            self.layout.addRow(lbl, text)
            self.values.append((lbl, text))
        
        self.setLayout(self.layout)
        
        accept_btn = QPushButton('Ok')
        accept_btn.clicked.connect(self.on_close)
        self.layout.addRow(accept_btn)
        
        self.show()
        
    def on_close(self):
        
        self.return_guess = {v[0].text(): float(v[1].text()) for v in self.values}
        self.return_guess['Ueff'] = self.return_guess['Ueff']*kB
        self.accept()

class SimulationDialog(QDialog):

    def __init__(self,
                 parent=None,
                 fitted_parameters=None,
                 plot_type_list=[],
                 plot_parameters={'tQT': 0.01, 'Cr': 0.00, 'n': 0.00, 't0': 0.00, 'Ueff': 0.00},
                 min_and_max_temps=[0]*2):
    
        super(SimulationDialog, self).__init__()
        
        self.setWindowTitle('Add simulation')
        
        # Containers for objects
        self.lineedit_inputs = {}
        
        self.headline_font = QFont()
        self.headline_font.setBold(True)
        
        self.plot_type_list = plot_type_list
        self.plot_parameters = plot_parameters
        
        self.min_and_max_temps = min_and_max_temps
        
        self.layout = QVBoxLayout()
        
        # Abstracting the validator for the QLineEdits
        self.validator = QDoubleValidator()
        self.validator.setNotation(QDoubleValidator.ScientificNotation)
        
        if fitted_parameters is not None:
            self.use_fit_cb = QCheckBox('Use fitted parameters')
            self.use_fit_cb.stateChanged.connect(self.fit_take_control)
            self.layout.addWidget(self.use_fit_cb)
            
            self.fitted_parameters = {}
            for quant in fitted_parameters['quantities']:
                try:
                    idx = fitted_parameters['quantities'].index(quant)
                    param = fitted_parameters['params'][idx]
                except IndexError:
                    param = None
                finally:
                    if quant=='Ueff':
                        param/=sc.Boltzmann
                    self.fitted_parameters[quant] = param
            
        # Controls to play with temperature
        self.temp_headline = QLabel('Temperature')
        self.temp_headline.setFont(self.headline_font)
        self.layout.addWidget(self.temp_headline)
        
        self.temp_hbl = QHBoxLayout()
        
        self.temp_min = QDoubleSpinBox()
        self.temp_min.setValue(min_and_max_temps[0])
        self.temp_min.editingFinished.connect(self.temp_interval_changed)
        self.temp_hbl.addWidget(self.temp_min)
        
        self.temp_max = QDoubleSpinBox()
        self.temp_max.setValue(min_and_max_temps[1])
        self.temp_max.editingFinished.connect(self.temp_interval_changed)
        self.temp_hbl.addWidget(self.temp_max)
        
        self.layout.addLayout(self.temp_hbl)
        
        # Controls for which type of plot to consider
        self.plot_headline = QLabel('Plot type to make')
        self.plot_headline.setFont(self.headline_font)
        self.layout.addWidget(self.plot_headline)
        
        self.plot_type_hbl = QHBoxLayout()
        
        self.use_orbach = QCheckBox('Orbach')
        if 'O' in self.plot_type_list: self.use_orbach.setChecked(True)
        self.use_orbach.clicked.connect(self.plot_type_changed)
        self.plot_type_hbl.addWidget(self.use_orbach)
        
        self.use_raman = QCheckBox('Raman')
        if 'R' in self.plot_type_list: self.use_raman.setChecked(True)
        self.use_raman.clicked.connect(self.plot_type_changed)
        self.plot_type_hbl.addWidget(self.use_raman)
        
        self.use_qt = QCheckBox('QT')
        if 'QT' in self.plot_type_list: self.use_qt.setChecked(True)
        self.use_qt.clicked.connect(self.plot_type_changed)
        self.plot_type_hbl.addWidget(self.use_qt)
        
        self.layout.addLayout(self.plot_type_hbl)
        
        self.sim_vals_layout = QFormLayout()
        
        self.tqt_val = QLineEdit()
        self.lineedit_inputs['tQT'] = self.tqt_val
        self.tqt_val.setValidator(self.validator)
        self.tqt_val.setText(str(self.plot_parameters['tQT']))
        self.sim_vals_layout.addRow('t_QT', self.tqt_val)
        
        self.cr_val = QLineEdit()
        self.lineedit_inputs['Cr'] = self.cr_val
        self.cr_val.setValidator(self.validator)
        self.cr_val.setText(str(self.plot_parameters['Cr']))
        self.sim_vals_layout.addRow('C_R', self.cr_val)
        
        self.n_val = QLineEdit()
        self.lineedit_inputs['n'] = self.n_val
        self.n_val.setValidator(self.validator)
        self.n_val.setText(str(self.plot_parameters['n']))
        self.sim_vals_layout.addRow('n', self.n_val)
        
        self.t0_val = QLineEdit()
        self.lineedit_inputs['t0'] = self.t0_val
        self.t0_val.setValidator(self.validator)
        self.t0_val.setText(str(self.plot_parameters['t0']))
        self.sim_vals_layout.addRow('t0', self.t0_val)
        
        self.Ueff_val = QLineEdit()
        self.lineedit_inputs['Ueff'] = self.Ueff_val
        self.Ueff_val.setValidator(self.validator)
        self.Ueff_val.setText(str(self.plot_parameters['Ueff']))
        self.sim_vals_layout.addRow('U_eff', self.Ueff_val)
        
        self.layout.addLayout(self.sim_vals_layout)
        
        # Making control buttons at the end
        self.button_layout = QHBoxLayout()
        
        self.cancel_btn = QPushButton('Cancel')
        self.cancel_btn.setAutoDefault(False)
        self.cancel_btn.clicked.connect(self.reject)
        self.button_layout.addWidget(self.cancel_btn)
        
        self.accept_btn = QPushButton('Ok')
        self.accept_btn.clicked.connect(self.replace_and_accept)
        self.button_layout.addWidget(self.accept_btn)
        
        self.layout.addLayout(self.button_layout)
        
        self.setLayout(self.layout)
        
        self.show()
    
    def fit_take_control(self):
        
        if self.use_fit_cb.isChecked():
            
            for key, val in self.fitted_parameters.items():
                if val is not None:
                    self.lineedit_inputs[key].setText('{}'.format(val))
                self.lineedit_inputs[key].setReadOnly(True)
            
        else:
            self.tqt_val.setReadOnly(False)
            self.cr_val.setReadOnly(False)
            self.n_val.setReadOnly(False)
            self.t0_val.setReadOnly(False)
            self.Ueff_val.setReadOnly(False)
        
        self.param_values_changed()
    
    def param_values_changed(self):
        
        self.plot_parameters['tQT'] = float(self.tqt_val.text())
        self.plot_parameters['Cr'] = float(self.cr_val.text())
        self.plot_parameters['n'] = float(self.n_val.text())
        self.plot_parameters['t0'] = float(self.t0_val.text())
        self.plot_parameters['Ueff'] = float(self.Ueff_val.text())
        
    def plot_type_changed(self):
        
        self.plot_type_list = []
        if self.use_qt.isChecked(): self.plot_type_list.append('QT')
        if self.use_raman.isChecked(): self.plot_type_list.append('R')
        if self.use_orbach.isChecked(): self.plot_type_list.append('O')
        
    def temp_interval_changed(self):
        
        try:
            self.min_and_max_temps[0] = self.temp_min.value()
            self.min_and_max_temps[1] = self.temp_max.value()
            assert self.min_and_max_temps[0]<=self.min_and_max_temps[1]
        except AssertionError:
            pass
            
    def replace_and_accept(self):
        
        self.param_values_changed()
        self.accept()

class AboutDialog(QDialog):
    
    def __init__(self, info):
    
        super(AboutDialog, self).__init__()
        
        self.layout = QVBoxLayout()
        
        self.setWindowTitle('About')
        
        self.author_lbl = QLabel('Written by {}'.format(info['author']))
        self.layout.addWidget(self.author_lbl)
        
        self.web_lbl = QLabel('<a href={}>Molecular magnetism at AU</a>'.format(info['webpage']))
        self.web_lbl.setOpenExternalLinks(True)
        self.layout.addWidget(self.web_lbl)
        
        self.pers_lbl = QLabel('Personal <a href={}>webpage</a>'.format(info['personal']))
        self.pers_lbl.setOpenExternalLinks(True)
        self.layout.addWidget(self.pers_lbl)
        
        self.setLayout(self.layout)
        self.show()

class ParamDialog(QDialog):

    def __init__(self,
                 parent=None,
                 param_dict=None):
                 
        super(ParamDialog, self).__init__()
        
        self.setWindowTitle('Fitted parameters')
        self.dialog_layout = QVBoxLayout()
        self.param_labels = {}
        
        for val in param_dict['quantities']:
            multiplier = 1
            idx = param_dict['quantities'].index(val)
            
            current_label = QLabel()
            if val == 'Ueff': multiplier = sc.Boltzmann
            current_label.setText('{}: {:6.3e} +/- {:6.3e}'.format(val,
                                                                   param_dict['params'][idx]/multiplier,
                                                                   param_dict['sigmas'][idx]/multiplier))
                                                                   
            self.dialog_layout.addWidget(current_label)
        
        self.setLayout(self.dialog_layout)
        self.show()



class FitResultPlotStatus(QDialog):

    def __init__(self, list_input=None):
    
        super(FitResultPlotStatus, self).__init__()
        
        self.layout = QVBoxLayout()
        
        self.scroll = QScrollArea(self)
        self.scroll.setWidgetResizable(True)
        self.layout.addWidget(self.scroll)
        
        self.content = QWidget(self.scroll)
        self.cont_lo = QVBoxLayout(self.content)
        self.content.setLayout(self.cont_lo)
        self.scroll.setWidget(self.content)
        
        self.checked_items = []
        
        num_of_temps = list_input.count()
        for idx in range(num_of_temps):
            item = list_input.item(idx)
            item_lo = QHBoxLayout()
            item_data = item.data(32)
            
            item_fit_bool = item_data['fit']
            item_raw_bool = item_data['raw']
            item_txt = item_data['temp']
            
            raw_checked = QCheckBox('R')
            fit_checked = QCheckBox('F')
            temp = QLabel('{:5.2f}K'.format(item_data['temp']))
            
            item_lo.addWidget(temp)
            item_lo.addWidget(raw_checked)
            item_lo.addWidget(fit_checked)
            
            self.checked_items.append([raw_checked, fit_checked])
            
            raw_checked.setChecked(item_raw_bool)
            fit_checked.setChecked(item_fit_bool)
            
            self.cont_lo.addLayout(item_lo)
        
        self.state_btn_lo = QHBoxLayout()
        
        self.check_all_btn = QPushButton('Check all')
        self.check_all_btn.clicked.connect(self.check_all_function)
        
        self.uncheck_all_btn = QPushButton('Uncheck all')
        self.uncheck_all_btn.clicked.connect(self.uncheck_all_function)
        
        self.state_btn_lo.addWidget(self.uncheck_all_btn)
        self.state_btn_lo.addWidget(self.check_all_btn)
        
        self.layout.addLayout(self.state_btn_lo)
        
        self.judge_btn_lo = QHBoxLayout()
        
        self.states_reject_btn = QPushButton('Cancel')
        self.states_reject_btn.clicked.connect(self.reject)
        self.judge_btn_lo.addWidget(self.states_reject_btn)
        
        self.states_accept_btn = QPushButton('Ok')
        self.states_accept_btn.clicked.connect(self.accept)
        self.judge_btn_lo.addWidget(self.states_accept_btn)
        
        self.layout.addLayout(self.judge_btn_lo)
        
        self.setLayout(self.layout)
        self.show()
        
    def check_all_function(self):
    
        for sublist in self.checked_items:
            sublist[0].setChecked(True)
            sublist[1].setChecked(True)
        
    def uncheck_all_function(self):
    
        for sublist in self.checked_items:
            sublist[0].setChecked(False)
            sublist[1].setChecked(False)