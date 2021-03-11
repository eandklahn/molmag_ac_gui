#std packages
from collections import OrderedDict

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
            self.cax.set_yticklabels([])
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
                 guess=None,
                 fit_history=None):
        super(GuessDialog, self).__init__()
        
        self.layout = QVBoxLayout()
        self.setWindowTitle('Guess parameters')

        self.validator = QDoubleValidator()
        self.validator.setNotation(QDoubleValidator.ScientificNotation)
        
        self.fit_history = fit_history
        self.init_guess = guess
        self.values = []
        
        self.parameter_inputs = OrderedDict()
        self.parameter_inputs['tQT']=None
        self.parameter_inputs['Cr']=None
        self.parameter_inputs['n']=None
        self.parameter_inputs['t0']=None
        self.parameter_inputs['Ueff']=None
        
        self.fit_history_lbl = QLabel('Fit history (latest first)')
        self.layout.addWidget(self.fit_history_lbl)
        
        self.fit_history_combo = QComboBox()
        for e in self.fit_history:
            rep = self.fit_history_element_repr(e)
            self.fit_history_combo.addItem(rep)
        self.layout.addWidget(self.fit_history_combo)
        self.fit_history_combo.activated.connect(self.fit_take_control)
        
        self.sim_vals_layout = QFormLayout()
        for key in self.parameter_inputs.keys():
            self.parameter_inputs[key] = QLineEdit()
            self.parameter_inputs[key].setValidator(self.validator)
            if key=='Ueff':
                self.parameter_inputs[key].setText(str(self.init_guess[key]/kB))
            else:
                self.parameter_inputs[key].setText(str(self.init_guess[key]))
            self.sim_vals_layout.addRow(key, self.parameter_inputs[key])
        self.layout.addLayout(self.sim_vals_layout)
        
        accept_btn = QPushButton('Fit')
        accept_btn.clicked.connect(self.on_close)
        self.layout.addWidget(accept_btn)

        self.setLayout(self.layout)
        self.show()
    
    def fit_history_element_repr(self, e):
        
        fit_type = e[0]
        fit_dict = e[1]
        params = fit_dict['params']
        quants = fit_dict['quantities']
        
        rep = []
        for key in self.parameter_inputs.keys():
            if key in quants:
                idx = quants.index(key)
                param_val = params[idx]
                if key=='Ueff':
                    param_val /= kB
                rep.append(f'{key}: {param_val:.2e}')
            else:
                rep.append(f'{key}: None')
                
        rep = f'Fit type: {fit_type}'+'\n'+'\n'.join(rep)    
        
        return rep    
    
    def fit_take_control(self):
        
        idx = self.fit_history_combo.currentIndex()
        fit = self.fit_history[idx]
        
        fit_dict = fit[1]
        params = fit_dict['params']
        quants = fit_dict['quantities']
        
        for key, val in self.parameter_inputs.items():
            if key in quants:
                key_idx = quants.index(key)
                new_val = params[key_idx]
                if key == 'Ueff':
                    new_val /= kB
                val.setText(str(new_val))
                
    def on_close(self):
        
        self.return_guess = {key: float(val.text()) for key, val in self.parameter_inputs.items()}
        self.return_guess['Ueff'] = self.return_guess['Ueff']*kB
        self.accept()

class SimulationDialog(QDialog):

    def __init__(self,
                 parent=None,
                 fit_history=[],
                 plot_type_list=[],
                 plot_parameters={'tQT': 0.1, 'Cr': 0.1, 'n': 0.1, 't0': 0.1, 'Ueff': 0.1},
                 min_and_max_temps=[0]*2):
    
        super(SimulationDialog, self).__init__()
        
        self.setWindowTitle('Add simulation')
        self.headline_font = QFont()
        self.headline_font.setBold(True)
        
        # Storing input values
        self.plot_type_list = plot_type_list
        self.plot_parameters = plot_parameters
        self.min_and_max_temps = min_and_max_temps
        self.fit_history = fit_history

        # Abstracting the validator for the QLineEdits
        self.validator = QDoubleValidator()
        self.validator.setNotation(QDoubleValidator.ScientificNotation)
        
        # Containers for objects
        self.parameter_inputs = OrderedDict()
        self.parameter_inputs['tQT']=None
        self.parameter_inputs['Cr']=None
        self.parameter_inputs['n']=None
        self.parameter_inputs['t0']=None
        self.parameter_inputs['Ueff']=None
        
        self.possible_functions = OrderedDict()
        self.possible_functions['QT']=['QT',None]
        self.possible_functions['R']=['Raman',None]
        self.possible_functions['O']=['Orbach',None]
        
        self.layout = QVBoxLayout()
        
        self.fit_history_lbl = QLabel('Fit history (latest first)')
        self.fit_history_lbl.setFont(self.headline_font)
        self.layout.addWidget(self.fit_history_lbl)
        
        self.fit_history_combo = QComboBox()
        for fit in self.fit_history:
            rep = self.fit_history_element_repr(fit)
            self.fit_history_combo.addItem(rep)
        
        self.fit_history_combo.activated.connect(self.fit_take_control)
        self.layout.addWidget(self.fit_history_combo)
        
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
        
        self.temp_hbl.addStretch()
        
        self.layout.addLayout(self.temp_hbl)
        
        # Controls for which type of plot to consider
        self.plot_headline = QLabel('Plot type to make')
        self.plot_headline.setFont(self.headline_font)
        self.layout.addWidget(self.plot_headline)
        
        self.plot_type_hbl = QHBoxLayout()
        for key, val in reversed(self.possible_functions.items()):
            shorthand = key
            fullname = val[0]
            val[1] = QCheckBox(fullname)
            val[1].clicked.connect(self.plot_type_changed)
            if key in self.plot_type_list: val[1].setChecked(True)
            self.plot_type_hbl.addWidget(val[1])
            
        self.plot_type_hbl.addStretch()
        self.layout.addLayout(self.plot_type_hbl)
        
        # Values to use
        self.sim_vals_layout = QFormLayout()
        for key in self.parameter_inputs.keys():
            self.parameter_inputs[key] = QLineEdit()
            self.parameter_inputs[key].setValidator(self.validator)
            self.parameter_inputs[key].setText(str(self.plot_parameters[key]))
            self.sim_vals_layout.addRow(key, self.parameter_inputs[key])
        self.layout.addLayout(self.sim_vals_layout)
        
        # Making control buttons at the end
        self.button_layout = QHBoxLayout()
        
        self.cancel_btn = QPushButton('Cancel')
        self.cancel_btn.setAutoDefault(False)
        self.cancel_btn.clicked.connect(self.reject)
        self.button_layout.addWidget(self.cancel_btn)
        
        self.accept_btn = QPushButton('Ok')
        self.accept_btn.setAutoDefault(True)
        self.accept_btn.clicked.connect(self.replace_and_accept)
        self.button_layout.addWidget(self.accept_btn)
        
        self.layout.addLayout(self.button_layout)
        
        self.setLayout(self.layout)
        
        self.show()
    
    def fit_history_element_repr(self, e):
        
        fit_type = e[0]
        fit_dict = e[1]
        params = fit_dict['params']
        quants = fit_dict['quantities']
        
        rep = []
        for key in self.parameter_inputs.keys():
            if key in quants:
                idx = quants.index(key)
                param_val = params[idx]
                if key=='Ueff':
                    param_val /= kB
                rep.append(f'{key}: {param_val:.2e}')
            else:
                rep.append(f'{key}: None')
                
        rep = f'Fit type: {fit_type}'+'\n'+'\n'.join(rep)    
        
        return rep
    
    def fit_take_control(self):
        
        idx = self.fit_history_combo.currentIndex()
        fit = self.fit_history[idx]
        
        fit_dict = fit[1]
        params = fit_dict['params']
        quants = fit_dict['quantities']
        
        for key, val in self.parameter_inputs.items():
            if key in quants:
                key_idx = quants.index(key)
                new_val = params[key_idx]
                if key == 'Ueff':
                    new_val /= kB
                val.setText(str(new_val))
            
    def param_values_changed(self):

        for key, val in self.parameter_inputs.items():
            self.plot_parameters[key] = float(val.text())
    
    def plot_type_changed(self):

        self.plot_type_list = []
        for key in self.possible_functions.keys():
            val = self.possible_functions[key]
            if val[1].isChecked(): self.plot_type_list.append(key)

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
                 fit_history,
                 parent=None):
                 
        super(ParamDialog, self).__init__()
        
        self.setWindowTitle('Fitted parameters')
        self.layout = QVBoxLayout()

        self.fit_history = fit_history
        
        self.parameter_labels = OrderedDict()
        self.parameter_labels['tQT']=QLabel()
        self.parameter_labels['Cr']=QLabel()
        self.parameter_labels['n']=QLabel()
        self.parameter_labels['t0']=QLabel()
        self.parameter_labels['Ueff']=QLabel()
        
        self.fit_history_combo = QComboBox()
        for e in self.fit_history:
            rep = self.fit_history_element_repr(e)
            self.fit_history_combo.addItem(rep)
        self.fit_history_combo.activated.connect(self.show_fit)
        self.layout.addWidget(self.fit_history_combo)
        
        for key, val in self.parameter_labels.items():
            self.layout.addWidget(val)        
            
        self.fit_history_combo.setCurrentIndex(0)
        self.show_fit()
        self.setLayout(self.layout)
        self.show()

    def show_fit(self):
        
        fit_idx = self.fit_history_combo.currentIndex()
        fit = self.fit_history[fit_idx]
        
        quants = fit[1]['quantities']
        params = fit[1]['params']
        sigmas = fit[1]['sigmas']
        
        for key, val in self.parameter_labels.items():
            if key in quants:
                key_idx = quants.index(key)
                key_param = params[key_idx]
                key_sigma = sigmas[key_idx]
                val.setText(f'{key} = {key_param:.6e} +- {key_sigma:.6e}')
            else:
                val.setText(f'{key} = None')
        
    def fit_history_element_repr(self, e):
        
        fit_type = e[0]
        fit_dict = e[1]
        params = fit_dict['params']
        quants = fit_dict['quantities']
        
        rep = []
        for key in self.parameter_labels.keys():
            if key in quants:
                idx = quants.index(key)
                param_val = params[idx]
                if key=='Ueff':
                    param_val /= kB
                rep.append(f'{key}: {param_val:.2e}')
            else:
                rep.append(f'{key}: None')
                
        rep = f'Fit type: {fit_type}'+'\n'+'\n'.join(rep)    
        
        return rep
    


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