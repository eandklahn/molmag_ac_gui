#std packages
from collections import OrderedDict

#third-party packages
import scipy.constants as sc

from PyQt5.QtGui import QIcon, QFont, QDoubleValidator
from PyQt5.QtWidgets import (QFrame, QMainWindow, QWidget, QApplication, QPushButton, QLabel, QAction, QComboBox, QStackedWidget,
                             QDoubleSpinBox, QFormLayout, QCheckBox, QVBoxLayout, QMessageBox, QSplitter, QGridLayout,
                             QHBoxLayout, QFileDialog, QDialog, QLineEdit, QListWidget, QListWidgetItem, QTabWidget,
                             QScrollArea, QStatusBar)

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT as NavigationToolbar

from .process_ac import default_parameters

#set constants
kB = sc.Boltzmann

class MagMessage(QMessageBox):

    def __init__(self, title, message):

        super(MagMessage, self).__init__()

        self.setWindowTitle(title)
        self.setText(message)

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
           
            lines_to_manage = []
            for line in self.ax.lines:
                if len(line.get_xdata())<1: pass
                elif not line._visible: pass
                else: lines_to_manage.append(line)
                
            x = lines_to_manage[0].get_xdata()
            y = lines_to_manage[0].get_ydata()
            
            new_x = [x.min(), x.max()]
            new_y = [y.min(), y.max()]
            
            for line in lines_to_manage:
                x = line.get_xdata()
                y = line.get_ydata()
                
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
                 parent,
                 guess,
                 fitwith):
        
        super(GuessDialog, self).__init__()
        
        self.setWindowTitle('Initial fit parameters')
        self.validator = QDoubleValidator()
        self.validator.setNotation(QDoubleValidator.ScientificNotation)
        
        self.parent = parent
        self.fit_history = parent.fit_history
        self.fitwith = fitwith
        self.init_guess = guess
        
        self.init_params = default_parameters(self.fitwith)
        self.set_init_params()

        self.final_params = default_parameters(self.fitwith)
        
        self.current_choice = self.init_params

        self.initUI()

    def initUI(self):

        self.layout = QVBoxLayout()

        self.fit_history_lbl = QLabel('Fit history')
        self.fit_history_lbl.setFont(self.parent.headline_font)
        self.layout.addWidget(self.fit_history_lbl)

        self.choose_fit_combo = QComboBox()
        self.update_fit_combo()
        self.choose_fit_combo.activated.connect(self.use_fit)
        self.layout.addWidget(self.choose_fit_combo)
        
        self.make_value_frame()
        self.set_params()

        self.accept_btn = QPushButton('Accept')
        self.accept_btn.clicked.connect(self.onclose)
        self.layout.addWidget(self.accept_btn)

        self.setLayout(self.layout)
    
    def make_value_frame(self):

        w = QFrame()
        layout = QHBoxLayout()
        w.setStyleSheet('background-color: rgb(240,240,240)')
        
        name = ['0']
        lbls = [QLabel('Name')]
        vals = [QLabel('Value')]
        mins = [QLabel('Min.')]
        maxs = [QLabel('Max.')]

        for p in [p for p in self.init_params if not 'use' in p]:
            newlbl = QLabel(p)
            newval = QLineEdit()
            newmin = QLineEdit()
            newmax = QLineEdit()

            newval.setValidator(self.validator)
            newmin.setValidator(self.validator)
            newmax.setValidator(self.validator)

            lbls.append(newlbl)
            vals.append(newval)
            mins.append(newmin)
            maxs.append(newmax)

        columns = [lbls, vals, mins, maxs]
        for items in columns:
            lo = QVBoxLayout()
            for item in items:
                lo.addWidget(item)
            layout.addLayout(lo)

        w.setLayout(layout)
        self.layout.addWidget(w)

        self.value_frame = {'lbls': lbls,
                            'vals': vals,
                            'mins': mins,
                            'maxs': maxs}

    def set_init_params(self):
        
        for key, val in self.init_guess.items():
            self.init_params[key].set(value=val)

    def set_params(self):
        
        params = self.current_choice

        for i, p in enumerate(self.value_frame['lbls'][1:]):
            i += 1
            name = p.text()
            param = params[name]
            vary = param.vary
            val = param.value
            min_val = param.min
            max_val = param.max

            if vary:
                self.value_frame['vals'][i].setText(str(val))
                self.value_frame['mins'][i].setText(str(min_val))
                self.value_frame['maxs'][i].setText(str(max_val))

    def use_fit(self):
        
        idx = self.choose_fit_combo.currentIndex()
        name, res, time = self.fit_history[idx]
        self.current_choice = res.params
        self.set_params()

    def update_fit_combo(self):

        for fit in self.fit_history:
            name, res, time = fit
            params = res.params
            repr = f'{time}: {name}'
            L = [p for p in params if not ('use' in p or params[p].vary==False)]
            for p in L:
                param = params[p]
                repr += f'\n{param.name}: {param.value}'
            self.choose_fit_combo.addItem(repr)
    
    def set_final_params(self):

        names = [lbl.text() for lbl in self.value_frame['lbls']]
        
        for p in self.final_params:
            param = self.final_params[p]
            if param.vary:
                idx = names.index(p)
                val = float(self.value_frame['vals'][idx].text())
                val_min = float(self.value_frame['mins'][idx].text())
                val_max = float(self.value_frame['maxs'][idx].text())
                param.set(value=val,
                          min=val_min,
                          max=val_max)

    def onclose(self):

        self.set_final_params()
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
        
        #self.show()
    
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
        #self.show()

class ParamDialog(QDialog):

    def __init__(self,
                 parent,
                 fit_history):
                 
        super(ParamDialog, self).__init__()
        
        self.parent = parent
        self.fit_history = fit_history
        
        self.setWindowTitle('Fitted parameters')
        
        self.initUI()
        self.update_fit_combo()
        
    def initUI(self):
        
        self.layout = QVBoxLayout()

        self.choose_fit_combo = QComboBox()
        self.choose_fit_combo.currentIndexChanged.connect(
                                                  self.show_MinimizerResult)
        self.layout.addWidget(self.choose_fit_combo)

        self.fit_title = QLabel()
        self.fit_title.setFont(self.parent.headline_font)
        self.layout.addWidget(self.fit_title)

        self.setLayout(self.layout)

    def update_fit_combo(self):

        for fit in self.fit_history:
            name, res, time = fit
            self.choose_fit_combo.addItem(f'{time}: {name}')

    def show_MinimizerResult(self):
        
        fit_idx = self.choose_fit_combo.currentIndex()
        name, res, time = self.fit_history[fit_idx]
        title = f'{time}: {name}'

        self.fit_title.setText(title)

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
        #self.show()
        
    def check_all_function(self):
    
        for sublist in self.checked_items:
            sublist[0].setChecked(True)
            sublist[1].setChecked(True)
        
    def uncheck_all_function(self):
    
        for sublist in self.checked_items:
            sublist[0].setChecked(False)
            sublist[1].setChecked(False)