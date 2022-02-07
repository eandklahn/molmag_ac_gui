#std packages
from collections import OrderedDict

#third-party packages
import scipy.constants as sc

from PyQt5.QtGui import QIcon, QFont, QDoubleValidator
from PyQt5.QtWidgets import (QFrame, QMainWindow, QTextEdit, QWidget, QApplication, QPushButton, QLabel, QAction, QComboBox, QStackedWidget,
                             QDoubleSpinBox, QFormLayout, QCheckBox, QVBoxLayout, QMessageBox, QSplitter, QGridLayout,
                             QHBoxLayout, QFileDialog, QDialog, QLineEdit, QListWidget, QListWidgetItem, QTabWidget,
                             QScrollArea, QStatusBar, QGridLayout)

from lmfit import fit_report

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
        self.current_guess = guess
        self.param_names = [p for p in self.current_guess if not 'use' in p]
        self.valueedits = None
        
        self.initUI()

    def initUI(self):

        self.layout = QVBoxLayout()

        self.fit_history_lbl = QLabel('Fit history')
        self.fit_history_lbl.setFont(self.parent.headline_font)
        self.layout.addWidget(self.fit_history_lbl)

        self.choose_fit_combo = QComboBox()
        self.update_fit_combo()
        self.choose_fit_combo.activated.connect(self.get_values_from_fit)
        self.layout.addWidget(self.choose_fit_combo)
        
        self.make_value_frame()
        self.write_params()

        self.accept_btn = QPushButton('Accept')
        self.accept_btn.clicked.connect(self.onclose)
        self.layout.addWidget(self.accept_btn)

        self.setLayout(self.layout)
    
    def make_value_frame(self):

        w = QFrame()
        layout = QGridLayout()
        w.setStyleSheet('background-color: rgb(240,240,240)')
        
        for idx, rowname in enumerate(self.param_names):
            lbl = QLabel(rowname)
            lbl.setFont(self.parent.headline_font)
            layout.addWidget(lbl, idx+1, 0)
        for idx, colname in enumerate(['Value', 'Min.', 'Max.']):
            lbl = QLabel(colname)
            lbl.setFont(self.parent.headline_font)
            layout.addWidget(lbl, 0, idx+1)

        self.valueedits = [None]
        for row in range(1,6):
            self.valueedits.append([None])
            for col in range(1,4):
                ledit = QLineEdit()
                ledit.setValidator(self.validator)
                self.valueedits[-1].append(ledit)
                layout.addWidget(ledit, row, col)

        self.value_frame = w
        w.setLayout(layout)
        self.layout.addWidget(w)

    def write_params(self):
        
        params = self.current_guess

        for idx, name in enumerate(self.param_names):
            idx += 1
            param = params[name]
            if self.current_guess[name].vary:
                self.valueedits[idx][1].setText(str(param.value))
                self.valueedits[idx][2].setText(str(param.min))
                self.valueedits[idx][3].setText(str(param.max))
                
    def read_params(self):

        for p in self.current_guess:
            param = self.current_guess[p]
            if param.vary:
                idx = self.param_names.index(p)+1
                set_val = float(self.valueedits[idx][1].text())
                val_min = float(self.valueedits[idx][2].text())
                val_max = float(self.valueedits[idx][3].text())
                param.set(value=set_val,
                          min=val_min,
                          max=val_max)

    def get_values_from_fit(self):
        
        idx = self.choose_fit_combo.currentIndex()
        name, res, time = self.fit_history[idx]
        potential_guess = res.params
        for p in potential_guess:
            current_param = self.current_guess[p]
            potential_param = potential_guess[p]
            if (current_param.vary and potential_param.vary):
                current_param.value = potential_param.value
        self.write_params()

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
    
    def onclose(self):

        self.read_params()
        self.accept()

class SimulationDialog(QDialog):

    def __init__(self,
                 parent=None,
                 fit_history=[],
                 params=default_parameters(),
                 min_max_T=[0,0]):

        super(SimulationDialog, self).__init__()

        self.setWindowTitle('Add/edit simulation')
        
        self.parent = parent
        self.fit_history = fit_history
        self.params = params
        self.min_max_T = min_max_T
        
        self.use_function_checkboxes = {}
        self.use_values_edits = {}

        self.validator = QDoubleValidator()
        self.validator.setNotation(QDoubleValidator.ScientificNotation)
        
        self.initUI()

    def initUI(self):

        self.layout = QVBoxLayout()
        
        self.fit_history_lbl = QLabel('Fit history')
        self.fit_history_lbl.setFont(self.parent.headline_font)
        self.layout.addWidget(self.fit_history_lbl)
        
        self.choose_fit_combo = QComboBox()
        self.choose_fit_combo.activated.connect(self.use_fitted_values)
        self.layout.addWidget(self.choose_fit_combo)
        self.update_fit_combo()
        
        # Controls to play with temperature
        self.temp_headline = QLabel('Temperature')
        self.temp_headline.setFont(self.parent.headline_font)
        self.layout.addWidget(self.temp_headline)
        
        self.temp_hbl = QHBoxLayout()
        
        self.temp_min = QDoubleSpinBox()
        self.temp_min.setValue(self.min_max_T[0])
        self.temp_hbl.addWidget(self.temp_min)
        
        self.temp_max = QDoubleSpinBox()
        self.temp_max.setValue(self.min_max_T[1])
        self.temp_hbl.addWidget(self.temp_max)
        
        self.temp_hbl.addStretch()
        self.layout.addLayout(self.temp_hbl)
        
        # Controls for which type of plot to consider
        self.plot_headline = QLabel('Plot type')
        self.plot_headline.setFont(self.parent.headline_font)
        self.layout.addWidget(self.plot_headline)
        
        self.plot_type_hbl = QHBoxLayout()
        functions = [p for p in self.params if 'use' in p]
        for p in functions:
            name = p[3:]
            use = bool(self.params[p].value)
            self.use_function_checkboxes[p] = QCheckBox(name)
            self.use_function_checkboxes[p].setChecked(use)
            self.plot_type_hbl.addWidget(self.use_function_checkboxes[p])

        self.plot_type_hbl.addStretch()
        self.layout.addLayout(self.plot_type_hbl)
        
        # Values to use
        self.parameter_headline = QLabel('Parameter values')
        self.parameter_headline.setFont(self.parent.headline_font)
        self.layout.addWidget(self.parameter_headline)

        self.sim_vals_layout = QFormLayout()
        params = [p for p in self.params if not 'use' in p]
        for p in params:
            self.use_values_edits[p] = QLineEdit()
            self.use_values_edits[p].setValidator(self.validator)
            self.use_values_edits[p].setText(str(self.params[p].value))
            self.sim_vals_layout.addRow(p, self.use_values_edits[p])
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

    def use_fitted_values(self):
        
        idx = self.choose_fit_combo.currentIndex()
        name, res, time = self.fit_history[idx]
        
        new_params = res.params
        param_names = [p for p in new_params if not 'use' in p]
        for p in param_names:
            self.use_values_edits[p].setText(str(new_params[p].value))

    def read_param_values(self):
        
        param_names = [p for p in self.params if not 'use' in p]
        for p in param_names:
            textval = self.use_values_edits[p].text()
            self.params[p].value = float(textval)
    
    def read_plot_type(self):
        
        function_names = [p for p in self.params if 'use' in p]
        for p in function_names:
            self.params[p].value = int(self.use_function_checkboxes[p].isChecked())

    def check_temperature(self):
        
        self.min_max_T[0] = self.temp_min.value()
        self.min_max_T[1] = self.temp_max.value()
        try:
            assert self.min_max_T[0]>0
            assert self.min_max_T[0]<self.min_max_T[1]
        except AssertionError:
            w = MagMessage("Wrong temperatures",
                           "The minimum must be larger than 0 and lower than the maximum")
            w.exec_()
            return False
        else:
            return True
            
    def replace_and_accept(self):
        
        try:
            assert self.check_temperature()
        except AssertionError:
            pass
        else:
            self.read_plot_type()
            self.read_param_values()
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

        self.fit_summary = QTextEdit()
        self.layout.addWidget(self.fit_summary)

        self.setLayout(self.layout)

    def update_fit_combo(self):

        for fit in self.fit_history:
            name, res, time = fit
            self.choose_fit_combo.addItem(f'{time}: {name}')

    def show_MinimizerResult(self):
        
        fit_idx = self.choose_fit_combo.currentIndex()
        name, res, time = self.fit_history[fit_idx]
        title = f'{time}: {name}'

        self.fit_summary.setText(fit_report(res))

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