#std packages
from collections import OrderedDict
from multiprocessing.sharedctypes import Value
import os
from re import X 

#third-party packages
import numpy as np 
import scipy.constants as sc

from PyQt5.QtGui import QFont, QDoubleValidator
from PyQt5.QtWidgets import (QInputDialog, QWidget, QPushButton, QLabel, QComboBox, 
                             QDoubleSpinBox, QFormLayout, QCheckBox, QVBoxLayout, QMessageBox,
                             QHBoxLayout, QFileDialog, QDialog, QLineEdit,
                             QScrollArea, QFrame, QGridLayout, QTextEdit)
from PyQt5.QtCore import Qt
from lmfit import fit_report
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT as NavigationToolbar

from .process_ac import default_parameters
from .layout import make_headline, headline_font

#set constants
kB = sc.Boltzmann

class MagMessage(QMessageBox):

    def __init__(self, title, message):

        super(MagMessage, self).__init__()

        self.setWindowTitle(title)
        self.setText(message)

class PlottingWindow(QWidget):

    def __init__(self, make_ax = False):
    
        super(PlottingWindow, self).__init__()
        
        self.layout = QVBoxLayout()
        
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.tools = NavigationToolbar(self.canvas, self)
        if make_ax == "cax":
            self.grid = plt.GridSpec(20,1)
            self.ax = self.fig.add_subplot(self.grid[:17,0])
            self.cax = self.fig.add_subplot(self.grid[-1,0])
            self.cax.set_yticklabels([])
        elif make_ax == "z": 
            self.grid = plt.GridSpec(2,1, height_ratios=[13,1])
            self.ax = self.fig.add_subplot(self.grid[0], projection = '3d')
            self.cax = self.fig.add_subplot(self.grid[1])
            self.cax.set_yticklabels([])
            self.cax.get_yaxis().labelpad = 15
            self.cax.set_ylabel("Temperature (K)")
        else:
            self.ax = self.fig.add_subplot(111)
        
        self.fig.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.95)
        
        self.layout.addWidget(self.canvas)
        
        self.tool_lo = QHBoxLayout()
        self.tool_lo.addWidget(self.tools)
        self.tool_lo.addStretch()
        
        if make_ax != "z": #Reset axes btn does not work for the 3D plot and is therefore not shown in this case
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

        make_headline(self, "Fit history", self.layout)    

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
            lbl.setFont(headline_font)
            layout.addWidget(lbl, idx+1, 0)
        for idx, colname in enumerate(['Value', 'Min.', 'Max.']):
            lbl = QLabel(colname)
            lbl.setFont(headline_font)
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
        self.fit_history_lbl.setFont(headline_font)
        self.layout.addWidget(self.fit_history_lbl)
        
        self.choose_fit_combo = QComboBox()
        self.choose_fit_combo.activated.connect(self.use_fitted_values)
        self.layout.addWidget(self.choose_fit_combo)
        self.update_fit_combo()
        
        # Controls to play with temperature
        self.temp_headline = QLabel('Temperature')
        self.temp_headline.setFont(headline_font)
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
        self.plot_headline.setFont(headline_font)
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
        self.parameter_headline.setFont(headline_font)
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
        
        self.author_lbl = QLabel('Written by: {}'.format(info['author']))
        self.layout.addWidget(self.author_lbl)
        
        self.link_headline = QLabel("Relevant websites:")
        self.layout.addWidget(self.link_headline)

        self.web_lbl = QLabel('<a href={}>Link to research group website </a>'.format(info['webpage']))
        self.web_lbl.setOpenExternalLinks(True)
        self.layout.addWidget(self.web_lbl)
        
        self.pers_lbl = QLabel('<a href={}>Link to github repository with source code</a>'.format(info['personal']))
        self.pers_lbl.setOpenExternalLinks(True)
        self.layout.addWidget(self.pers_lbl)
        
        self.setLayout(self.layout)
        

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
        self.fit_title.setFont(headline_font)
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
        self.setWindowTitle("Temperature subsets to be shown")
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
            
            raw_checked = QCheckBox('Raw data points')
            fit_checked = QCheckBox('Fitted line')
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
        
        self.state_btn_lo.addWidget(self.check_all_btn)

        self.layout.addLayout(self.state_btn_lo)
        

        self.uncheck_lo = QHBoxLayout() 

        self.uncheck_raw_btn = QPushButton('Uncheck all raw data')
        self.uncheck_raw_btn.clicked.connect(self.uncheck_all_raw)

        self.uncheck_fit_btn = QPushButton('Uncheck all fitted lines')
        self.uncheck_fit_btn.clicked.connect(self.uncheck_all_fit)

        self.uncheck_lo.addWidget(self.uncheck_raw_btn)
        self.uncheck_lo.addWidget(self.uncheck_fit_btn)
        
        self.layout.addLayout(self.uncheck_lo)



        self.judge_btn_lo = QHBoxLayout()
        
        self.states_reject_btn = QPushButton('Cancel')
        self.states_reject_btn.clicked.connect(self.reject)
        self.judge_btn_lo.addWidget(self.states_reject_btn)
        
        self.states_accept_btn = QPushButton('Ok')
        self.states_accept_btn.clicked.connect(self.accept)
        self.judge_btn_lo.addWidget(self.states_accept_btn)
        
        self.layout.addLayout(self.judge_btn_lo)
        self.resize(500,700)

        self.setLayout(self.layout)
        #self.show()
        
    def check_all_function(self):
    
        for sublist in self.checked_items:
            sublist[0].setChecked(True)
            sublist[1].setChecked(True)


    def uncheck_all_raw(self): 
        
        for sublist in self.checked_items: 
            sublist[0].setChecked(False)
    
    def uncheck_all_fit(self): 
        
        for sublist in self.checked_items: 
            sublist[1].setChecked(False)

class SampleInformation(QDialog): 
    def __init__(self, parent):
    
        super(SampleInformation, self).__init__()
        self.parent = parent 
        self.layout = QVBoxLayout()
        self.setWindowTitle("Input sample information (used for diamagnetic correction)")

        #Hides Help button
        self.setWindowFlag(Qt.WindowContextHelpButtonHint,False)  

        ## SAMPLE INFO
        self.sample_info_layout = QVBoxLayout()
        self.sample_info_header = QLabel('Sample information')
        self.sample_info_header.setFont(headline_font)
        self.sample_info_layout.addWidget(self.sample_info_header)
        
        ## Sample mass edit
        self.sample_mass_layout = QHBoxLayout()
        self.sample_info_layout.addLayout(self.sample_mass_layout)
        
        self.sample_mass_lbl = QLabel('m (sample) [mg]')
        self.sample_mass_layout.addWidget(self.sample_mass_lbl)
        
        self.sample_mass_inp = QLineEdit()
        self.sample_mass_inp.setValidator(QDoubleValidator())
        self.sample_mass_layout.addWidget(self.sample_mass_inp)

        ## Sample molar mass edit
        self.molar_mass_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.molar_mass_lo)
        
        self.molar_mass_lbl = QLabel('M (sample) [g/mol]')
        self.molar_mass_lo.addWidget(self.molar_mass_lbl)
        
        self.molar_mass_inp = QLineEdit()
        self.molar_mass_inp.setValidator(QDoubleValidator())
        self.molar_mass_lo.addWidget(self.molar_mass_inp)


        ## Sample Xd edit
        self.sample_xd_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.sample_xd_lo)
        

        self.sample_xd_lbl = QLabel(u"<a href={}>X\u1D05</a>".format(
                                    self.parent.tooltips_dict['English']['Xd_link'])
                                    +' (sample) [emu/(Oe*mol)]')
        self.sample_xd_lbl.setOpenExternalLinks(True)
        self.sample_xd_lo.addWidget(self.sample_xd_lbl)
        
        self.sample_xd_inp = QLineEdit()
        self.sample_xd_inp.setValidator(QDoubleValidator())
        self.sample_xd_lo.addWidget(self.sample_xd_inp)
        
        # Constant terms edit
        self.constant_terms_layout = QHBoxLayout()
        self.sample_info_layout.addLayout(self.constant_terms_layout)
        
        self.constant_terms_lbl = QLabel('Constant terms')
        self.constant_terms_layout.addWidget(self.constant_terms_lbl)
             
        self.constant_terms_inp = QLineEdit()
        self.constant_terms_layout.addWidget(self.constant_terms_inp)
        
        # Variable amount edit
        self.var_amount_layout = QHBoxLayout()
        self.sample_info_layout.addLayout(self.var_amount_layout)
        
        self.var_amount_lbl = QLabel('Variable amounts')
        self.var_amount_layout.addWidget(self.var_amount_lbl)
        
        self.var_amount_inp = QLineEdit()
        self.var_amount_layout.addWidget(self.var_amount_inp)
        
        #Insert values to the QLineEdits if known: 
        self.insert_values() 

        # Mass load button
        self.sample_data_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.sample_data_lo)
        
        self.load_sample_data_btn = QPushButton('Load sample data from file')
        self.load_sample_data_btn.clicked.connect(self.load_sample_data)
        self.sample_data_lo.addWidget(self.load_sample_data_btn)
        
        self.save_sample_data_btn = QPushButton('Save sample data in file')
        self.save_sample_data_btn.clicked.connect(self.save_sample_data)
        self.sample_data_lo.addWidget(self.save_sample_data_btn)
        
        self.done_sample_data_btn = QPushButton('Ok')
        self.done_sample_data_btn.clicked.connect(self.save_sample_data_for_diamag)
        self.sample_data_lo.addWidget(self.done_sample_data_btn)

        #Sets gui langugage for tooltips
        self.set_gui_language() 

        #Finalizing layout
        self.sample_data_lo.addStretch()
        self.setLayout(self.sample_info_layout)


    def load_sample_data(self):

        filename_info = QFileDialog().getOpenFileName(self, 'Open file', self.parent.last_loaded_file)
        filename = filename_info[0]
        try:
            f = open(filename, 'r')
            d = f.readlines()
            f.close()
            
            assert all([len(line.split())>=2 for line in d])
            
        except FileNotFoundError:
            MagMessage("Error", "File was not selected").exec_() 
        except UnicodeDecodeError:
            MagMessage("Error", "Cannot open binary file").exec_() 
        except AssertionError:
            MagMessage("Error", "Some lines have length less than two").exec_() 
        else:
            
            # These are the default values that are "read" if nothing else is
            # seen in the file
            m_sample = '0'
            M_sample = '0'
            Xd_sample = '0'
            constant_terms = '0'
            var_amount = '0,0'
            
            self.parent.last_loaded_file = os.path.split(filename)[0]
            for line in d:
                line = line.split()
                if line[0] == 'm_sample':
                    m_sample = line[1]
                elif line[0] == 'M_sample':
                    M_sample = line[1]
                elif line[0] == 'Xd_sample':
                    Xd_sample = line[1]
                elif line[0] == 'constants':
                    constant_terms = line[1]
                elif line[0] == 'var_amount':
                    var_amount = line[1]
            
            self.sample_mass_inp.setText(m_sample)
            self.molar_mass_inp.setText(M_sample)
            self.sample_xd_inp.setText(Xd_sample)
            self.constant_terms_inp.setText(constant_terms)
            self.var_amount_inp.setText(var_amount)
    
    def save_sample_data(self):
    
        filename_info = QFileDialog.getSaveFileName(self,
                                                   'Save sample file',
                                                   self.parent.last_loaded_file)
        filename = filename_info[0]
        
        try:
            assert filename != ''
            self.parent.last_loaded_file = os.path.split(filename)[0]
            filename, ext = os.path.splitext(filename)
            if ext == '':
                ext = '.dat'
            
            comment, ok = QInputDialog.getText(self,
                                              'Comment',
                                              'Comment for saved sample data')

            fc = ''
            fc += '# ' + comment + '\n'
            fc += 'm_sample ' + self.sample_mass_inp.text() + '\n'
            fc += 'M_sample ' + self.molar_mass_inp.text() + '\n'
            fc += 'Xd_sample ' + self.sample_xd_inp.text() + '\n'
            fc += 'constants ' + self.constant_terms_inp.text() + '\n'
            fc += 'var_amount ' + self.var_amount_inp.text() + '\n'
            
            f = open(filename+ext, 'w')
            f.write(fc)
            f.close()
            
        except AssertionError:
            pass
            
    def save_sample_data_for_diamag(self): 
        self.parent.data_treat.m_sample = self.sample_mass_inp.text()
        self.parent.data_treat.M_sample = self.molar_mass_inp.text()
        self.parent.data_treat.Xd_sample = self.sample_xd_inp.text()
        self.parent.data_treat.constant_terms = self.constant_terms_inp.text()
        self.parent.data_treat.var_am = self.var_amount_inp.text()
        self.accept()
    
    def insert_values(self): 
        try: 
            self.sample_mass_inp.setText(str(self.parent.data_treat.m_sample))
        except AttributeError: 
            pass
        try: 
            self.molar_mass_inp.setText(str(self.parent.data_treat.M_sample))
        except AttributeError: 
            pass
        try: 
            self.sample_xd_inp.setText(str(self.parent.data_treat.Xd_sample))
        except AttributeError: 
            pass
        try: 
            self.constant_terms_inp.setText(str(self.parent.data_treat.constant_terms))
        except AttributeError: 
            pass
        try: 
            self.var_amount_inp.setText(str(self.parent.data_treat.var_am))
        except AttributeError: 
            pass
    
    def set_gui_language(self): 
        language = self.parent.gui_language 
        try: 
            self.sample_mass_inp.setToolTip(self.parent.tooltips_dict[language]['m_sample'])
            self.molar_mass_inp.setToolTip(self.parent.tooltips_dict[language]['M_sample'])
            self.sample_xd_inp.setToolTip(self.parent.tooltips_dict[language]['Xd_sample'])
            self.constant_terms_inp.setToolTip(self.parent.tooltips_dict[language]['const_terms'])
            self.var_amount_inp.setToolTip(self.parent.tooltips_dict[language]['var_amounts'])
        except: 
            pass 




