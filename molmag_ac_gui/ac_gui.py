#std packages
import ctypes
import sys
import re
import os
import time
from subprocess import Popen, PIPE

#third-party packages
import numpy as np
import pandas as pd

from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
#import matplotlib.pyplot as plt
#from matplotlib.figure import Figure
#from matplotlib.backends.qt_compat import QtCore, QtWidgets
#from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT as NavigationToolbar
#if int(QtCore.qVersion().split('.')[0])>=5:
#    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

import scipy.constants as sc
from scipy.optimize import minimize, curve_fit
from lmfit import Parameters, minimize

from PyQt5.QtWinExtras import QWinTaskbarButton
from PyQt5.QtGui import QIcon, QFont, QDoubleValidator
from PyQt5.QtWidgets import (QMainWindow, QWidget, QApplication, QPushButton, QLabel, QAction, QComboBox, QStackedWidget,
                             QDoubleSpinBox, QFormLayout, QCheckBox, QVBoxLayout, QMessageBox, QSplitter, QGridLayout,
                             QHBoxLayout, QFileDialog, QDialog, QLineEdit, QListWidget, QListWidgetItem, QTabWidget,
                             QScrollArea, QStatusBar)

#local imports
from .process_ac import (Xp_, Xpp_, Xp_dataset, Xpp_dataset, getParameterGuesses,
                         getStartParams, getFittingFunction, readPopt, addPartialModel, tau_err_RC)
from .dialogs import GuessDialog, SimulationDialog, AboutDialog, ParamDialog, FitResultPlotStatus, PlottingWindow
from .utility import read_ppms_file, get_ppms_column_name_matches, update_data_names
from .exceptions import FileFormatError

#set constants
kB = sc.Boltzmann

"""
MAIN GUI WINDOW
"""

class ACGui(QMainWindow):

    def __init__(self):
    
        super().__init__()
        
        self.initUI()
        
    def initUI(self):
        
        """About"""
        self.about_information = {'author': 'Emil A. Klahn (eklahn@chem.au.dk)',
                                  'webpage': 'https://chem.au.dk/en/molmag',
                                  'personal': 'https://eandklahn.github.io'}
        
        self.startUp = True
        
        """ Things to do with how the window is shown """
        self.setWindowTitle('AC Processing')
        self.setWindowIcon(QIcon('double_well_potential_R6p_icon.ico'))
        
        """ FONT STUFF """
        self.headline_font = QFont()
        self.headline_font.setBold(True)
        
        """ Setting up the main tab widget """
        self.all_the_tabs = QTabWidget()
        self.setCentralWidget(self.all_the_tabs)
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        
        """ Constructing the data analysis tab ---------------------------------------------------------------------------------- """
        # Data containers for analysis
        self.last_loaded_file = os.getcwd()
        
        self.data_file_name = None
        self.data_file_dir = os.getcwd()
        
        self.data_T = None
        self.data_tau = None
        self.data_dtau = None
        
        self.plotted_data_pointers = None
        self.data_used_pointer = None
        self.data_not_used_pointer = None
        
        self.used_T = None
        self.not_used_T = None
        
        self.used_tau = None
        self.not_used_tau = None
        
        self.used_dtau = None
        self.not_used_dtau = None
        
        self.fitted_parameters = None
        self.used_indices = None
        
        # Data containers for treatment
        
        self.read_options = {
                             'AC Frequency (Hz)': [
                                'Frequency (Hz)', 'AC Frequency (Hz)'],
                             'AC Amplitude (Oe)': [
                                'Amplitude (Oe)', 'AC Drive (Oe)'],
                             'Magnetic Field (Oe)': [
                                'Magnetic Field (Oe)'],
                             'Mp (emu)': [
                                "M' (emu)"],
                             'Mpp (emu)': [
                                "M'' (emu)"],
                             'Xp (emu/Oe)': [
                                "AC X'  (emu/Oe)"],
                             'Xpp (emu/Oe)': [
                                'AC X" (emu/Oe)']
                            }
        self.raw_df = None
        self.raw_df_header = None
        self.ppms_data_file = None
        self.meas_info = {}
        self.data_header_idx = 0
        self.num_meas_freqs = 0
        self.num_meas_temps = 0
        self.temp_subsets = []
        self.meas_temps = []
        self.Tmin, self.Tmax = 0,0
        
        self.temp_colormap = self.gui_default_colormap()
        
        self.Xd_capsule = -1.8*10**-8 # unit: emu/Oe
        self.Xd_film = 6.47*10**-10 # unit: emu/(Oe*mg)
        
        self.raw_data_fit = None
        
        """ Creating data analysis tab ----------------------------------------------------- """
        self.data_analysis_tab = QSplitter()
        self.all_the_tabs.addTab(self.data_analysis_tab, 'Analysis')
        
        ## Adding load controls
        self.load_wdgt = QWidget()
        self.load_layout = QVBoxLayout()
        self.load_layout.addStretch()
        
        self.see_fit_btn = QPushButton('Fitted params')
        self.see_fit_btn.clicked.connect(self.print_fitted_params)
        self.load_layout.addWidget(self.see_fit_btn)
        
        self.load_btn = QPushButton('Load')
        self.load_btn.clicked.connect(self.load_t_tau_data)
        self.load_layout.addWidget(self.load_btn)
        
        self.load_wdgt.setLayout(self.load_layout)
        
        """ Adding plotting controls """
        self.ana_plot = PlottingWindow()
        self.ana_plot.ax.set_xlabel('Temperature [$K^{-1}$]')
        self.ana_plot.ax.set_ylabel(r'$\ln{\tau}$ [$\ln{s}$]')
        
        # Adding fit controls
        self.fit_wdgt = QWidget()
        self.fit_layout = QVBoxLayout()
        
        self.cb_headline = QLabel('Fit types to consider')
        self.cb_headline.setFont(self.headline_font)
        self.fit_layout.addWidget(self.cb_headline)
        
        self.orbach_cb = QCheckBox('Orbach')
        self.orbach_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_layout.addWidget(self.orbach_cb)
        
        self.raman_cb = QCheckBox('Raman')
        self.raman_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_layout.addWidget(self.raman_cb)
        
        self.qt_cb = QCheckBox('QT')
        self.qt_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_layout.addWidget(self.qt_cb)
        
        # Adding temperature controls
        self.temp_headline = QLabel('Temperature interval')
        self.temp_headline.setFont(self.headline_font)
        self.fit_layout.addWidget(self.temp_headline)
        
        self.temp_horizontal_layout = QHBoxLayout()
        self.temp_line = [QLabel('('), QDoubleSpinBox(), QLabel(','), QDoubleSpinBox(), QLabel(')')]
        
        self.temp_line[1].setRange(0,self.temp_line[3].value())
        self.temp_line[1].setSingleStep(0.1)
        self.temp_line[3].setRange(self.temp_line[1].value(),1000)
        self.temp_line[3].setSingleStep(0.1)
        
        self.temp_line[1].editingFinished.connect(self.set_new_temp_ranges)
        self.temp_line[3].editingFinished.connect(self.set_new_temp_ranges)
        for w in self.temp_line:
            self.temp_horizontal_layout.addWidget(w)
        
        self.fit_layout.addLayout(self.temp_horizontal_layout)
        
        # Adding a button to run a fit
        self.run_fit_btn = QPushButton('Run fit!')
        self.run_fit_btn.clicked.connect(self.make_the_fit)
        self.fit_layout.addWidget(self.run_fit_btn)
        
        # Adding a list to hold information about simulations
        self.simulations_headline = QLabel('Simulations')
        self.simulations_headline.setFont(self.headline_font)
        self.fit_layout.addWidget(self.simulations_headline)
        
        self.list_of_simulations = QListWidget()
        self.list_of_simulations.doubleClicked.connect(self.edit_a_simulation)
        self.fit_layout.addWidget(self.list_of_simulations)
        
        # Adding buttons to control simulation list
        self.sim_btn_layout = QHBoxLayout()
        
        self.delete_sim_btn = QPushButton('Delete')
        self.delete_sim_btn.clicked.connect(self.delete_sim)
        self.sim_btn_layout.addWidget(self.delete_sim_btn)
        
        self.edit_sim_btn = QPushButton('Edit')
        self.edit_sim_btn.clicked.connect(self.edit_a_simulation)
        self.sim_btn_layout.addWidget(self.edit_sim_btn)
        
        self.new_sim_btn = QPushButton('New')
        self.new_sim_btn.clicked.connect(self.add_new_simulation)
        self.sim_btn_layout.addWidget(self.new_sim_btn)
        
        self.fit_layout.addLayout(self.sim_btn_layout)
        
        self.fit_wdgt.setLayout(self.fit_layout)
        
        # Finalizing layout of the data analysis tab
        self.data_analysis_tab.addWidget(self.load_wdgt)
        self.data_analysis_tab.addWidget(self.ana_plot)
        self.data_analysis_tab.addWidget(self.fit_wdgt)
        
        """ Creating the data treatment tab  ------------------------------------------------------------------------------- """
        self.data_treatment_tab = QSplitter()
        self.all_the_tabs.addTab(self.data_treatment_tab, 'Data treatment')
        
        ### Making the left column (data loading and visualization controls)
        self.data_loading_wdgt = QWidget()
        self.data_layout = QVBoxLayout()
        
        self.raw_data_load_btn = QPushButton('(1) Load')
        self.raw_data_load_btn.clicked.connect(self.load_ppms_data)
        self.data_layout.addWidget(self.raw_data_load_btn)
        
        self.diamag_correction_btn = QPushButton("(2) Diamagnetic correction")
        self.diamag_correction_btn.clicked.connect(self.calculate_diamag_correction)
        self.data_layout.addWidget(self.diamag_correction_btn)
        
        self.fit_Xp_Xpp_btn = QPushButton("(3) Fit X', X''")
        self.fit_Xp_Xpp_btn.clicked.connect(self.fit_Xp_Xpp_standalone)
        self.data_layout.addWidget(self.fit_Xp_Xpp_btn)
        
        self.copy_fit_to_ana_btn = QPushButton('(4) Copy for analysis')
        self.copy_fit_to_ana_btn.clicked.connect(self.copy_fit_to_analysis)
        self.data_layout.addWidget(self.copy_fit_to_ana_btn)
        
        self.save_fit_to_file_btn = QPushButton('(5) Save fit to file')
        self.save_fit_to_file_btn.clicked.connect(self.save_fit_to_file)
        self.data_layout.addWidget(self.save_fit_to_file_btn)
        
        ## Constructing data plotting layout
        self.raw_data_plot_lo = QVBoxLayout()
        
        self.analysis_plot_type_header = QLabel('Axis content')
        self.analysis_plot_type_header.setFont(self.headline_font)
        self.raw_data_plot_lo.addWidget(self.analysis_plot_type_header)
        
        self.analysis_plot_type_combo = QComboBox()
        self.analysis_plot_type_combo.addItems(['Raw data', 'Fitted'])
        self.analysis_plot_type_combo.currentIndexChanged.connect(self.switch_analysis_view)
        self.raw_data_plot_lo.addWidget(self.analysis_plot_type_combo)
        
        self.raw_data_plot_header = QLabel('Raw data plotting')
        self.raw_data_plot_header.setFont(self.headline_font)
        self.raw_data_plot_lo.addWidget(self.raw_data_plot_header)
        
        # Constructing the x combobox
        self.data_ana_x_lo = QHBoxLayout()
        self.analysis_x_combo_lbl = QLabel('x')
        self.data_ana_x_lo.addWidget(self.analysis_x_combo_lbl)
        
        self.analysis_x_combo = QComboBox()
        self.analysis_x_combo.currentIndexChanged.connect(self.plot_from_combo)
        self.data_ana_x_lo.addWidget(self.analysis_x_combo)
        self.raw_data_plot_lo.addLayout(self.data_ana_x_lo)
        
        # Constructing the y combobox
        self.data_ana_y_lo = QHBoxLayout()
        self.analysis_y_combo_lbl = QLabel('y')
        self.data_ana_y_lo.addWidget(self.analysis_y_combo_lbl)
        
        self.analysis_y_combo = QComboBox()
        self.analysis_y_combo.currentIndexChanged.connect(self.plot_from_combo)
        self.data_ana_y_lo.addWidget(self.analysis_y_combo)
        self.raw_data_plot_lo.addLayout(self.data_ana_y_lo)
        
        # Constructing a combobox for plotting fitted data
        self.fitted_data_plot_header = QLabel('Fit data plotting')
        self.fitted_data_plot_header.setFont(self.headline_font)
        self.raw_data_plot_lo.addWidget(self.fitted_data_plot_header)
        
        self.fit_data_plot_combo = QComboBox()
        self.fit_data_plot_combo.addItems(['ColeCole', 'FreqVSXpp'])
        self.fit_data_plot_combo.currentIndexChanged.connect(self.plot_from_itemlist)
        self.raw_data_plot_lo.addWidget(self.fit_data_plot_combo)
        
        self.data_layout.addLayout(self.raw_data_plot_lo)
        
        ## Finalizing the data loading widget
        self.data_layout.addStretch()
        self.data_loading_wdgt.setLayout(self.data_layout)
        
        # Making the middle of the tab (data visualization)
        self.treat_sw = QStackedWidget()
        
        self.treat_raw_plot = PlottingWindow()
        self.treat_fit_plot = PlottingWindow()
        
        self.treat_sw.addWidget(self.treat_raw_plot)
        self.treat_sw.addWidget(self.treat_fit_plot)
        
        ## Making the right column (parameter controls)
        self.param_wdgt = QWidget()
        self.param_layout = QVBoxLayout()
        
        # Sample mass
        self.sample_info_layout = QVBoxLayout()
        
        self.sample_info_header = QLabel('Sample information')
        self.sample_info_header.setFont(self.headline_font)
        self.sample_info_layout.addWidget(self.sample_info_header)
        
        # Sample mass edit
        self.sample_mass_layout = QHBoxLayout()
        self.sample_info_layout.addLayout(self.sample_mass_layout)
        
        self.sample_mass_lbl = QLabel('m (sample) [mg]')
        self.sample_mass_layout.addWidget(self.sample_mass_lbl)
        
        self.sample_mass_inp = QLineEdit()
        self.sample_mass_inp.setValidator(QDoubleValidator())
        self.sample_mass_layout.addWidget(self.sample_mass_inp)
        
        # Film mass edit
        self.film_mass_layout = QHBoxLayout()
        self.sample_info_layout.addLayout(self.film_mass_layout)
        
        self.film_mass_lbl = QLabel('m (film) [mg]')
        self.film_mass_layout.addWidget(self.film_mass_lbl)
             
        self.film_mass_inp = QLineEdit()
        self.film_mass_inp.setValidator(QDoubleValidator())
        self.film_mass_layout.addWidget(self.film_mass_inp)
        
        # Molar mass edit
        self.molar_mass_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.molar_mass_lo)
        
        self.molar_mass_lbl = QLabel('M [g/mol]')
        self.molar_mass_lo.addWidget(self.molar_mass_lbl)
        
        self.molar_mass_inp = QLineEdit()
        self.molar_mass_inp.setValidator(QDoubleValidator())
        self.molar_mass_lo.addWidget(self.molar_mass_inp)
        
        # Sample X' edit
        self.sample_xd_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.sample_xd_lo)
        
        self.sample_xd_lbl = QLabel(u"X\u1D05"+' (sample) [emu/(Oe*mol)]')
        self.sample_xd_lo.addWidget(self.sample_xd_lbl)
        
        self.sample_xd_inp = QLineEdit()
        self.sample_xd_inp.setText('-6e-7')
        self.sample_xd_inp.setValidator(QDoubleValidator())
        self.sample_xd_lo.addWidget(self.sample_xd_inp)
        
        # Mass load button
        self.load_mass_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.load_mass_lo)
        
        self.load_mass_btn = QPushButton('Load mass from file')
        self.load_mass_btn.clicked.connect(self.load_sample_film_mass)
        self.load_mass_lo.addWidget(self.load_mass_btn)
        self.load_mass_lo.addStretch()
        
        self.param_layout.addLayout(self.sample_info_layout)
        
        # List of fitted raw data
        self.treat_raw_fit_headline = QLabel('Fitted parameters')
        self.treat_raw_fit_headline.setFont(self.headline_font)
        self.param_layout.addWidget(self.treat_raw_fit_headline)
        
        self.treat_raw_fit_list = QListWidget()
        self.param_layout.addWidget(self.treat_raw_fit_list)
        self.treat_raw_fit_list.doubleClicked.connect(self.update_raw_plot_status)
        
        ## Finalizing layout
        
        self.param_layout.addStretch()
        self.param_wdgt.setLayout(self.param_layout)
        
        self.data_treatment_tab.addWidget(self.data_loading_wdgt)
        self.data_treatment_tab.addWidget(self.treat_sw)
        self.data_treatment_tab.addWidget(self.param_wdgt)
        
        """ Making a menubar -------------------------------------------------------------- """
        self.menu_bar = self.menuBar()
        
        # File menu
        self.file_menu = self.menu_bar.addMenu('File')
        
        self.quit_action = QAction('&Quit', self)
        self.quit_action.setShortcut("Ctrl+Q")
        self.quit_action.triggered.connect(sys.exit)
        self.file_menu.addAction(self.quit_action)
        
        # Simulation menu
        self.sim_menu = self.menu_bar.addMenu('Simulation')
        
        self.add_sim_w_menu = QAction('&New', self)
        self.add_sim_w_menu.setShortcut("Ctrl+Shift+N")
        self.add_sim_w_menu.triggered.connect(self.add_new_simulation)
        self.sim_menu.addAction(self.add_sim_w_menu)
        
        # About menu
        self.help_menu = self.menu_bar.addMenu('Help')
        
        self.help_about_menu = QAction('About', self)
        self.help_about_menu.triggered.connect(self.show_about_dialog)
        self.help_about_menu.setShortcut("F10")
        self.help_menu.addAction(self.help_about_menu)
        
        # Showing the GUI
        self.load_t_tau_data()
        
        self.showMaximized()
        self.show()
    
    def show_about_dialog(self):
    
        w = AboutDialog(info=self.about_information)
        w.exec_()
    
    def copy_fit_to_analysis(self):
    
        try:
            D = np.array(list(zip(self.meas_temps,
                                  self.raw_data_fit['Tau'],
                                  self.raw_data_fit['dTau'])))
            self.set_new_t_tau(D)
            self.read_indices_for_used_temps()
            self.plot_t_tau_on_axes()
            self.ana_plot.reset_axes()
        except TypeError:
            print('When the fitted data does not yet exist')
        
    
    def switch_analysis_view(self):
        
        idx = self.analysis_plot_type_combo.currentIndex()
        self.treat_sw.setCurrentIndex(idx)
        
    def update_treat_raw_fit_list(self):
        
        self.treat_raw_fit_list.clear()
        
        for i in range(self.num_meas_temps):
            newitem = QListWidgetItem()
            newitem.setText('{}, {}, {}'.format(self.meas_temps[i],True, True))
            plotting_dict = {'temp': self.meas_temps[i],
                             'raw': True,
                             'fit': True}
            newitem.setData(32, plotting_dict)
            self.treat_raw_fit_list.addItem(newitem)
    
    def fit_Xp_Xpp_standalone(self):
        
        if self.raw_df is None:
            pass
        
        else:
            self.statusBar.showMessage('Running fit...')
            
            #### Defining the objective function for lmfit ####
            def objective(params, v, data):
        
                Xp_data = data[0]
                Xpp_data = data[1]
                
                Xp_model = Xp_dataset(params, v)
                Xpp_model = Xpp_dataset(params, v)
                
                Xp_resid = Xp_data - Xp_model
                Xpp_resid = Xpp_data - Xpp_model
        
                return np.concatenate([Xp_resid, Xpp_resid])
            ########
            
            T, Xs, Xt, tau, alpha, resid, tau_fit_err = [],[],[],[],[],[],[]
            if 'Xp (emu/Oe)' in self.raw_df.columns:
                for t_idx in range(self.num_meas_temps):
                    v = np.array(self.temp_subsets[t_idx]['AC Frequency (Hz)'])
                    Xp = np.array(self.temp_subsets[t_idx]['Xp (emu/Oe)'])
                    Xpp = np.array(self.temp_subsets[t_idx]['Xpp (emu/Oe)'])
                    
                    data = [Xp, Xpp]
                    
                    tau_init = (v[np.argmax(Xpp)]*2*np.pi)**-1
                    fit_params = Parameters()
                    fit_params.add('Xs', value=Xp[-1], min=0, max=np.inf)
                    fit_params.add('Xt', value=Xp[0], min=0, max=np.inf)
                    fit_params.add('tau', value=tau_init, min=0, max=np.inf)
                    fit_params.add('alpha', value=0.1, min=0, max=np.inf)
                    
                    out = minimize(objective, fit_params, args=(v, data))
                    if 'Could not estimate error-bars.' in out.message:
                        tau_err = 0
                    else:
                        tau_idx = out.var_names.index('tau')
                        tau_err = np.sqrt(out.covar[tau_idx, tau_idx])
                    T.append(self.meas_temps[t_idx])
                    Xs.append(out.params['Xs'].value)
                    Xt.append(out.params['Xt'].value)
                    tau.append(out.params['tau'].value)
                    alpha.append(out.params['alpha'].value)
                    resid.append(out.residual.sum())
                    tau_fit_err.append(tau_err)
            fit_result = pd.DataFrame(data={'Temp': T,
                                            'ChiS': Xs,
                                            'ChiT': Xt,
                                            'Tau': tau,
                                            'Alpha': alpha,
                                            'Residual': resid,
                                            'Tau_Err': tau_fit_err,
                                            'dTau': tau_err_RC(tau, tau_fit_err, alpha)})
            
            self.raw_data_fit = fit_result
            self.update_treat_raw_fit_list()
            self.plot_from_itemlist()
            set_idx = self.analysis_plot_type_combo.findText('Fitted')
            self.analysis_plot_type_combo.setCurrentIndex(set_idx)
            
            self.statusBar.showMessage("Fit of X' and X'' complete")
        
    def get_single_line(self, filename, line_number):

        f = open(filename, 'r')
        for n in range(line_number):
            line = f.readline()
        f.close()
        
        return line
    
    def save_fit_to_file(self):
        
        name = QFileDialog.getSaveFileName(self, 'Save File')
        filename = name[0]
        if self.raw_data_fit is None:
            msg = QMessageBox()
            msg.setText('There is nothing to save')
            msg.setDetailedText("""There is probably no fit yet...""")
            msg.exec_()
        elif name=='':
            pass
            print('No file selected')
        else:
            df_to_save = self.raw_data_fit.copy()
            df_to_save = df_to_save.reindex(columns=['Temp', 'Tau', 'dTau', 'Alpha','ChiS', 'ChiT', 'Residual', 'Tau_Err'])
            df_to_save.sort_values('Temp', inplace=True)
            
            name, ext = os.path.splitext(filename)
            if ext == '':
                ext = '.dat'
            df_to_save.to_csv(name+'{}'.format(ext),
                              sep=';',
                              index=False,
                              float_format='%20.10e')
    
    def load_sample_film_mass(self):
    
        filename_info = QFileDialog().getOpenFileName(self, 'Open file', self.last_loaded_file)
        filename = filename_info[0]
        
        self.last_loaded_file = os.path.split(filename)[0]
        
        if filename == '':
            return 0
        else:
            sample = self.get_single_line(filename, 9).strip()
            film = self.get_single_line(filename, 10).strip()
        
        try:
            sample = sample.split(',')
            film = film.split(',')
            
            assert sample[0] == 'INFO' and (sample[1] == 's' or sample[1] == 'sample')
            assert film[0] == 'INFO' and (film[1] == 'f' or film[1] == 'film')
        except AssertionError:
            msg = QMessageBox()
            msg.setText('File not as expected')
            msg.setDetailedText("""Expected:
line 9: INFO,s,<mass>mg
line 10: INFO,f,<mass>mg""")
            msg.exec_()
        else:
            sample = float(sample[2][:-2])
            film = float(film[2][:-2])
            
            self.sample_mass_inp.setText(str(sample))
            self.film_mass_inp.setText(str(film))
    
    def calculate_diamag_correction(self):
        
        if self.raw_df is None:
            # Don't do the calculation, if there is nothing to calculate on
            pass
        
        else:
            try:
                sample_mass = float(self.sample_mass_inp.text())
                film_mass = float(self.film_mass_inp.text())
                molar_mass = float(self.molar_mass_inp.text())
                Xd_sample = float(self.sample_xd_inp.text())
            except ValueError as error:
                msg = QMessageBox()
                msg.setWindowTitle('Error')
                msg.setText('Could not read sample information')
                msg.exec_()
            else:
                H = self.raw_df['AC Amplitude (Oe)']
                H0 = self.raw_df['Magnetic Field (Oe)']
                Mp = self.raw_df["Mp (emu)"]
                Mpp = self.raw_df["Mpp (emu)"]
                
                Xp = (Mp - self.Xd_capsule*H - self.Xd_film*film_mass*H)*molar_mass/(sample_mass*H) - Xd_sample*molar_mass
                Xpp = Mpp/(sample_mass*H)*molar_mass
                
                Xp_idx = self.raw_df.columns.get_loc("M' (emu)")+1
                self.raw_df.insert(Xp_idx, column="X' (emu/(Oe*mol))", value=Xp)
                
                Xpp_idx = self.raw_df.columns.get_loc("M'' (emu)")+1
                self.raw_df.insert(Xpp_idx, column="X'' (emu/(Oe*mol))", value=Xpp)
        
            self.update_temp_subsets()
            self.update_analysis_combos()
    
        #elif "X'' (emu/(Oe*mol))" in self.raw_df.columns:
        #    # Don't add an element that is already there
        #    pass
        #elif "AC X'  (emu/Oe)" in self.raw_df.columns:
        #    # Data was read from source where the column already exists.
        #    # Update data names
        #    self.data_names.update({'Xp': "AC X'  (emu/Oe)",
        #                            'Xpp': 'AC X" (emu/Oe)'})
        
    
    def update_itemdict(self, item, itemdict):
        
        item.setData(32, itemdict)
        item.setText('{}, {}, {}'.format(itemdict['temp'],
                                         itemdict['raw'],
                                         itemdict['fit']))
    
    def update_raw_plot_status(self):
        
        w = FitResultPlotStatus(list_input=self.treat_raw_fit_list)
        finished_value = w.exec_()
        if not finished_value:
            pass
        else:
            final_states = w.checked_items
            
            for i, boxes in enumerate(final_states):
                item = self.treat_raw_fit_list.item(i)
                item_data = item.data(32)
                item_data['raw'] = boxes[0].isChecked()
                item_data['fit'] = boxes[1].isChecked()
                self.update_itemdict(item, item_data)
            
        self.plot_from_itemlist()
        self.treat_fit_plot.canvas.draw()
    
    def plot_from_itemlist(self):
        
        if self.treat_raw_fit_list.count()==0:
            return
        else:
            self.treat_fit_plot.ax.clear()
            plot_type = self.fit_data_plot_combo.currentText()
            
            if plot_type == 'FreqVSXpp':
                x_name = 'AC Frequency (Hz)'
                y_name = 'Xpp (emu/Oe)'
                for row in range(self.num_meas_temps):
                
                    T = self.meas_temps[row]
                    rgb = self.temp_colormap((T-self.Tmin)/(self.Tmax-self.Tmin))
                    
                    item = self.treat_raw_fit_list.item(row)
                    itemdict = item.data(32)
                    
                    if itemdict['raw']:
                        self.treat_fit_plot.ax.plot(self.temp_subsets[row][x_name],
                                                    self.temp_subsets[row][y_name],
                                                    marker='o',
                                                    mec='k',
                                                    mfc='none',
                                                    linestyle='None')
                    if itemdict['fit']:
                        self.treat_fit_plot.ax.plot(self.temp_subsets[row][x_name],
                                                    Xpp_(self.temp_subsets[row]['AC Frequency (Hz)'],
                                                        self.raw_data_fit['ChiS'].iloc[row],
                                                        self.raw_data_fit['ChiT'].iloc[row],
                                                        self.raw_data_fit['Tau'].iloc[row],
                                                        self.raw_data_fit['Alpha'].iloc[row]),
                                                    c=rgb)
                self.treat_fit_plot.ax.set_xscale('log')
                self.treat_fit_plot.ax.set_xlabel(x_name)
                self.treat_fit_plot.ax.set_ylabel(y_name)
        
            elif plot_type == 'ColeCole':
                x_name = 'Xp (emu/Oe)'
                y_name = 'Xpp (emu/Oe)'
                for row in range(self.num_meas_temps):

                    T = self.meas_temps[row]
                    rgb = self.temp_colormap((T-self.Tmin)/(self.Tmax-self.Tmin))
                    
                    item = self.treat_raw_fit_list.item(row)
                    itemdict = item.data(32)
                    
                    if itemdict['raw']:
                        self.treat_fit_plot.ax.plot(self.temp_subsets[row][x_name],
                                                self.temp_subsets[row][y_name],
                                                marker='o',
                                                mec='k',
                                                mfc='none',
                                                linestyle='None')
                    if itemdict['fit']:
                        self.treat_fit_plot.ax.plot(Xp_(self.temp_subsets[row]['AC Frequency (Hz)'],
                                                        self.raw_data_fit['ChiS'].iloc[row],
                                                        self.raw_data_fit['ChiT'].iloc[row],
                                                        self.raw_data_fit['Tau'].iloc[row],
                                                        self.raw_data_fit['Alpha'].iloc[row]),
                                                    Xpp_(self.temp_subsets[row]['AC Frequency (Hz)'],
                                                        self.raw_data_fit['ChiS'].iloc[row],
                                                        self.raw_data_fit['ChiT'].iloc[row],
                                                        self.raw_data_fit['Tau'].iloc[row],
                                                        self.raw_data_fit['Alpha'].iloc[row]),
                                                    c=rgb)
                                                    
                self.treat_fit_plot.ax.set_xlabel(x_name)
                self.treat_fit_plot.ax.set_ylabel(y_name)
            
            self.treat_fit_plot.canvas.draw()
        
    def plot_from_combo(self):
        
        self.treat_raw_plot.ax.clear()
        
        idx_x = self.analysis_x_combo.currentIndex()
        idx_y = self.analysis_y_combo.currentIndex()
        
        x_label = self.raw_df.columns[idx_x]
        y_label = self.raw_df.columns[idx_y]
        
        self.treat_raw_plot.ax.plot(self.raw_df[x_label],
                                    self.raw_df[y_label],
                                    marker='o',
                                    mfc='none',
                                    mec='k',
                                    linestyle='-',
                                    c='k',
                                    linewidth=1,
                                    )
        self.treat_raw_plot.ax.set_xlabel(x_label)
        self.treat_raw_plot.ax.set_ylabel(y_label)
        
        self.treat_raw_plot.canvas.draw()
    
    def fill_df_data_values(self):
    
        if ('Xp (emu/Oe)' in self.raw_df.columns and not ('Mp (emu)' in self.raw_df.columns)):
            # Susceptibility exists in the data frame, but magnetisation does not
            Mp = self.raw_df['Xp (emu/Oe)']*self.raw_df['Magnetic Field (Oe)']
            Mpp = self.raw_df['Xpp (emu/Oe)']*self.raw_df['Magnetic Field (Oe)']
            Xp_idx = self.raw_df.columns.get_loc('Xp (emu/Oe)')
            self.raw_df.insert(Xp_idx, column='Mp (emu)', value=Mp)
            self.raw_df.insert(Xp_idx, column='Mpp (emu)', value=Mpp)
            
        elif (not 'Xp (emu/Oe)' in self.raw_df.columns and ('Mp (emu)' in self.raw_df.columns)):
            # Magnetisation exists in the data frame, but susceptibility does not
            Xp = self.raw_df['Mp (emu)']/self.raw_df['Magnetic Field (Oe)']
            Xpp = self.raw_df['Mpp (emu)']/self.raw_df['Magnetic Field (Oe)']
            Mp_idx = self.raw_df.columns.get_loc('Mp (emu)')
            self.raw_df.insert(Mp_idx+2, column='Xp (emu/Oe)', value=Xp)
            self.raw_df.insert(Mp_idx+3, column='Xpp (emu/Oe)', value=Xpp)
    
    def load_ppms_data(self):
        
        open_file_dialog = QFileDialog()
        filename_info = open_file_dialog.getOpenFileName(self, 'Open file', self.last_loaded_file)
        filename = filename_info[0]
        
        try:
            # FileNotFoundError and UnicodeDecodeError will be raised here
            potential_header, potential_df = read_ppms_file(filename)
            if potential_header is None:
                raise FileFormatError(filename)
            summary = update_data_names(potential_df, self.read_options)
            counts = [val>1 for key, val in summary.items()]
            # To make sure that none of the names in read_options were matched more than once.
            assert not any(counts)
            # To make sure that only Mp (and therefore Mpp) OR Xp (and therefore Xpp) can appear at once.
            # In the case that this is ever an error, self.fill_df_data_values will have to be changed.
            assert (summary['Mp (emu)']>0) != (summary['Xp (emu/Oe)']>0)
        
        except FileNotFoundError:
            # Did not read any file
            pass
        except UnicodeDecodeError:
            # File being read is binary, not a text file
            print('The file is not a text file')
        except FileFormatError as e:
            # File does not have correct header and data blocks
            print('Trying to read a file that does not look correct')
        except AssertionError:
            # The data names could not be mapped correctly
            print('A data name from self.read_options is showing up more than once in the columns')
            print('OR')
            print('both Mp and Xp are unexpectedly both showing up in the data names')
        
        else:
            # Now that everything has been seen to work,
            # save potential header and potential df as actual header and df
            self.last_loaded_file = os.path.split(filename)[0]
            self.ppms_data_file = filename
            self.raw_df = potential_df
            self.raw_df_header = potential_header
            
            # Clear old data and set new names
            self.treat_raw_fit_list.clear()
            self.fill_df_data_values()
            
            self.cleanup_loaded_ppms()
            self.num_meas_freqs = len(set(self.raw_df['AC Frequency (Hz)']))
            self.update_temp_subsets()
            self.update_meas_temps()
            
            # Clearing axes of "old" drawings and setting front widget to the raw data
            self.treat_raw_plot.clear_canvas()
            self.treat_fit_plot.clear_canvas()
            
            combo_idx = self.analysis_plot_type_combo.findText('Raw data')
            self.analysis_plot_type_combo.setCurrentIndex(combo_idx)
            
            # Updating analysis combos, which will automatically draw the new data
            self.update_analysis_combos()
 
    def update_meas_temps(self):
        
        meas_temps = []
        for sub in self.temp_subsets:
            meas_temps.append(sub['Temperature (K)'].mean())
        
        self.meas_temps = np.array(meas_temps)
        self.num_meas_temps = len(self.meas_temps)
        self.Tmin = self.meas_temps.min()
        self.Tmax = self.meas_temps.max()
        
    def update_temp_subsets(self):
        
        self.temp_subsets = []
        idx_list = [-1]
        # Splits based on where the frequency "restarts"
        # Assumes that the frequency is always increasing within a measurement
        i=0
        old_val = 0
        while i<self.raw_df.shape[0]:
            new_val = self.raw_df['AC Frequency (Hz)'].iloc[i]
            if new_val<old_val:
                idx_list.append(i)
            else:
                pass
            old_val = new_val
            i+=1
        idx_list.append(self.raw_df.shape[0])
        for n in range(len(idx_list)-1):
            self.temp_subsets.append(self.raw_df.iloc[idx_list[n]+1:idx_list[n+1]])
    
    def cleanup_loaded_ppms(self):
        
        # Drop columns where all values are NaN
        self.raw_df.dropna(axis=1, how='all', inplace=True)
        # Removing instrument comment lines
        # Drop "Comment" column
        if 'Comment' in self.raw_df.columns:
            self.raw_df.drop(['Comment'], axis='columns', inplace=True)
        # Drop all rows where there is still a NaN value
        self.raw_df.dropna(axis=0, inplace=True)
        
        # Make sure that the rows are named continuously
        old_indices = self.raw_df.index.values
        new_indices = list(range(len(old_indices)))
        self.raw_df.rename(index=dict(zip(old_indices, new_indices)),
                           inplace=True)
       
    def update_analysis_combos(self):
    
        self.analysis_x_combo.clear()
        self.analysis_x_combo.addItems(self.raw_df.columns)
        
        self.analysis_y_combo.clear()
        self.analysis_y_combo.addItems(self.raw_df.columns)
    
    def print_fitted_params(self):
        
        if self.fitted_parameters == None:
            pass
        else:
            dialog = ParamDialog(param_dict=self.fitted_parameters)
            finished = dialog.exec_()
    
    def add_new_simulation(self):
    
        sim_dialog = SimulationDialog(fitted_parameters=self.fitted_parameters,
                                      plot_type_list=[],
                                      plot_parameters={'tQT': 0.01, 'Cr': 0.01, 'n': 0.01, 't0': 0.01, 'Ueff': 0.01},
                                      min_and_max_temps=[0,0])
        
        finished_value = sim_dialog.exec_()
        
        if finished_value:
            
            plot_type = sim_dialog.plot_type_list

            if len(plot_type)<1:
                pass
            else:
                p_fit = sim_dialog.plot_parameters
                T_vals = sim_dialog.min_and_max_temps
                
                plot_to_make = ''.join(plot_type)
                new_item_text = '{}, ({},{}), tQT: {}, Cr: {}, n: {}, t0: {}, Ueff: {}'.format(
                                plot_type, T_vals[0], T_vals[1], p_fit['tQT'], p_fit['Cr'],
                                p_fit['n'], p_fit['t0'], p_fit['Ueff'])
                
                new_list_item = QListWidgetItem()
                
                line = addPartialModel(self.ana_plot.fig,
                                        T_vals[0],
                                        T_vals[1],
                                        self.prepare_sim_dict_for_plotting(p_fit),
                                        plotType=plot_to_make)
                                        
                list_item_data = {'plot_type': plot_type,
                                  'p_fit': p_fit,
                                  'T_vals': T_vals,
                                  'line': line}
                
                self.list_of_simulations.addItem(new_list_item)
                new_list_item.setText(new_item_text)
                new_list_item.setData(32, list_item_data)
                
                self.ana_plot.canvas.draw()
                
        else:
            pass
    
    def edit_a_simulation(self):
        
        try:
            sim_item = self.list_of_simulations.selectedItems()[0]
        except IndexError:
            pass
        else:
            # Reading off information from the selected item
            old_data = sim_item.data(32)
            old_plot_type_input = old_data['plot_type']
            old_p_fit = old_data['p_fit']
            old_T_vals = old_data['T_vals']
            old_line = old_data['line']
            
            
            
            # Opening simulation dialog with old parameters
            sim_dialog = SimulationDialog(fitted_parameters=self.fitted_parameters,
                                          plot_type_list=old_plot_type_input,
                                          plot_parameters=old_p_fit,
                                          min_and_max_temps=old_T_vals)
            
            finished_value = sim_dialog.exec_()
            
            if finished_value:
                # Reading new parameters of simulation
                new_plot_type = sim_dialog.plot_type_list
                
                
                if len(new_plot_type)<1:
                    pass
                else:
                    new_p_fit = sim_dialog.plot_parameters
                    new_T_vals = sim_dialog.min_and_max_temps
                    
                    plot_to_make = ''.join(new_plot_type)
                    new_item_text = '{}, ({},{}), tQT: {}, Cr: {}, n: {}, t0: {}, Ueff: {}'.format(
                                    new_plot_type, new_T_vals[0], new_T_vals[1], new_p_fit['tQT'], new_p_fit['Cr'],
                                    new_p_fit['n'], new_p_fit['t0'], new_p_fit['Ueff'])
                    
                    self.ana_plot.ax.lines.remove(old_line)
                    
                    new_line = addPartialModel(self.ana_plot.fig,
                                               new_T_vals[0],
                                               new_T_vals[1],
                                               self.prepare_sim_dict_for_plotting(new_p_fit),
                                               plotType=plot_to_make)
                    
                    list_item_data = {'plot_type': new_plot_type,
                                      'p_fit': new_p_fit,
                                      'T_vals': new_T_vals,
                                      'line': new_line}
                    
                    self.ana_plot.canvas.draw()
                    
                    sim_item.setData(32, list_item_data)
                    sim_item.setText(new_item_text)
            else:
                pass
        
    def delete_sim(self):
        
        try:
            sim_item = self.list_of_simulations.selectedItems()[0]
        except IndexError:
            pass
        else:
            line_pointer = sim_item.data(32)['line']
            self.ana_plot.ax.lines.remove(line_pointer)
            self.ana_plot.canvas.draw()
            
            item_row = self.list_of_simulations.row(sim_item)
            sim_item = self.list_of_simulations.takeItem(item_row)
            
            del sim_item
        
    def plot_t_tau_on_axes(self):
        
        if self.plotted_data_pointers is not None:
            for line in self.plotted_data_pointers:
                line.remove()
        self.plotted_data_pointers = []
        
        if self.data_dtau is None:
            used, = self.ana_plot.ax.plot(1/self.used_T, np.log(self.used_tau), 'bo')
            not_used, = self.ana_plot.ax.plot(1/self.not_used_T, np.log(self.not_used_tau), 'ro')
            self.plotted_data_pointers.append(used)
            self.plotted_data_pointers.append(not_used)
        else:
            err_used_point, caplines1, barlinecols1 = self.ana_plot.ax.errorbar(1/self.used_T, np.log(self.used_tau), yerr=self.used_dtau, fmt='bo', ecolor='b')
            err_not_used_point, caplines2, barlinecols2 = self.ana_plot.ax.errorbar(1/self.not_used_T, np.log(self.not_used_tau), yerr=self.not_used_dtau, fmt='ro', ecolor='r')

            self.plotted_data_pointers.append(err_used_point)
            self.plotted_data_pointers.append(err_not_used_point)
            for e in [caplines1, caplines2, barlinecols1, barlinecols2]:
                for line in e:
                    self.plotted_data_pointers.append(line)

        self.ana_plot.canvas.draw()
    
    def reset_analysis_containers(self):
        self.data_file_name = None
        self.data_T = None
        self.data_tau = None
        self.data_dtau = None
        
        self.used_T = None
        self.not_used_T = None
        self.used_tau = None
        self.not_used_tau = None
        self.used_dtau = None
        self.not_used_dtau = None
        
        self.fitted_parameters = None
        self.used_indices = None
        
    def gui_default_colormap(self):
        
        colormap = LinearSegmentedColormap.from_list('temp_colormap', [(0,   (0, 0, 0.5)),
                                                                       (1/8, (0, 0.0, 1)),
                                                                       (2/8, (0, 0.5, 1)),
                                                                       (3/8, (0.0, 1, 1.0)),
                                                                       (4/8, (0.5, 1, 0.5)),
                                                                       (5/8, (1, 1.0, 0)),
                                                                       (6/8, (1, 0.5, 0)),
                                                                       (7/8, (1.0, 0, 0)),
                                                                       (8/8, (0.5, 0, 0))])
        
        return colormap
        
    def load_t_tau_data(self):
        
        if self.startUp:
            try:
                filename = sys.argv[1]
            except IndexError:
                pass
            finally:
                self.startUp = False
                return 0
        else:
            filename_info = QFileDialog().getOpenFileName(self, 'Open file', self.last_loaded_file)
            filename = filename_info[0]
            
            self.last_loaded_file = os.path.split(filename)[0]
        
        if filename == '':
            pass
        else:
            self.reset_analysis_containers()
            self.data_file_name = filename
            self.data_file_dir = os.path.dirname(filename)
        
        try:
            D = np.loadtxt(self.data_file_name,
                           skiprows=1,
                           delimiter=';')
        except (ValueError, OSError) as error:
            sys.stdout.flush()
            if error_type == 'ValueError':
                msg = QMessageBox()
                msg.setWindowTitle('ValueError')
                msg.setIcon(QMessageBox.Warning)
                msg.setText('File format is not as expected')
                msg.exec_()
            elif error_type == 'OSError':
                pass
        else:
            self.set_new_t_tau(D)
            self.read_indices_for_used_temps()
            self.plot_t_tau_on_axes()
            self.ana_plot.reset_axes()
    
    def set_new_t_tau(self, D):
        """
        Uses the array D to set new values for T, tau, and alpha
        Assumes that the first column is temperatures, second column is tau-values
        and third column is error on tau.
        """
        
        T = D[:,0]
        tau = D[:,1]
        
        sort_indices = T.argsort()
        self.data_T = T[sort_indices]
        self.data_tau = tau[sort_indices]
        self.data_dtau = None
        
        if D.shape[1]>2:
            dtau = D[:,2]
            self.data_dtau = dtau[sort_indices]
            
    def prepare_sim_dict_for_plotting(self, p_fit_gui_struct):

        params = []
        quantities = []
        sigmas = [0]*5

        for key, val in p_fit_gui_struct.items():
            params.append(val)
            quantities.append(key)
        
        Ueff = params[quantities.index('Ueff')]
        params[quantities.index('Ueff')] = Ueff*sc.Boltzmann
        
        p_fit_script_type = {'params': params,
                             'quantities': quantities,
                             'sigmas': sigmas}
        
        return p_fit_script_type
    
    def read_fit_type_cbs(self):
    
        list_of_checked = []
        if self.qt_cb.isChecked(): list_of_checked.append('QT')
        if self.raman_cb.isChecked(): list_of_checked.append('R')
        if self.orbach_cb.isChecked(): list_of_checked.append('O')
        fitToMake = ''.join(list_of_checked)
        
        return fitToMake
    
    def read_indices_for_used_temps(self):
        
        min_t = self.temp_line[1].value()
        max_t = self.temp_line[3].value()
        
        try:
            self.used_indices = [list(self.data_T).index(t) for t in self.data_T if t>=min_t and t<=max_t]
            
            self.used_T = self.data_T[self.used_indices]
            self.used_tau = self.data_tau[self.used_indices]
            
            self.not_used_T = np.delete(self.data_T, self.used_indices)
            self.not_used_tau = np.delete(self.data_tau, self.used_indices)
            
            if self.data_dtau is not None:
                self.used_dtau = self.data_dtau[self.used_indices]
                self.not_used_dtau = np.delete(self.data_dtau, self.used_indices)
            
        except (AttributeError, TypeError):
            print('No data have been selected yet!')
    
    def fit_relaxation(self):
        
        guess = getParameterGuesses(self.used_T, self.used_tau)
        guess_dialog = GuessDialog(guess=guess)
        guess_dialog.exec_()
        try:
            guess = guess_dialog.return_guess
        except AttributeError:
            guess = guess
        perform_this_fit = self.read_fit_type_cbs()
        
        f = getFittingFunction(fitType=perform_this_fit)
        p0 = getStartParams(guess, fitType=perform_this_fit)
        
        if self.used_dtau is None:
            popt, pcov = curve_fit(f, self.used_T, np.log(self.used_tau), p0)
        else:
            popt, pcov = curve_fit(f, self.used_T, np.log(self.used_tau), p0, sigma=self.used_dtau)
        p_fit = readPopt(popt, pcov, fitType=perform_this_fit)
        
        return p_fit
        
    def set_new_temp_ranges(self):
    
        new_max_for_low = self.temp_line[3].value()
        new_min_for_high = self.temp_line[1].value()
        self.temp_line[1].setRange(0,new_max_for_low)
        self.temp_line[3].setRange(new_min_for_high,1000)
        
        self.read_indices_for_used_temps()
        if self.data_T is not None:
            self.plot_t_tau_on_axes()
        
    def make_the_fit(self):
        
        try:
            Tmin = self.temp_line[1].value()
            Tmax = self.temp_line[3].value()
            perform_this_fit = self.read_fit_type_cbs()
            assert Tmin != Tmax
            assert perform_this_fit != ''
            
            p_fit = self.fit_relaxation()
            self.fitted_parameters = p_fit
        
        except (AssertionError, RuntimeError, ValueError) as error:
            
            error_type = error.__class__.__name__
            if error_type == 'AssertionError':
                msg_text = 'Check your temperature and fit settings'
                msg_details = """Possible errors:
 - min and max temperatures are the same
 - no fit options have been selected"""
            elif error_type == 'RuntimeError':
                msg_text = 'This fit cannot be made within the set temperatures'
                msg_details = ''
            elif error_type == 'ValueError':
                msg_text = 'No file has been loaded'
                msg_details = ''
            
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle('Fit aborted')
            msg.setText(msg_text)
            msg.setDetailedText(msg_details)
            msg.exec_()
    
if __name__ == '__main__':
    
    myappid = 'AC Processing v1.0'
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    
    app = QApplication(sys.argv)
    w = ACGui()
    sys.exit(app.exec_())
        
        
        