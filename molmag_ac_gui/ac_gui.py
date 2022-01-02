#std packages
import sys
import os
import json
from importlib.resources import read_text
from multiprocessing import Pool
from collections import deque
import datetime

#third-party packages
import numpy as np
import pandas as pd
import names

import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, to_hex
from matplotlib._color_data import TABLEAU_COLORS

import scipy.constants as sc
from scipy.optimize import minimize, curve_fit
from lmfit import Parameter, Parameters, Minimizer, minimize

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon, QFont, QDoubleValidator, QColor
from PyQt5.QtWidgets import (QMainWindow, QWidget, QApplication, QPushButton,
                             QLabel, QAction, QComboBox, QStackedWidget,
                             QDoubleSpinBox, QFormLayout, QCheckBox,
                             QVBoxLayout, QMessageBox, QSplitter, QGridLayout,
                             QHBoxLayout, QFileDialog, QDialog, QLineEdit,
                             QListWidget, QListWidgetItem, QTabWidget,
                             QScrollArea, QStatusBar, QInputDialog,
                             QActionGroup)


#local imports
from .__init__ import __version__
from .process_ac import (Xp_, Xpp_, Xp_dataset, Xpp_dataset, add_partial_model,
                         getParameterGuesses, getStartParams,
                         getFittingFunction, readPopt, addPartialModel,
                         tau_err_RC, diamag_correction, fit_Xp_Xpp_genDebye,
                         _QT, _R, _O, general_relaxation, relaxation_dataset,
                         fit_relaxation, default_parameters)
from .dialogs import (GuessDialog, SimulationDialog, AboutDialog, ParamDialog,
                      FitResultPlotStatus, PlottingWindow, MagMessage)
from .utility import (read_ppms_file, get_ppms_column_name_matches,
                      update_data_names)
from .exceptions import FileFormatError, NoGuessExistsError
from . import data as pkg_static_data

#Datatabletab import:
from .DataTableTab import DataTableTab 


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
        self.about_information = {'author':
                                  'Emil A. Klahn (eklahn@chem.au.dk)',
                                  'webpage':
                                  'https://chem.au.dk/en/molmag',
                                  'personal':
                                  'https://eandklahn.github.io'
                                  }
        
        self.startUp = True #Some old stuff about reading a filee directly from the terminal. Does not work fully. 
        self.last_loaded_file = os.getcwd() #Remember the last used folder.
        self.current_file = ''

        
        """ Things to do with how the window is shown """
        self.setWindowTitle('Molmag AC GUI v{}'.format(__version__))
        self.setWindowIcon(QIcon('double_well_potential_R6p_icon.ico'))
        
        """ FONT STUFF """
        self.headline_font = QFont() #Defines headline font
        self.headline_font.setBold(True) #Makes it bold
        
        """ Setting up the main tab widget """
        self.all_the_tabs = QTabWidget() 
        self.setCentralWidget(self.all_the_tabs)
        self.statusBar = QStatusBar() #A statusbar in the bottom of the window
        self.setStatusBar(self.statusBar)
        
        """ Constructing the data analysis tab """
        # Data containers for analysis
        
        self.simulation_colors = [x for x in TABLEAU_COLORS]
        self.simulation_colors.remove('tab:gray') #Hvorfor fjerne grå og hvorfor tab: 
        self.simulation_colors = deque(self.simulation_colors) #Hvorfor deque i stedet for liste? 
        
        #Creating all the needed attributes or whatever it is called
        self.fit_history = list()

        self.data_T = None
        self.data_tau = None
        self.data_dtau = None
        
        self.plotted_data_pointers = None
        self.data_used_pointer = None
        self.data_not_used_pointer = None
        
        self.used_indices = None
        
        self.used_T = None
        self.not_used_T = None
        
        self.used_tau = None
        self.not_used_tau = None
        
        self.used_dtau = None
        self.not_used_dtau = None
        
        # Data containers for treatment
        self.read_options = json.loads(read_text(pkg_static_data,
                                                'read_options.json'))
        #JSON: JavaScript Object Notation. Vil du forklare lidt om det? Hvordan det bruges
        
        self.diamag_constants = json.loads(read_text(pkg_static_data,
                                                    'diamag_constants.json'))
        
        self.temperature_cmap = LinearSegmentedColormap.from_list(
            'temp_colormap',
            json.loads(read_text(pkg_static_data, 'default_colormap.json')))
        
        self.tooltips_dict = json.loads(read_text(pkg_static_data,
                                                  'tooltips.json'))

        #self.setStyleSheet(read_text(pkg_static_data, 'styles.qss'))

        #Creates dataframe with header and raw_df
        self.raw_df = None
        self.raw_df_header = None
        self.num_meas_freqs = 0
        self.num_meas_temps = 0
        self.temp_subsets = []
        self.meas_temps = []
        self.Tmin, self.Tmax = 0,0
        
        self.raw_data_fit = None
        
        """ Creating data analysis tab """
        self.data_analysis_tab = QSplitter()
        self.all_the_tabs.addTab(self.data_analysis_tab, 'Analysis')
        
        ## Adding load controls
        self.load_wdgt = QWidget()
        self.load_layout = QVBoxLayout()
        self.load_layout.addStretch()
        
        self.see_fit_btn = QPushButton('Fitted params')
        self.see_fit_btn.clicked.connect(self.show_fitted_params)
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
        self.temp_line = [QLabel('('), QDoubleSpinBox(), QLabel(','),
                          QDoubleSpinBox(), QLabel(')')]
        
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
        self.list_of_simulations.setDragDropMode(self.list_of_simulations.InternalMove)
        
        self.list_of_simulations.doubleClicked.connect(self.edit_simulation_from_list)
        """https://stackoverflow.com/questions/41353653/how-do-i-get-the-checked-items-in-a-qlistview"""
        self.list_of_simulations.itemChanged.connect(self.redraw_simulation_lines)
        
        self.fit_layout.addWidget(self.list_of_simulations)
        
        # Adding buttons to control simulation list
        self.sim_btn_layout = QHBoxLayout()
        
        self.delete_sim_btn = QPushButton('Delete')
        self.delete_sim_btn.clicked.connect(self.delete_sim)
        self.sim_btn_layout.addWidget(self.delete_sim_btn)
        
        self.edit_sim_btn = QPushButton('Edit')
        self.edit_sim_btn.clicked.connect(self.edit_simulation_from_list)
        self.sim_btn_layout.addWidget(self.edit_sim_btn)
        
        self.new_sim_btn = QPushButton('New')
        self.new_sim_btn.clicked.connect(self.edit_simulation_from_list)
        self.sim_btn_layout.addWidget(self.new_sim_btn)
        
        self.fit_layout.addLayout(self.sim_btn_layout)
        
        self.fit_wdgt.setLayout(self.fit_layout)
        
        # Finalizing layout of the data analysis tab
        self.data_analysis_tab.addWidget(self.load_wdgt)
        self.data_analysis_tab.addWidget(self.ana_plot)
        self.data_analysis_tab.addWidget(self.fit_wdgt)
        
        """ Creating the data treatment tab """
        self.data_treatment_tab = QSplitter()
        self.all_the_tabs.addTab(self.data_treatment_tab, 'Data treatment')
        
        ### Making the left column (data loading and visualization controls)
        self.data_loading_wdgt = QWidget()
        self.data_layout = QVBoxLayout()
        
        self.raw_data_load_btn = QPushButton('(1) Load')
        self.raw_data_load_btn.clicked.connect(self.load_ppms_data)
        self.data_layout.addWidget(self.raw_data_load_btn)
        
        self.diamag_correction_btn = QPushButton("(2) Diamagnetic correction")
        self.diamag_correction_btn.clicked.connect(self.make_diamag_correction_calculation)
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
        self.fit_data_plot_combo.addItems(['ColeCole', 'FreqVSXp', 'FreqVSXpp'])
        self.fit_data_plot_combo.currentIndexChanged.connect(self.plot_from_itemlist)
        self.raw_data_plot_lo.addWidget(self.fit_data_plot_combo)
        
        # Checkbox for coloring measured points
        self.fit_data_color_cb_lo = QHBoxLayout()
        self.fit_data_color_cb_header = QLabel('Use temp. colors')
        self.fit_data_color_cb_lo.addWidget(self.fit_data_color_cb_header)
        
        self.fit_data_color_cb = QCheckBox()
        self.fit_data_color_cb.stateChanged.connect(self.plot_from_itemlist)
        self.fit_data_color_cb_lo.addWidget(self.fit_data_color_cb)
        
        self.raw_data_plot_lo.addLayout(self.fit_data_color_cb_lo)
        
        # Finalizing the raw data layout
        self.data_layout.addLayout(self.raw_data_plot_lo)
        
        ## Finalizing the data loading widget
        self.data_layout.addStretch()
        self.data_loading_wdgt.setLayout(self.data_layout)
        
        # Making the middle of the tab (data visualization)
        self.treat_sw = QStackedWidget()
        
        self.treat_raw_plot = PlottingWindow()
        self.treat_fit_plot = PlottingWindow(make_cax=True)
        
        self.treat_sw.addWidget(self.treat_raw_plot)
        self.treat_sw.addWidget(self.treat_fit_plot)
        
        ## Making the right column (parameter controls)
        self.param_wdgt = QWidget()
        self.param_layout = QVBoxLayout()
        
        # SAMPLE INFO
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

        # Sample molar mass edit
        self.molar_mass_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.molar_mass_lo)
        
        self.molar_mass_lbl = QLabel('M (sample) [g/mol]')
        self.molar_mass_lo.addWidget(self.molar_mass_lbl)
        
        self.molar_mass_inp = QLineEdit()
        self.molar_mass_inp.setValidator(QDoubleValidator())
        self.molar_mass_lo.addWidget(self.molar_mass_inp)

        # Sample Xd edit
        self.sample_xd_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.sample_xd_lo)
        
        self.sample_xd_lbl = QLabel(u"<a href={}>X\u1D05</a>".format(
                                    self.tooltips_dict['English']['Xd_link'])
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
        
        
        # Mass load button
        self.sample_data_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.sample_data_lo)
        
        self.load_sample_data_btn = QPushButton('Load sample data')
        self.load_sample_data_btn.clicked.connect(self.load_sample_data)
        self.sample_data_lo.addWidget(self.load_sample_data_btn)
        
        self.save_sample_data_btn = QPushButton('Save sample data')
        self.save_sample_data_btn.clicked.connect(self.save_sample_data)
        self.sample_data_lo.addWidget(self.save_sample_data_btn)
        
        self.sample_data_lo.addStretch()
        
        self.param_layout.addLayout(self.sample_info_layout)
        
        # List of fitted raw data
        self.treat_raw_fit_headline = QLabel('Fitted parameters')
        self.treat_raw_fit_headline.setFont(self.headline_font)
        self.param_layout.addWidget(self.treat_raw_fit_headline)
        
        self.treat_raw_fit_list = QListWidget()
        self.param_layout.addWidget(self.treat_raw_fit_list)
        self.treat_raw_fit_list.doubleClicked.connect(self.update_raw_plot)
        
        ## Finalizing layout
        
        self.param_layout.addStretch()
        self.param_wdgt.setLayout(self.param_layout)
        
        self.data_treatment_tab.addWidget(self.data_loading_wdgt)
        self.data_treatment_tab.addWidget(self.treat_sw)
        self.data_treatment_tab.addWidget(self.param_wdgt)
        
        """ Making a menubar """
        self.menu_bar = self.menuBar()
        
        # File menu
        self.file_menu = self.menu_bar.addMenu('File')
        
        self.settings_action = QAction('&Settings', self)
        self.settings_action.setShortcut("Ctrl+I")
        
        self.settings_action.triggered.connect(lambda: os.system(os.path.join(
            os.path.dirname(__file__),
            'data',
            'read_options.json')))
        self.file_menu.addAction(self.settings_action)

        self.quit_action = QAction('&Quit', self)
        self.quit_action.setShortcut("Ctrl+Q")
        self.quit_action.triggered.connect(sys.exit)
        self.file_menu.addAction(self.quit_action)
        
        # Simulation menu
        self.sim_menu = self.menu_bar.addMenu('Simulation')
        
        self.add_sim_w_menu = QAction('&New', self)
        self.add_sim_w_menu.setShortcut("Ctrl+Shift+N")
        self.add_sim_w_menu.triggered.connect(self.edit_simulation_from_list)
        self.sim_menu.addAction(self.add_sim_w_menu)
        
        # About menu
        self.help_menu = self.menu_bar.addMenu('Help')
        
        self.help_lang_menu = self.help_menu.addMenu('Language')
        self.help_lang_actiongrp = QActionGroup(self)
        
        self.help_lang_eng = QAction('English', self)
        self.help_lang_eng.setCheckable(True)
        self.help_lang_eng.setChecked(True)
        self.help_lang_dan = QAction('Dansk', self)
        self.help_lang_dan.setCheckable(True)

        self.help_lang_eng.triggered.connect(self.set_gui_language)
        self.help_lang_dan.triggered.connect(self.set_gui_language)
        
        self.help_lang_menu.addAction(self.help_lang_eng)
        self.help_lang_actiongrp.addAction(self.help_lang_eng)
        self.help_lang_menu.addAction(self.help_lang_dan)
        self.help_lang_actiongrp.addAction(self.help_lang_dan)

        self.help_about_menu = QAction('About', self)
        self.help_about_menu.triggered.connect(self.show_about_dialog)
        self.help_about_menu.setShortcut("F10")
        self.help_menu.addAction(self.help_about_menu)
        
        #Makes "Table of Data" tab
        self.widget_table = DataTableTab(self)
        self.all_the_tabs.addTab(self.widget_table, "Table of Data")

        # Showing the GUI
        self.load_t_tau_data()
        
        self.help_lang_eng.trigger()
        self.showMaximized()
        self.show()
    
    def set_gui_language(self):

        language = self.sender().text()
        self.sample_mass_inp.setToolTip(self.tooltips_dict[language]['m_sample'])
        self.molar_mass_inp.setToolTip(self.tooltips_dict[language]['M_sample'])
        self.sample_xd_inp.setToolTip(self.tooltips_dict[language]['Xd_sample'])
        self.constant_terms_inp.setToolTip(self.tooltips_dict[language]['const_terms'])
        self.var_amount_inp.setToolTip(self.tooltips_dict[language]['var_amounts'])

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
            T = self.meas_temps[i]
            newitem = QListWidgetItem()
            newitem.setText('{:<6.2f} K, {}, {}'.format(round(T,2),True, True))
            plotting_dict = {'temp': self.meas_temps[i],
                             'raw': True,
                             'fit': True}
            newitem.setData(32, plotting_dict)
            t_float = (T-self.Tmin)/(self.Tmax-self.Tmin)
            newitem.setBackground(QColor(to_hex(self.temperature_cmap(t_float))))
            self.treat_raw_fit_list.addItem(newitem)
    
    def set_fit_message_on_statusbar(self, n, N):

        self.statusBar.showMessage('Fitted {}/{} curves'.format(n,N))

    def fit_Xp_Xpp_standalone(self):
        
        try:
            # Check to see if there has been loaded a data frame
            self.raw_df.columns
            # Check to see if the data to work on is in the data frame
            assert 'Xp_m (emu/(Oe*mol))' in self.raw_df.columns
            
        except AttributeError:
            print("There hasn't been loaded a data frame to work on")
        except AssertionError:
            msg = QMessageBox()
            msg.setWindowTitle('Error')
            msg.setText('Calculate diamagnetic correction first\nto make Xp_m and Xpp_m for the algorithm')
            msg.exec_()
        
        else:
            self.statusBar.showMessage('Running fit...')
            
            T = [x for x in self.meas_temps]
            Xs, Xt, tau, alpha, resid, tau_fit_err = [],[],[],[],[],[]
            
            v_all = [np.array(self.temp_subsets[t_idx]['AC Frequency (Hz)']) for t_idx in range(self.num_meas_temps)]
            Xp_all = [np.array(self.temp_subsets[t_idx]['Xp_m (emu/(Oe*mol))']) for t_idx in range(self.num_meas_temps)]
            Xpp_all = [np.array(self.temp_subsets[t_idx]['Xpp_m (emu/(Oe*mol))']) for t_idx in range(self.num_meas_temps)]
            
            inputs = tuple(zip(v_all, Xp_all, Xpp_all))
            
            #thread = QThread()
            #worker = FitWorker(inputs)
            #worker.moveToThread(thread)
            #thread.started.connect(worker.run)
            #worker.finished.connect(thread.quit)
            #worker.finished.connect(worker.deleteLater)
            #thread.finished.connect(thread.deleteLater)
            #worker.progress.connect(self.set_fit_message_on_statusbar)

            #thread.start()
            #res = worker.res

            with Pool() as pool:
                res = pool.starmap(fit_Xp_Xpp_genDebye, inputs)
            
            tau = [e[0] for e in res]
            tau_fit_err = [e[1] for e in res]
            alpha = [e[2] for e in res]
            Xs = [e[3] for e in res]
            Xt = [e[4] for e in res]
            resid = [e[5] for e in res]
            
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
    
    def load_sample_data(self):
    
        filename_info = QFileDialog().getOpenFileName(self, 'Open file', self.last_loaded_file)
        filename = filename_info[0]

        try:
            f = open(filename, 'r')
            d = f.readlines()
            f.close()
            
            assert all([len(line.split())>=2 for line in d])
            
        except FileNotFoundError:
            print('File was not selected')
        except UnicodeDecodeError:
            print('Cant open a binary file')
        except AssertionError:
            print('Some of the lines have lengths less than two')
        else:
            
            # These are the default values that are "read" if nothing else is
            # seen in the file
            m_sample = '0'
            M_sample = '0'
            Xd_sample = '0'
            constant_terms = '0'
            var_amount = '0,0'
            
            self.last_loaded_file = os.path.split(filename)[0]
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
                                                   self.last_loaded_file)
        filename = filename_info[0]
        
        try:
            assert filename != ''
            self.last_loaded_file = os.path.split(filename)[0]
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
        
    def make_diamag_correction_calculation(self):
        
        if self.raw_df is None:
            # Don't do the calculation, if there is nothing to calculate on
            pass
        
        else:
            try:
                m_sample = float(self.sample_mass_inp.text())
                M_sample = float(self.molar_mass_inp.text())
                Xd_sample = float(self.sample_xd_inp.text())
                constant_terms = [float(x) for x in self.constant_terms_inp.text().split(',')]
                var_am = [float(x) for x in self.var_amount_inp.text().split(',')]
                
                assert len(var_am)%2==0
                paired_terms = [(var_am[n], var_am[n+1]) for n in range(0,len(var_am),2)]
                
                if Xd_sample == 0:
                    Xd_sample = -6e-7*M_sample
                
            except (ValueError, AssertionError):
                msg = QMessageBox()
                msg.setWindowTitle('Error')
                msg.setText('Something wrong in "Sample information"\n')
                msg.exec_()
            
            else:
                H = self.raw_df['AC Amplitude (Oe)']
                H0 = self.raw_df['Magnetic Field (Oe)']
                Mp = self.raw_df["Mp (emu)"]
                Mpp = self.raw_df["Mpp (emu)"]
                
                # Get molar, corrected values from function in process_ac
                Mp_molar, Mpp_molar, Xp_molar, Xpp_molar = diamag_correction(
                    H, H0, Mp, Mpp, m_sample, M_sample, Xd_sample,
                    constant_terms = constant_terms, paired_terms = paired_terms)
            
                # PUT THE CODE HERE TO INSERT CORRECTED VALUES INTO DATA FRAME
                if "Mp_m (emu/mol)" in self.raw_df.columns:
                    self.raw_df.replace(to_replace="Mp_m (emu/mol)", value=Mp_molar)
                    self.raw_df.replace(to_replace="Mpp_m (emu/mol)", value=Mpp_molar)
                    self.raw_df.replace(to_replace="Xp_m (emu/(Oe*mol))", value=Xp_molar)
                    self.raw_df.replace(to_replace="Xpp_m (emu/(Oe*mol))", value=Xpp_molar)
                else:
                    Mp_idx = self.raw_df.columns.get_loc('Mp (emu)')
                    self.raw_df.insert(Mp_idx+1, column="Mp_m (emu/mol)", value=Mp_molar)
                    
                    Mpp_idx = self.raw_df.columns.get_loc('Mpp (emu)')
                    self.raw_df.insert(Mpp_idx+1, column="Mpp_m (emu/mol)", value=Mpp_molar)
                    
                    Xp_idx = self.raw_df.columns.get_loc('Xp (emu/Oe)')
                    self.raw_df.insert(Xp_idx+1, column="Xp_m (emu/(Oe*mol))", value=Xp_molar)
                    
                    Xpp_idx = self.raw_df.columns.get_loc('Xpp (emu/Oe)')
                    self.raw_df.insert(Xpp_idx+1, column="Xpp_m (emu/(Oe*mol))", value=Xpp_molar)
            
            self.update_temp_subsets()
            self.update_analysis_combos()
            self.widget_table.updatetable()

            w = MagMessage('Diamagnetic correction',
                          'Diamagnetic correction successful!').exec_()
    
    def update_itemdict(self, item, itemdict):
        
        item.setData(32, itemdict)
        item.setText('{:<6.2f} K, {}, {}'.format(round(itemdict['temp'],2),
                                         itemdict['raw'],
                                         itemdict['fit']))
    
    def update_raw_plot(self):
        
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
        
        self.treat_fit_plot.ax.clear()
        plot_type = self.fit_data_plot_combo.currentText()
        
        if plot_type == 'FreqVSXp':
            x_name = 'AC Frequency (Hz)'
            y_name = 'Xp_m (emu/(Oe*mol))'
            fcn_y = Xp_
            x_scale = 'log'
        elif plot_type == 'FreqVSXpp':
            x_name = 'AC Frequency (Hz)'
            y_name = 'Xpp_m (emu/(Oe*mol))'
            fcn_y = Xpp_
            x_scale = 'log'
        elif plot_type == 'ColeCole':
            x_name = 'Xp_m (emu/(Oe*mol))'
            y_name = 'Xpp_m (emu/(Oe*mol))'
            fcn_y = Xpp_
            x_scale = 'linear'
            
        for row in range(self.num_meas_temps):
        
            T = self.meas_temps[row]
            rgb = self.temperature_cmap((T-self.Tmin)/(self.Tmax-self.Tmin))
            markercolor = 'k'
            if self.fit_data_color_cb.isChecked():
                markercolor = rgb
            
            if plot_type == 'ColeCole':
                x_data = Xp_(self.temp_subsets[row]['AC Frequency (Hz)'],
                             self.raw_data_fit['ChiS'].iloc[row],
                             self.raw_data_fit['ChiT'].iloc[row],
                             self.raw_data_fit['Tau'].iloc[row],
                             self.raw_data_fit['Alpha'].iloc[row])
            else:
                x_data = self.temp_subsets[row][x_name]
                
            item = self.treat_raw_fit_list.item(row)
            itemdict = item.data(32)
            
            if itemdict['raw']:
                self.treat_fit_plot.ax.plot(self.temp_subsets[row][x_name],
                                            self.temp_subsets[row][y_name],
                                            marker='o',
                                            mec=markercolor,
                                            mfc='none',
                                            linestyle='None')
            if itemdict['fit']:
                self.treat_fit_plot.ax.plot(x_data,
                                            fcn_y(self.temp_subsets[row]['AC Frequency (Hz)'],
                                                  self.raw_data_fit['ChiS'].iloc[row],
                                                  self.raw_data_fit['ChiT'].iloc[row],
                                                  self.raw_data_fit['Tau'].iloc[row],
                                                  self.raw_data_fit['Alpha'].iloc[row]),
                                            c=rgb)
            
        self.treat_fit_plot.ax.set_xscale(x_scale)
        self.treat_fit_plot.ax.set_xlabel(x_name)
        self.treat_fit_plot.ax.set_ylabel(y_name)
        
        norm = mpl.colors.Normalize(vmin=self.Tmin, vmax=self.Tmax)
        self.treat_fit_plot.fig.colorbar(
            mpl.cm.ScalarMappable(norm=norm,
                                  cmap=self.temperature_cmap),
                                        orientation='horizontal',
            cax=self.treat_fit_plot.cax)
        
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
            self.raw_df.insert(Xp_idx+1, column='Mpp (emu)', value=Mpp)
            
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
            self.current_file = filename
            
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
            self.treat_fit_plot.cax.clear()
            
            combo_idx = self.analysis_plot_type_combo.findText('Raw data')
            self.analysis_plot_type_combo.setCurrentIndex(combo_idx)
            
            # Updating analysis combos, which will automatically draw the new data
            self.update_analysis_combos()
            #Updates "Table of Data" tab with the loaded data
            self.widget_table.updatetable()

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
        idx_list = [0]
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
            self.temp_subsets.append(self.raw_df.iloc[idx_list[n]:idx_list[n+1]])
    
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
    
    def show_fitted_params(self):
        
        try:
            fit = self.fit_history[0]
        except (IndexError, TypeError):
            w = MagMessage('Fit history error', 'There is no fit history yet!')
        else:
            w = ParamDialog(self, self.fit_history)
        finally:
            w.exec_()
        
    def edit_simulation_from_list(self):
    
        try:
            sender = self.sender().text()
        except AttributeError:
            # Sent here because of double-click on QListWidget
            action = 'Edit'
        else:
            if sender == 'Edit':
                action = 'Edit'
            elif sender in ('New', '&New'):
                action = 'New'

        if action == 'Edit':
            try:
                sim_item = self.list_of_simulations.selectedItems()[0]
            except IndexError:
                w = MagMessage("Did not find any selected line",
                               "Select a line first to edit it")
                w.exec_()
                return
            else:
                data = sim_item.data(32)
                params = data['params']
                T_vals = data['T_vals']
                line = data['line']
                color = line._color
                label = line._label
        
        elif action == 'New':
            
            params = default_parameters()
            T_vals = [1,3]
            line = False
            label = None
            color = None
            if len(self.simulation_colors)<1:
                self.statusBar.showMessage("ERROR: can't make any more simulations")
                return
            
        sim_dialog = SimulationDialog(parent=self,
                                      fit_history=self.fit_history,
                                      params = params,
                                      min_max_T=T_vals)
        finished_value = sim_dialog.exec_()
        functions = [bool(sim_dialog.params[p].value)
                     for p in sim_dialog.params if 'use' in p]
        
        try:
            assert finished_value
            assert(any(functions))
        except AssertionError:
            pass
        else:
            params = sim_dialog.params
            T_vals = sim_dialog.min_max_T
            
            if line:
                self.ana_plot.ax.lines.remove(line)
            else:
                # In this case, there was no old line and therefore also no sim_item
                """https://stackoverflow.com/questions/55145390/pyqt5-qlistwidget-with-checkboxes-and-drag-and-drop"""
                sim_item = QListWidgetItem()
                sim_item.setFlags( sim_item.flags() | Qt.ItemIsUserCheckable )
                sim_item.setCheckState(Qt.Checked)
                
                self.list_of_simulations.addItem(sim_item)
                color = self.simulation_colors.pop()
                label = names.get_first_name()
            
            line = add_partial_model(self.ana_plot.fig,
                                     T_vals[0],
                                     T_vals[1],
                                     params,
                                     c=color,
                                     label=label)
            
            list_item_data = {'params': params,
                              'T_vals': T_vals,
                              'line': line,
                              'color': color}
            
            new_item_text = self.represent_simulation(T_vals, params)
            
            sim_item.setData(32, list_item_data)
            sim_item.setText(new_item_text)
            sim_item.setBackground(QColor(to_hex(color)))
            
            self.redraw_simulation_lines()

    def represent_simulation(self, T_vals, params):
        
        fs = [p for p in params if 'use' in p]
        used = [bool(params[p].value) for p in fs]
        text = ['Using QT: {}, Raman: {}, Orbach: {}\n'.format(*used),
                'to plot between {} K and {} K\n'.format(*T_vals)]
        return ''.join(text)

    def delete_sim(self):
        
        try:
            sim_item = self.list_of_simulations.selectedItems()[0]
        except IndexError:
            pass
        else:
            line_pointer = sim_item.data(32)['line']
            line_color = line_pointer._color
            
            self.ana_plot.ax.lines.remove(line_pointer)
            self.ana_plot.canvas.draw()
            
            item_row = self.list_of_simulations.row(sim_item)
            sim_item = self.list_of_simulations.takeItem(item_row)
            
            self.simulation_colors.append(line_color)
            
            del sim_item
    
    def redraw_simulation_lines(self):
        
        for idx in range(self.list_of_simulations.count()):
            item = self.list_of_simulations.item(idx)
            data = item.data(32)
            
            if item.checkState() == Qt.Checked:
                data['line']._visible = True
            elif item.checkState() == Qt.Unchecked:
                data['line']._visible = False
        
        self.ana_plot.canvas.draw()
        
    def plot_t_tau_on_axes(self):
        
        if self.plotted_data_pointers is not None:
            for line in self.plotted_data_pointers:
                line.remove()
        self.plotted_data_pointers = []
        
        if self.data_dtau is None:
            used, = self.ana_plot.ax.plot(1/self.used_T, np.log(self.used_tau), 'bo', zorder=0.1)
            not_used, = self.ana_plot.ax.plot(1/self.not_used_T, np.log(self.not_used_tau), 'ro', zorder=0.1)
            self.plotted_data_pointers.append(used)
            self.plotted_data_pointers.append(not_used)
        else:
            err_used_point, caplines1, barlinecols1 = self.ana_plot.ax.errorbar(1/self.used_T,
                                                                                np.log(self.used_tau),
                                                                                yerr=self.used_dtau,
                                                                                fmt='bo',
                                                                                ecolor='b',
                                                                                label='Data',
                                                                                zorder=0.1)
            err_not_used_point, caplines2, barlinecols2 = self.ana_plot.ax.errorbar(1/self.not_used_T,
                                                                                    np.log(self.not_used_tau),
                                                                                    yerr=self.not_used_dtau,
                                                                                    fmt='ro',
                                                                                    ecolor='r',
                                                                                    label='Data',
                                                                                    zorder=0.1)

            self.plotted_data_pointers.append(err_used_point)
            self.plotted_data_pointers.append(err_not_used_point)
            for e in [caplines1, caplines2, barlinecols1, barlinecols2]:
                for line in e:
                    self.plotted_data_pointers.append(line)
        
        self.ana_plot.canvas.draw()
    
    def reset_analysis_containers(self):
        self.data_T = None
        self.data_tau = None
        self.data_dtau = None
        
        self.used_T = None
        self.not_used_T = None
        self.used_tau = None
        self.not_used_tau = None
        self.used_dtau = None
        self.not_used_dtau = None
        
        self.used_indices = None
        
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
        
        try:
            D = np.loadtxt(filename,
                           skiprows=1,
                           delimiter=';')
        except (ValueError, OSError) as error_type:
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
        If three columns in D: assume the third is dtau
        If four columns in D: assume third is alpha, fourth is tau_fit_error
            dtau will then be calculated from these values
        """
        
        T = D[:,0]
        tau = D[:,1]
        
        sort_indices = T.argsort()
        self.data_T = T[sort_indices]
        self.data_tau = tau[sort_indices]
        self.data_dtau = None
        
        if D.shape[1]==3:
            # Three columns in the array loaded, assume the third is error
            dtau = D[:,2]
            dtau = dtau[sort_indices]
            
        elif D.shape[1]==4:
            # Four columns in the array loaded, assume the third is alpha
            # and that the fourth is the fitting error on tau
            alpha = D[:,2]
            tau_fit_err = D[:,3]
            dtau = tau_err_RC(tau, tau_fit_err, alpha)
            dtau = dtau[sort_indices]
        else:
            dtau = None
            
        self.data_dtau = dtau
        
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
    
    def make_fitting_function(self, perform_this_fit):
    
        use_QT = 0
        use_R = 0
        use_O = 0
        
        if 'QT' in perform_this_fit: use_QT = 1
        if 'R' in perform_this_fit: use_R = 1
        if 'O' in perform_this_fit: use_O = 1
        
        fn = lambda T, tQT, Cr, n, t0, Ueff: np.log(1/(use_QT*1/_QT(T, tQT) + use_R*1/_R(T, Cr, n) + use_O*1/_O(T, t0, Ueff)))
        
        return fn
        
    def fit_relaxation(self, guess, perform_this_fit):
        """
        First function used to fit t-tau-data, but this function uses curve_fit,
        where it is more difficult to control parameters. Trying to switch to
        lmfit, where control of refinable parameters is more extensive.
        """
        
        f = getFittingFunction(fitType=perform_this_fit)
        p0 = getStartParams(guess, fitType=perform_this_fit)
        
        if self.used_dtau is None:
            popt, pcov = curve_fit(f, self.used_T, np.log(self.used_tau), p0)
        else:
            popt, pcov = curve_fit(f, self.used_T, np.log(self.used_tau), p0, sigma=np.log(self.used_dtau))
        
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
        window_title = 'Fit aborted'
        msg_text = ''
        msg_details = ''
            
        try:
            
            # This will raise TypeError and IndexError first
            # to warn that no data was loaded
            fitwith = self.read_fit_type_cbs()
            assert fitwith != ''
            guess = getParameterGuesses(self.used_T, self.used_tau, fitwith)
            
            Tmin = self.temp_line[1].value()
            Tmax = self.temp_line[3].value()
            assert Tmin != Tmax
            
            guess_dialog = GuessDialog(self,
                                       guess,
                                       fitwith)
            accepted = guess_dialog.exec_()
            if not accepted: raise NoGuessExistsError
            
            # If both fit and temperature setting are good,
            # and the GuessDialog was accepted, get the
            # guess and perform fitting
            
            params = guess_dialog.current_guess
            minimize_res = fit_relaxation(self.used_T, self.used_tau, params)
            
        except (AssertionError, IndexError):
            msg_text = 'Bad temperature or fit settings'
            msg_details = """Possible errors:
 - min and max temperatures are the same
 - no fit options have been selected
 - can't fit only one data point"""
        except RuntimeError:
            msg_text = 'This fit cannot be made within the set temperatures'
        except ValueError as e:
            msg_text = 'No file has been loaded'
        except TypeError as e:
            msg_text = 'No data has been selected'
            print(e)
        except NoGuessExistsError:
            msg_text = 'Made no guess for initial parameters'
        
        else:
            window_title = 'Fit successful!'
            msg_text = 'Congratulations'
            
            self.add_to_history(minimize_res, fitwith)
            
        finally:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle(window_title)
            msg.setText(msg_text)
            msg.setDetailedText(msg_details)
            msg.exec_()        
    
    def add_to_history(self, p_fit, perform_this_fit):
        
        now = datetime.datetime.now()
        now = now.strftime("%d.%m %H:%M:%S")

        if len(self.fit_history)>9:
            self.fit_history.pop()
        self.fit_history.insert(0, (perform_this_fit, p_fit, now))
        