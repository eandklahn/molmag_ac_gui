#std packages
from importlib.resources import read_text
import sys
import os
import json

#third-party packages
from PyQt5.QtGui import QIcon, QFont
from PyQt5.QtWidgets import (QMainWindow, QAction, QTabWidget, QStatusBar,
                             QActionGroup)
from matplotlib.colors import LinearSegmentedColormap

#local imports
from .__init__ import __version__
from .dialogs import  AboutDialog
from .Datatabletab import Datatabletab 
from .DataAnalysisTab import DataAnalysisTab
from .DataTreatmentTab import DataTreatmentTab
from . import data as pkg_static_data


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
                                  'https://eandklahn.github.io', 
                                  'Contributor': 
                                  'Sofie Stampe Leiszner (sofiesl@chem.au.dk)'
                                  }
        
        self.last_loaded_file = os.getcwd() #Remember the last used folder.
        self.current_file = ''


        """ Things to do with how the window is shown """
        self.setWindowTitle('Molmag AC GUI v{}'.format(__version__))
        self.setWindowIcon(QIcon('double_well_potential_R6p_icon.ico'))
        
        """ FONT STUFF """
        self.headline_font = QFont() #Defines headline font
        self.headline_font.setBold(True) #Makes it bold
        
        # Data containers for treatment
        self.read_options = json.loads(read_text(pkg_static_data,
                                                'read_options.json'))
        
        self.diamag_constants = json.loads(read_text(pkg_static_data,
                                                    'diamag_constants.json'))
        
        self.temperature_cmap = LinearSegmentedColormap.from_list(
            'temp_colormap',
            json.loads(read_text(pkg_static_data, 'default_colormap.json')))
        
        self.tooltips_dict = json.loads(read_text(pkg_static_data,
                                                  'tooltips.json'))

        """ Setting up the main tab widget """
        self.all_the_tabs = QTabWidget() 
        self.setCentralWidget(self.all_the_tabs)
        self.statusBar = QStatusBar() #A statusbar in the buttom of the window
        self.setStatusBar(self.statusBar)
        #self.setStyleSheet(read_text(pkg_static_data, 'styles.qss')
        
        """ Adding all the tabs"""    



        #Makes "Data treatment" tab    
        self.data_treat = DataTreatmentTab(self)        
        self.all_the_tabs.addTab(self.data_treat, "Data treatment")

        #Makes "Table of Data" tab
        self.widget_table = Datatabletab(self)
        self.all_the_tabs.addTab(self.widget_table, "Table of data")       
        
        #Makes "Data analysis" tab
        self.data_ana = DataAnalysisTab(self)        
        self.all_the_tabs.addTab(self.data_ana, "Data analysis")   
        
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
        self.add_sim_w_menu.triggered.connect(self.data_ana.edit_simulation_from_list)
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

        # Showing the GUI
        self.data_ana.load_t_tau_data()
        self.help_lang_eng.trigger()        
        self.showMaximized()
        self.show()

    def set_gui_language(self): 
        self.gui_language = self.sender().text()

    def show_about_dialog(self):
        w = AboutDialog(info=self.about_information)
        w.exec_()
    

    

        


