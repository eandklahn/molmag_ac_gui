#std packages
import sys
import os
import json
from importlib.resources import read_text
import datetime

#third-party packages
from PyQt5.QtGui import QIcon, QFont
from PyQt5.QtWidgets import (QMainWindow, QAction, QTabWidget, QStatusBar,
                             QActionGroup)
from matplotlib.colors import LinearSegmentedColormap

#local imports
from .__init__ import __version__
from .dialogs import  AboutDialog
from .DataTableTab import DataTableTab 
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
        
        # Initilization: 
        self.set_about_and_current_file() 
        self.set_data_containers() 

        # Setting up GUI: 
        self.setWindowTitle('Molmag AC GUI v{}'.format(__version__))
        self.setWindowIcon(QIcon('double_well_potential_R6p_icon.ico'))
        self.make_tab_wdgt() 
        self.make_menubar()        
        
        # Finalizing GUI
        self.help_lang_eng.trigger()        
        self.showMaximized()
        self.show()



    def set_about_and_current_file(self):
        """ Sets the about information, remembers the last loaded file and sets current file to an empty string """
        
        self.about_information = {'author':
                                  '\nEmil A. Klahn (eklahn@chem.au.dk) \nSofie Stampe Leiszner (sofiesl@chem.au.dk)',
                                  'webpage':
                                  'https://chem.au.dk/en/research/research-areas-and-research-groups/inorganicchemistrymaterialschemistry/molecular-magnetism',
                                  'personal':
                                  'https://github.com/eandklahn/molmag_ac_gui'
                                  }
        
        self.last_loaded_file = os.getcwd() #Remember the last used folder.
        self.current_file = ''

    def set_data_containers(self): 
        """Sets data containers (read options, diamagnetic constants, temperature colormap 
        and the tooltips dictionary) """

        self.read_options = json.loads(read_text(pkg_static_data,
                                                'read_options.json'))
        
        self.diamag_constants = json.loads(read_text(pkg_static_data,
                                                    'diamag_constants.json'))
        
        self.temperature_cmap = LinearSegmentedColormap.from_list(
            'temp_colormap',
            json.loads(read_text(pkg_static_data, 'default_colormap.json')))
        
        self.tooltips_dict = json.loads(read_text(pkg_static_data,
                                                  'tooltips.json'))

    def make_tab_wdgt(self): 
        """ Make the main window with the three tabs: Data treatment, Table of Data
        and Data analysis. """

        # Make the overall tab widget that contain all tabs 
        self.all_the_tabs = QTabWidget() 
        self.setCentralWidget(self.all_the_tabs)

        # Make a statusbar in the buttom of the window
        self.statusBar = QStatusBar() 
        self.setStatusBar(self.statusBar)

        # Make "Data treatment", "Table of data" and "Data analysis" tabs
        self.widget_table = DataTableTab(self)
        self.data_treat = DataTreatmentTab(self)    
        self.data_ana = DataAnalysisTab(self)        

        #Adds the tabs to the tabwidget 
        self.all_the_tabs.addTab(self.data_treat, "Data treatment")
        self.all_the_tabs.addTab(self.widget_table, "Table of data")       
        self.all_the_tabs.addTab(self.data_ana, "Data analysis (AC only)")   


    def add_file_menu(self): 
        """ Adds the file menu in the top menu with options settings and quit their 
        corresponding shortcuts """

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


    def add_simulation_menu(self): 
        """Adds the simulation option to the top menubar. This is not fully functional yet. """

        self.sim_menu = self.menu_bar.addMenu('Simulation')
        
        self.add_sim_w_menu = QAction('&New', self)
        self.add_sim_w_menu.setShortcut("Ctrl+Shift+N")
        self.add_sim_w_menu.triggered.connect(self.data_ana.edit_simulation_from_list)
        self.sim_menu.addAction(self.add_sim_w_menu)

    def add_help_menu(self): 
        """ Adds a help menu, where the language can be changed (this only works for the 
        tooltips in the sample information menu, not for the rest of the GUI). 
        Also adds an option to view the about information. """
        
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

    def make_menubar(self): 
        """ Makes the menubar in the top of the GUI window"""

        self.menu_bar = self.menuBar()
        self.add_file_menu() 
        self.add_simulation_menu() 
        self.add_help_menu()
        

    def set_gui_language(self): 
        """ Sets the GUI language """

        self.gui_language = self.sender().text()

    def show_about_dialog(self):
        """ Shows the about dialog"""

        w = AboutDialog(info=self.about_information)
        w.exec_()
    


        
