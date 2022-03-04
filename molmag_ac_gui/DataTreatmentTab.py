#std packages
import os
from multiprocessing import Pool

#third-party packages 
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import to_hex
import matplotlib.ticker as mticker
from PyQt5.QtWidgets import (QFileDialog, QFormLayout, QListWidgetItem, QWidget, 
                             QVBoxLayout, QHBoxLayout, QComboBox, QStackedWidget, 
                             QCheckBox, QListWidget, QSplitter, QPushButton, QLabel)
from PyQt5.QtGui import QColor

from molmag_ac_gui.layout import make_btn, make_headline, make_line, headline_font

from math import isclose

#local imports
from .dialogs import PlottingWindow, MagMessage, FitResultPlotStatus, SampleInformation
from .utility import read_ppms_file, update_data_names, formatlabel
from .exceptions import FileFormatError
from .process_ac import (diamag_correction, diamag_correction_dc, fit_Xp_Xpp_genDebye, tau_err_RC, 
                         Xp_, Xpp_)


class DataTreatmentTab(QSplitter): 
    def __init__(self, parent): 
        super(DataTreatmentTab,self).__init__() #Før __init__() peger på QWidget
        self.parent = parent
        self.initUI() 
    
    def initUI(self): 
        #Creates dataframe with header and raw_df
        self.initialize_attributes()

        ## Making the left column: data loading, fitting, visualization controls and show/hide option
        self.layout = QVBoxLayout()
        self.options_wdgt = QWidget()
        
        # Constructing data import and treatment buttons
        make_headline(self, "Data import and treatment", self.layout)        
        self.add_data_treatment_buttons() 

        # Constructing axis content combobox
        make_headline(self, "Axis content", self.layout)
        self.add_plot_type_combobox()

        # Constructing the x and y comboboxes
        make_headline(self, "Raw data: Choose x- and y-contents", self.layout)
        self.add_xy_comboboxes() 

        # Constructing a combobox for plotting fitted data
        make_headline(self, "Fitted: Choose plot type (AC only)", self.layout)
        self.add_fit_combobox() 

        # Constructing parameter controls
        self.add_parameter_controls() 
        
        # Finalizing layout on the left side of the GUI
        self.layout.addStretch()
        self.options_wdgt.setLayout(self.layout)
        self.addWidget(self.options_wdgt)

        ## Making the right columns: data visualization stacking widget  
        self.add_stackwidget()

        ## Finalization 
        self.setSizes([1,1000])
        self.show()

    """ Functions for setting up GUI layout"""

    def initialize_attributes(self): 
        """Initializes attributes such as the raw_df containing all data, temperatures etc."""

        self.raw_df = None
        self.raw_df_header = None
        self.num_meas_freqs = 0
        self.num_meas_temps = 0
        self.temp_subsets = []
        self.meas_temps = []
        self.Tmin, self.Tmax = 0,0
        self.raw_data_fit = None
        self.data_type = None

    def add_plot_type_combobox(self): 
        """ Adds a plot_type combobox such that it is possible to switch between viewing 
        the raw data, fitted data and the 3D plot in the plotting window on the right """

        self.plot_type_combo = QComboBox()
        self.plot_type_combo.currentIndexChanged.connect(self.switch_view)
        self.layout.addWidget(self.plot_type_combo)

    def add_fit_combobox(self): 
        """ Adds a combobox to choose which fitted data plot you want to view.
        Also adds an option to use temperature colors for datapoints."""
        
        self.fit_combo = QComboBox()
        #self.fit_combo.addItems(['Cole-Cole', 'Freq VS Xp', 'Freq VS Xpp'])
        self.fit_combo.currentIndexChanged.connect(self.plot_from_itemlist)
        self.layout.addWidget(self.fit_combo)
        
        self.fit_color_cb_lo = QHBoxLayout()


    def add_data_treatment_buttons(self): 
        """ Adds the first 6 buttons in the option widget on the left side of the gui window. 
        and connects these to the corresponding functions. """

        make_btn(self, "(1) Load datafile", self.load_data, self.layout)
        make_btn(self, "(2) Make diamagnetic correction", self.update_sample_info_and_make_dia_correction, self.layout)

        self.fit_X_btn = QPushButton("(3) Fit X' and X'' (AC only)")
        self.fit_X_btn.clicked.connect(self.fit_Xp_Xpp_standalone)
        self.fit_X_btn.setEnabled(False)
        self.layout.addWidget(self.fit_X_btn)
        
        #Copy fit to data analysis and save fit to file layout
        self.copy_save_layout = QHBoxLayout()

        self.copy_fit_btn = QPushButton("Copy fit to Data analysis (AC only)")
        self.copy_fit_btn.clicked.connect(self.copy_fit_to_analysis)
        self.copy_fit_btn.setEnabled(False)
        self.copy_save_layout.addWidget(self.copy_fit_btn)

        self.save_fit_X_btn = QPushButton("Save fit to file (AC only)")
        self.save_fit_X_btn.clicked.connect(self.save_fit_to_file)
        self.save_fit_X_btn.setEnabled(False)
        self.copy_save_layout.addWidget(self.save_fit_X_btn)

        self.layout.addLayout(self.copy_save_layout)

    def add_xy_comboboxes(self): 
        """Adds an option to show more option for x and y when plotting raw data. 
        More options will give you all the non-empty datatypes from your datafile. """

        self.xy_lo = QFormLayout() 

        self.x_combo = QComboBox()
        self.x_combo.currentIndexChanged.connect(self.plot_raw_data)        
        self.y_combo = QComboBox()
        self.y_combo.currentIndexChanged.connect(self.plot_raw_data)
        
        self.xy_lo.addRow("x: ", self.x_combo) 
        self.xy_lo.addRow("y: ", self.y_combo) 

        self.layout.addLayout(self.xy_lo)

        self.xy_options_lo = QHBoxLayout()
        make_line(self, "Show additional options for x and y ", self.xy_options_lo)
        self.xy_options_cb = QCheckBox()
        self.xy_options_cb.stateChanged.connect(self.update_xy_combos)
        self.xy_options_lo.addWidget(self.xy_options_cb)

        self.layout.addLayout(self.xy_options_lo)

    def add_stackwidget(self): 
        """ Adds the stackwidget which is the plotting window on the right side of the gui window"""

        self.sw = QStackedWidget()
        
        self.raw_plot = PlottingWindow()
        self.fit_plot = PlottingWindow(make_ax="cax")
        self.threeD_plot = PlottingWindow(make_ax = "z") 
        self.HM_plot = PlottingWindow(make_ax="cax")
        self.TXT_plot = PlottingWindow()

        self.sw.addWidget(self.raw_plot)
        self.sw.addWidget(self.fit_plot)
        self.sw.addWidget(self.threeD_plot)
        self.sw.addWidget(self.HM_plot)
        self.sw.addWidget(self.TXT_plot)

        self.addWidget(self.sw)


    def add_parameter_controls(self): 
        """ Adds controls for the data plotting where one can pick which datapoints
        and fits to be visualized. """
        
        self.headline = QLabel( "Hide/show raw data and lines")
        self.headline.setFont(headline_font)
        self.layout.addWidget(self.headline)
        
        make_line(self, "Double click to edit list", self.layout)
        
        #List of fitted raw data
        self.raw_fit_list = QListWidget()
        self.layout.addWidget(self.raw_fit_list)
        self.raw_fit_list.doubleClicked.connect(self.update_fitted_and_3D_plot)
        
    """ Other functions"""
    
    def fill_df_data_values(self):
        """ If Xp is in the loaded data, and Mp is not: Mp will be calculated and added to the
        dataframe. If Mp is in the loaded data and Xp is not, Xp will be calculated and added 
        to the dataframe. """

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

    def cleanup_loaded_ppms(self):
        """Cleans up the data that is loaded in by removing all empty columns,
        removes Comment column and renames the columns if necessary"""

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

    def update_temp_subsets(self):
        """ Splits datasets into one dataset for each temperature. 
        The dataset it splitted based on when the frequency "restarts". 
        The function assumes that the frequency is always increasing within a measurement. """

        self.temp_subsets = []
        idx_list = [0]
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

    def update_meas_temps(self):
        """ Updates the measured temperatures based on the temperatures on self.temp_subsets.
        Does this by taking the average of all measured temperatures in that given subset.
        Also sets a minimum and maximum temperature used for colormapping """

        meas_temps = []
        for sub in self.temp_subsets:
            meas_temps.append(sub['Temperature (K)'].mean())
        
        self.meas_temps = np.array(meas_temps)
        self.num_meas_temps = len(self.meas_temps)
        self.Tmin = self.meas_temps.min()
        self.Tmax = self.meas_temps.max()



    def update_xy_combos(self):
        """ Updates the x and y comboboxes for raw data plotting such that they contain only the desired datatypes. 
        If the xy_options checkbox is not chekced, fewer datatypes are possible to choose from in the combobox. 
        If the xy_options checkbox is checked, all datatypes in the dataframe are possible to choose from the combobox. """

        self.x_combo.clear()
        self.y_combo.clear()
        if self.data_type == "AC":
            chosenlabels = ["Temperature (K)", "AC Frequency (Hz)", "AC Amplitude (Oe)", "Magnetic Field (Oe)", "Mp (emu)",\
                           "Mpp (emu)", "Xp (emu/Oe)", "Xpp (emu/Oe)"]
            molarlabels = ["Mp_m (emu/mol)", "Mpp_m (emu/mol)", "Xp_m (emu/(Oe*mol))", "Xpp_m (emu/(Oe*mol))"]
        elif self.data_type == "DC": 
            chosenlabels = ["Time Stamp (sec)", "Temperature (K)", "Moment (emu)","Magnetic Field (Oe)"]
            molarlabels = ["Moment_m (emu/mol)", "X_m (emu/(Oe*mol))", "XT_m (emu*K/Oe)"]
        
        if not self.xy_options_cb.isChecked():
            if all([label in self.raw_df for label in chosenlabels + molarlabels]): #if molar properties have been calculated
                self.x_combo.addItems(chosenlabels + molarlabels)
                self.y_combo.addItems(chosenlabels + molarlabels)
            elif all([label in self.raw_df for label in chosenlabels]): #If no molar properties has been calculated
                self.x_combo.addItems(chosenlabels)
                self.y_combo.addItems(chosenlabels)
            else: 
                self.y_combo.addItems(self.raw_df.columns)
                self.x_combo.addItems(self.raw_df.columns)
        else:  #If all properties in chosen labels are not there, show all labels from input file 
            self.y_combo.addItems(self.raw_df.columns)
            self.x_combo.addItems(self.raw_df.columns)


    def determine_data_type(self): 
        try: #Checks if it is AC data
            self.raw_df['AC Amplitude (Oe)']
            self.data_type = "AC"

        except: #Otherwise expects DC data 
            self.data_type = "DC"


    def reset_btns(self):
        self.fit_X_btn.setEnabled(False)
        self.copy_fit_btn.setEnabled(False)
        self.save_fit_X_btn.setEnabled(False)

    def load_data(self):
        """ This function uses subfunctions to perform the following tasks:  
        Loads the ppms data, clears the dataframe, fills in the dataframe with the loaded raw data, 
        adds Mp or Xp to the dataframe, cleans up the dataframe, updates temperatures subsets, 
        clears the raw data and fitted plots for any old data, plots raw data and updates the data table. """

        open_file_dialog = QFileDialog()
        filename_info = open_file_dialog.getOpenFileName(self, 'Open file', self.parent.last_loaded_file)
        #filename_info =  ('C:/Users/au592011/OneDrive - Aarhus Universitet/Skrivebord/TestData_MAG/dc-data/DC_DyBrTHF_ODPM2blank_1.dat', 'All Files (*)')
        #filename_info =  ("C:/Users/au592011/OneDrive - Aarhus Universitet/Skrivebord/TestData_MAG/ac-data/ac-data/dy-dbm/20180209DyII_1000.dat", 'All Files (*)')
        filename = filename_info[0]
        try:
            # FileNotFoundError and UnicodeDecodeError will be raised here
            potential_header, potential_df = read_ppms_file(filename)
            if potential_header is None:
                raise FileFormatError(filename)
            summary = update_data_names(potential_df, self.parent.read_options)
            counts = [val>1 for key, val in summary.items()]

        except FileNotFoundError:
            # Did not read any file
            pass
        except UnicodeDecodeError:
            # File being read is binary, not a text file
            MagMessage("Unicode error", "The file is not a text file").exec_() 
        except FileFormatError as e:
            # File does not have correct header and data blocks
            MagMessage("Wrong header and data blocks", "Trying to read a file that does not look correct").exec_() 
        
        else:
            # Now that everything has been seen to work,
            # save potential header and potential df as actual header and df
            self.reset_btns() 
            self.parent.last_loaded_file = os.path.split(filename)[0]
            self.parent.current_file = filename
            self.raw_df = potential_df
            self.raw_df_header = potential_header


            # Clear old data and set new names
            self.raw_fit_list.clear()
            self.fill_df_data_values()
            self.plot_type_combo.clear() 
            self.fit_combo.clear() 
            self.plot_type_combo.addItems(['Raw data'])
            
            self.cleanup_loaded_ppms()
            update_data_names(self.raw_df, self.parent.read_options)
            self.determine_data_type() 

            if self.data_type == "AC":

                self.num_meas_freqs = len(set(self.raw_df['AC Frequency (Hz)']))
                self.update_temp_subsets()
                self.update_meas_temps()

            if self.data_type == "DC": 

                self.update_subsets() 
            
            
            #Clearing all the plots
            self.raw_plot.clear_canvas()
            self.fit_plot.clear_canvas()
            self.fit_plot.cax.clear()
            self.threeD_plot.clear_canvas() 
            self.raw_plot.clear_canvas()
            self.HM_plot.clear_canvas()
            self.TXT_plot.clear_canvas() 

            # Clearing axes of "old" drawings and setting front widget to the raw data
            combo_idx = self.plot_type_combo.findText('Raw data')
            self.plot_type_combo.setCurrentIndex(combo_idx)
            
            # Updating analysis combos, which will automatically draw the new data
            
            self.update_xy_combos()
            
            #Updates "Table of Data" tab with the loaded data
            self.parent.widget_table.updatetable()

    def update_subsets(self): 
            old_time = 0.0
            i = 0
            split_idx = []
            subsets = []
            self.HM_subsets = []
            self.TXT_subset = None

            #Finds indicies to split data and asigns TXT or HM to each subset
            for i in range(self.raw_df.shape[0]):  
                current_time = self.raw_df["Time Stamp (sec)"].values[i]
                if current_time > old_time + 500:
                    split_idx.append(i)
                    if isclose(self.raw_df["Magnetic Field (Oe)"].values[i], self.raw_df["Magnetic Field (Oe)"].values[i+10], rel_tol = 0.05): 
                        subsets.append((i,"TXT"))
                    if isclose(self.raw_df["Temperature (K)"].values[i],self.raw_df["Temperature (K)"].values[i+10],rel_tol = 0.05): 
                        if self.raw_df["Temperature (K)"].values[i] < 250: #Assumption: You would never measure T vs. XT data above 250 K  
                            subsets.append((i,"HM"))
                old_time = current_time
            split_idx.append(self.raw_df.shape[0])


            #Makes subsets
            for i in range(len(split_idx)): 
                if (split_idx[i], "HM") in subsets: 
                    self.HM_subsets.append(self.raw_df.iloc[split_idx[i]:split_idx[i+1]])
                if (split_idx[i], "TXT") in subsets: 
                    self.TXT_subset = self.raw_df.iloc[split_idx[i]:split_idx[i+1]]
            #Ads each HM subset temperature to measured temperature and finds min and max temperature of all HM subsets
            if self.HM_subsets != []: 
                meas_temps = []
                for sub in self.HM_subsets:
                    meas_temps.append(sub['Temperature (K)'].mean())

                self.meas_temps = np.array(meas_temps)
                self.num_meas_temps = len(self.meas_temps)
                self.Tmin = self.meas_temps.min()
                self.Tmax = self.meas_temps.max()


    def make_diamag_correction_calculation(self):
        """ Calculates the diamagnetic correction if data is loaded, and 
        inserts these into the dataframe in their molar form: Mp_m, Mpp_m, Xp_m, Xpp_m.
        Afterwards, the temperature subsets, xy_combos and the datatable are updated. """

        if self.raw_df is None:
            # Don't do the calculation, if there is nothing to calculate on
            pass
        
        else:
            try:
                m_sample = float(self.m_sample)
                M_sample = float(self.M_sample)
                Xd_sample = float(self.Xd_sample)
                constant_terms = [float(x) for x in self.constant_terms.split(',')]
                var_am = [float(x) for x in self.var_am.split(',')]
                
                assert len(var_am)%2==0
                paired_terms = [(var_am[n], var_am[n+1]) for n in range(0,len(var_am),2)]
                
                if Xd_sample == 0:
                    Xd_sample = -6e-7*M_sample
                
            except (ValueError, AssertionError, AttributeError):
                MagMessage('Error', 'Something wrong in "Sample information"\n').exec_()
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
                self.update_xy_combos()
                self.parent.widget_table.updatetable()
                MagMessage('Diamagnetic correction',
                               'Diamagnetic correction successful!').exec_()


    
    def update_raw_fit_list(self):   
        """Updates the raw_fit_list in the ListWidget, where it can be chosen which fits and raw 
        data you want to visualize in the plotting window. Uses colormap. """

        self.raw_fit_list.clear()
        for i in range(self.num_meas_temps):
            T = self.meas_temps[i]
            newitem = QListWidgetItem()
            newitem.setText('T = {:<6.2f} K, Show raw data: {}, Show fit/line: {}'.format(round(T,2),True, True)) #Text for Fitted Parameters box
            plotting_dict = {'temp': self.meas_temps[i],
                             'raw': True,
                             'fit': True}
            newitem.setData(32, plotting_dict)
            t_float = (T-self.Tmin)/(self.Tmax-self.Tmin) #Makes t_float scaled by temperature from 0 to 1 for colormap
            newitem.setBackground(QColor(to_hex(self.parent.temperature_cmap(t_float))))
            self.raw_fit_list.addItem(newitem)



    def plot_from_itemlist(self):
        """ Plots the fitted plot depending on what is chosen in the fitted combobox 
        (either Freq vs. χ', Freq. vs. χ'' or Cole-Cole) """
        
        if self.data_type == "DC": 
            return
        if self.raw_fit_list.count()==0:
            return
        
        self.fit_plot.ax.clear()
        plot_type = self.fit_combo.currentText()
        
        if plot_type == 'Freq VS Xp':
            x_name = 'AC Frequency (Hz)'
            y_name = 'Xp_m (emu/(Oe*mol))'
            fcn_y = Xp_
            x_scale = 'log'
        elif plot_type == 'Freq VS Xpp':
            x_name = 'AC Frequency (Hz)'
            y_name = 'Xpp_m (emu/(Oe*mol))'
            fcn_y = Xpp_
            x_scale = 'log'
        elif plot_type == 'Cole-Cole':
            x_name = 'Xp_m (emu/(Oe*mol))'
            y_name = 'Xpp_m (emu/(Oe*mol))'
            fcn_y = Xpp_
            x_scale = 'linear'
            
        for row in range(self.num_meas_temps):
        
            T = self.meas_temps[row]
            rgb = self.parent.temperature_cmap((T-self.Tmin)/(self.Tmax-self.Tmin))

            
            if plot_type == 'Cole-Cole':
                x_data = Xp_(self.temp_subsets[row]['AC Frequency (Hz)'],
                             self.raw_data_fit['ChiS'].iloc[row],
                             self.raw_data_fit['ChiT'].iloc[row],
                             self.raw_data_fit['Tau'].iloc[row],
                             self.raw_data_fit['Alpha'].iloc[row])
            else:
                x_data = self.temp_subsets[row][x_name]
                
            item = self.raw_fit_list.item(row)
            itemdict = item.data(32)
            if itemdict['raw']:
                self.fit_plot.ax.plot(self.temp_subsets[row][x_name],
                                            self.temp_subsets[row][y_name],
                                            marker='o',
                                            mec=rgb,
                                            mfc='none',
                                            linestyle='None')
            if itemdict['fit']:
                self.fit_plot.ax.plot(x_data,
                                            fcn_y(self.temp_subsets[row]['AC Frequency (Hz)'],
                                                  self.raw_data_fit['ChiS'].iloc[row],
                                                  self.raw_data_fit['ChiT'].iloc[row],
                                                  self.raw_data_fit['Tau'].iloc[row],
                                                  self.raw_data_fit['Alpha'].iloc[row]),
                                            c=rgb)
            
        self.fit_plot.ax.set_xscale(x_scale)
        self.fit_plot.ax.set_xlabel(formatlabel(x_name))
        self.fit_plot.ax.set_ylabel(formatlabel(y_name))

        norm = mpl.colors.Normalize(vmin=self.Tmin, vmax=self.Tmax)
        self.fit_plot.fig.colorbar(
            mpl.cm.ScalarMappable(norm=norm,
                                  cmap=self.parent.temperature_cmap),
                                        orientation='horizontal',
            cax=self.fit_plot.cax)
        
        self.fit_plot.canvas.draw()
        
        combo_idx = self.plot_type_combo.findText('Fitted (AC)')
        self.plot_type_combo.setCurrentIndex(combo_idx)

    def fit_Xp_Xpp_standalone(self):
        """ Fits χ' and χ'' for all temperature subsets. This gives a tau, tau fit error, alpha, 
        χs, χt and residual to every temperature subset. 
        It also updates the raw_fit_list, plots the fitted data, plots the 3D data
        and updates the xy comboboxes with the new variables introduced in the dataframe """
        
        if 'Fitted (AC)' not in [self.plot_type_combo.itemText(i) for i in range(self.plot_type_combo.count())]: 
            self.plot_type_combo.addItems(['Fitted (AC)', '3D plot (AC)'])
        if "Cole-Cole" not in [self.fit_combo.itemText(i) for i in range(self.fit_combo.count())]: 
            self.fit_combo.addItems(['Cole-Cole', 'Freq VS Xp', 'Freq VS Xpp'])

        try:
            # Check to see if there has been loaded a data frame
            self.raw_df.columns
            # Check to see if the data to work on is in the data frame
            assert 'Xp_m (emu/(Oe*mol))' in self.raw_df.columns
            
        except AttributeError:
            MagMessage("Error", "There hasn't been loaded a data frame to work on").exec_() 
        except AssertionError:
            MagMessage('Error','Calculate diamagnetic correction first\nto make Xp_m and Xpp_m for the algorithm').exec_()
        
        else:
            self.parent.statusBar.showMessage('Running fit...')
            
            T = [x for x in self.meas_temps]
            Xs, Xt, tau, alpha, resid, tau_fit_err = [],[],[],[],[],[]
            
            v_all = [np.array(self.temp_subsets[t_idx]['AC Frequency (Hz)']) for t_idx in range(self.num_meas_temps)]
            Xp_all = [np.array(self.temp_subsets[t_idx]['Xp_m (emu/(Oe*mol))']) for t_idx in range(self.num_meas_temps)]
            Xpp_all = [np.array(self.temp_subsets[t_idx]['Xpp_m (emu/(Oe*mol))']) for t_idx in range(self.num_meas_temps)]
            
            inputs = tuple(zip(v_all, Xp_all, Xpp_all))
            
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
            self.update_raw_fit_list()
            self.plot_from_itemlist()
            self.plot3D()
            self.update_xy_combos() 
            set_idx = self.plot_type_combo.findText('Fitted')
            self.plot_type_combo.setCurrentIndex(set_idx)
            
            self.parent.statusBar.showMessage("Fit of X' and X'' complete")
            self.save_fit_X_btn.setEnabled(True)
            self.copy_fit_btn.setEnabled(True)
            self.headline.setText("Fitted/3D plot: Hide/show raw data and fitted lines")


    def save_fit_to_file(self):
        """ The fit of χ' and χ'' is saved to a file with columns: 
        Temp, Tau, dTay, Alpha, ChiS, ChiT, Residual and Tau_err. 
        This file can be loaded in the data analysis tab. """
        
        name = QFileDialog.getSaveFileName(self, 'Save File')
        filename = name[0]
        if self.raw_data_fit is None:
            MagMessage("There is nothing to save", "There is probably no fit yet...").exec_()
        elif filename=='':
            pass 
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

            MagMessage("File succesfully saved", "A file with your fit of tau was succesfully saved!").exec_()


    def plot3D(self):
        """Plots 3D plot with temperature vs. frequency vs. X'' (molar) """

        if self.raw_fit_list.count()==0:
            return
        
        self.threeD_plot.ax.clear()

        self.add_labels_3Dplot()
        self.plot_each_temp_3D() 
        self.add_colorbar_3D() 

        self.threeD_plot.canvas.draw()

    def add_labels_3Dplot(self):
        """ Adds labels to the 3D plot. Logaritmic scale is used for the frequency. """
        
        def log_tick_formatter(val, pos=None):
            return f"$10^{{{int(val)}}}$" 
        
        self.x_label_3D = 'Temperature (K)'
        self.y_label_3D = "AC Frequency (Hz)"
        self.z_label_3D = "Xpp_m (emu/(Oe*mol))"
        
        self.threeD_plot.ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
        self.threeD_plot.ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        self.threeD_plot.ax.set_xlabel(formatlabel(self.x_label_3D))
        self.threeD_plot.ax.set_ylabel(formatlabel(self.y_label_3D))
        self.threeD_plot.ax.set_zlabel(formatlabel(self.z_label_3D))             
        
    def plot_each_temp_3D(self): 
        """ Plots each temperature subset in colors according to the colorbar"""
        for row in range(self.num_meas_temps):      
            T = self.meas_temps[row]
            rgb = self.parent.temperature_cmap((T-self.Tmin)/(self.Tmax-self.Tmin))           

            item = self.raw_fit_list.item(row)
            itemdict = item.data(32)
            if itemdict['raw']:
                self.threeD_plot.ax.scatter3D(self.temp_subsets[row][self.x_label_3D],
                                        np.log10(self.temp_subsets[row][self.y_label_3D]),
                                        self.temp_subsets[row][self.z_label_3D], 
                                        color = rgb, s = 7
                                        )
            if itemdict['fit']:       
                self.threeD_plot.ax.plot(self.temp_subsets[row][self.x_label_3D], 
                                            np.log10(self.temp_subsets[row][self.y_label_3D]),
                                            Xpp_(self.temp_subsets[row]['AC Frequency (Hz)'],
                                                  self.raw_data_fit['ChiS'].iloc[row],
                                                  self.raw_data_fit['ChiT'].iloc[row],
                                                  self.raw_data_fit['Tau'].iloc[row],
                                                  self.raw_data_fit['Alpha'].iloc[row]),
                                            c=rgb)  

    def add_colorbar_3D(self): 
        """ Adds colorbar to the 3D plot"""

        norm = mpl.colors.Normalize(vmin=self.Tmin, vmax=self.Tmax)
        self.threeD_plot.fig.colorbar(
            mpl.cm.ScalarMappable(norm=norm,
                                  cmap=self.parent.temperature_cmap),
                                        orientation='horizontal',
            cax=self.threeD_plot.cax)

    def switch_view(self):
        """ Switches the view of the stackingwidget according to what is chosen 
        in the axis content combobox """

        idx = self.plot_type_combo.currentIndex()
        if self.data_type != None: 
            if idx != 0 and self.data_type== "DC": 
                idx += 2
        self.sw.setCurrentIndex(idx)


    def plot_raw_data(self):
        """ Plots the raw data based on what is chosen in the xy comboxes"""

        self.raw_plot.ax.clear()

        x_label, y_label = self.x_combo.currentText(), self.y_combo.currentText()
        if x_label == '' or y_label == '': 
            x_label = self.raw_df.columns[0]
            y_label = self.raw_df.columns[0]

        self.raw_plot.ax.plot(self.raw_df[x_label],
                                    self.raw_df[y_label],
                                    marker='o',
                                    mfc='none',
                                    mec='k',
                                    linestyle='-',
                                    c='k',
                                    linewidth=1,
                                    )

        self.raw_plot.ax.set_xlabel(formatlabel(x_label))
        self.raw_plot.ax.set_ylabel(formatlabel(y_label))
        self.raw_plot.canvas.draw()
        
        combo_idx = self.plot_type_combo.findText('Raw data')
        self.plot_type_combo.setCurrentIndex(combo_idx)


    def update_itemdict(self, item, itemdict):
        """ Updates a given item in the ListWidget where on can chose which 
        raw data and fitted data to show. """

        item.setData(32, itemdict)
        item.setText('T = {:<6.2f} K, Show raw data: {}, Show fit: {}'.format(round(itemdict['temp'],2),
                                         itemdict['raw'],
                                         itemdict['fit']))

    def update_fitted_and_3D_plot(self): #Previously names update_raw_plot
        """ Updates the fitted plot with the data that has been chosen to be visualized in 
        the ListWidget """

        w = FitResultPlotStatus(list_input=self.raw_fit_list)
        finished_value = w.exec_()
        if not finished_value:
            pass
        else:
            final_states = w.checked_items
            for i, boxes in enumerate(final_states):
                item = self.raw_fit_list.item(i)
                item_data = item.data(32)
                item_data['raw'] = boxes[0].isChecked()
                item_data['fit'] = boxes[1].isChecked()
                self.update_itemdict(item, item_data)
        if self.data_type == "AC": 
            self.plot_from_itemlist()
            self.plot3D() 
        elif self.data_type == "DC": 
            try:
                self.plot_HM() 
            except: 
                pass
            try: 
                self.plot_XT
            except: 
                pass
        self.fit_plot.canvas.draw()
        self.threeD_plot.canvas.draw()
        self.HM_plot.canvas.draw() 


    def copy_fit_to_analysis(self):
        """ Copies the fit of tau obtained after fitting χ' and χ'' 
        to the data analysis tab. """
    
        try:
            D = np.array(list(zip(self.meas_temps,
                                  self.raw_data_fit['Tau'],
                                  self.raw_data_fit['dTau'])))
            self.parent.data_ana.set_new_t_tau(D)
            self.parent.data_ana.read_indices_for_used_temps()
            self.parent.data_ana.plot_t_tau_on_axes()
            self.parent.data_ana.add_T_axis() 
            self.parent.data_ana.plot_wdgt.reset_axes()
            self.parent.data_ana.update_T_axis() 
        except TypeError:
            MagMessage("Fitted data does not exist", "Fitted data does not yet exist in the Data Treatment tab").exec_() 
        else: 
            MagMessage("Succes", "Your fit of tau was succesfully copied to the Data Analysis tab").exec_()

    def make_diamag_correction_calculation_dc(self):
        """ Calculates the diamagnetic correction if data is loaded, and 
        inserts these into the dataframe in their molar form: Mp_m, Mpp_m, Xp_m, Xpp_m.
        Afterwards, the temperature subsets, xy_combos and the datatable are updated. """

        if self.raw_df is None:
            # Don't do the calculation, if there is nothing to calculate on
            pass
        #elif "Moment_m (emu/mol)" not in self.raw_df: 
        #    print("ERROR, no molar moment")

        else:
            if "H VS M" not in [self.plot_type_combo.itemText(i) for i in range(self.plot_type_combo.count())]: 
                self.plot_type_combo.addItems(["H VS M", "T VS XT"])
            try:
                m_sample = float(self.m_sample)
                M_sample = float(self.M_sample)
                Xd_sample = float(self.Xd_sample)
                constant_terms = [float(x) for x in self.constant_terms.split(',')]
                var_am = [float(x) for x in self.var_am.split(',')]

                assert len(var_am)%2==0
                paired_terms = [(var_am[n], var_am[n+1]) for n in range(0,len(var_am),2)]

                if Xd_sample == 0:
                    Xd_sample = -6e-7*M_sample
                
            except (ValueError, AssertionError, AttributeError):
                MagMessage('Error', 'Something wrong in "Sample information"\n').exec_()
            else: 
                H0 = self.raw_df['Magnetic Field (Oe)']
                M = self.raw_df["Moment (emu)"]
                
                # Get molar, corrected values from function in process_ac
                M_molar, X_molar = diamag_correction_dc(
                    H0, M, m_sample, M_sample, Xd_sample,
                    constant_terms = constant_terms, paired_terms = paired_terms)
                          # PUT THE CODE HERE TO INSERT CORRECTED VALUES INTO DATA FRAME
                if "Moment_m (emu/mol)" in self.raw_df.columns:
                    self.raw_df.replace(to_replace="Moment_m (emu/mol)", value=M_molar)
                if "X_m (emu/(Oe*mol))" in self.raw_df.columns: 
                    self.raw_df.replace(to_replace="X_m (emu/(Oe*mol))", value=X_molar)
                else:
                    try: 
                        M_idx = self.raw_df.columns.get_loc('Moment (emu)')
                        self.raw_df.insert(M_idx+1, column = "X_m (emu/(Oe*mol))", value=M_molar)   
                        self.raw_df.insert(M_idx+2, column = "XT_m (emu*K/Oe)", value = self.raw_df["X_m (emu/(Oe*mol))"]*self.raw_df["Temperature (K)"])
                    except KeyError: 
                            pass
                    if self.HM_subsets != []: 
                        M_idx = self.raw_df.columns.get_loc('Moment (emu)')
                        self.raw_df.insert(M_idx+1, column="Moment_m (emu/mol)", value=M_molar)


                self.update_subsets()
                self.update_xy_combos()
                self.update_raw_fit_list()

                try: 
                    self.plot_XT()
                except: 
                    pass 

                if self.HM_subsets != []: 
                    self.plot_HM() 

                if self.HM_subsets != []: 
                    combo_idx = self.plot_type_combo.findText('H VS M')
                else: 
                    combo_idx = self.plot_type_combo.findText('T VS XT')
                self.plot_type_combo.setCurrentIndex(combo_idx)
                self.headline.setText("H VS M: Hide/show raw data and fitted lines") 
                MagMessage('Diamagnetic correction',
                            'Diamagnetic correction successful!').exec_()



    def plot_XT(self): #Maybe it should be just XT instead of TXT

        x_name = "Temperature (K)"
        y_name = "XT_m (emu*K/Oe)"
        self.TXT_plot.ax.plot(self.TXT_subset[x_name],
                            self.TXT_subset[y_name],
                            marker='o',
                            mfc='none',
                            markersize = 5) 

        self.TXT_plot.ax.set_xlabel(formatlabel(x_name))
        self.TXT_plot.ax.set_ylabel(formatlabel(y_name))
        self.TXT_plot.canvas.draw() 


    def plot_HM(self): 

        self.HM_plot.ax.clear()

        try: 
            for row in range(self.num_meas_temps):
                
                T = self.meas_temps[row]
                if len(self.HM_subsets) == 1: 
                    rgb = 'k'
                else:
                    rgb = self.parent.temperature_cmap((T-self.Tmin)/(self.Tmax-self.Tmin))
    
                x_name = "Magnetic Field (Oe)"
                y_name = "Moment_m (emu/mol)"
                
                item = self.raw_fit_list.item(row)
                itemdict = item.data(32)
                if itemdict['raw']:
                    self.HM_plot.ax.plot(self.HM_subsets[row][x_name],
                                        self.HM_subsets[row][y_name],
                                        marker='o',
                                        color=rgb,
                                        mfc='none',
                                        markersize = 3, 
                                        linestyle = "none") 
                if itemdict['fit']:       
                    self.HM_plot.ax.plot(self.HM_subsets[row][x_name],
                                        self.HM_subsets[row][y_name],
                                        linestyle = '-',
                                        color=rgb,
                                        mfc='none',
                                        markersize = 3) 
    
                
            self.HM_plot.ax.set_xlabel(formatlabel(x_name))
            self.HM_plot.ax.set_ylabel(formatlabel(y_name))

            if len(self.HM_subsets) > 1: 
                norm = mpl.colors.Normalize(vmin=self.Tmin, vmax=self.Tmax)
                self.HM_plot.fig.colorbar(
                    mpl.cm.ScalarMappable(norm=norm,
                                          cmap=self.parent.temperature_cmap),
                                                orientation='horizontal',
                    cax=self.HM_plot.cax)

            self.HM_plot.canvas.draw()      
        except KeyError: #If errors in sample information  
            pass





    def update_sample_info_and_make_dia_correction(self): 
        """ Updates the sample information """
        w = SampleInformation(self.parent, self)
        w.exec_()
        
