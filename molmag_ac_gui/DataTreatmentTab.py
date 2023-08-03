#std packages
import os
from multiprocessing import Pool

#third-party packages 
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.colors import to_hex
import matplotlib.ticker as mticker
from PyQt5.QtWidgets import (QFileDialog, QFormLayout, QListWidgetItem, QWidget, 
                             QVBoxLayout, QHBoxLayout, QComboBox, QStackedWidget, 
                             QListWidget, QSplitter, QMessageBox, QApplication)
from PyQt5.QtGui import QColor
from torch import multiply
from PyQt5.QtCore import Qt

from molmag_ac_gui.layout import make_btn, make_checkbox, make_headline, make_line, headline_font

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
        # Initilization: Sets attributes to empty, zero or None 
        self.initialize_attributes()

        # Make the left column: data loading, fitting, visualization controls and show/hide option
        self.layout = QVBoxLayout()
        self.options_wdgt = QWidget()
        
        # Constructing data import and treatment buttons
        make_headline(self, "Data import and treatment", self.layout)        
        self.add_data_treatment_buttons() 

        # Constructing axis content combobox
        make_headline(self, "Content of plotting window", self.layout)
        self.add_plot_type_combobox()

        # Constructing the x and y comboboxes
        make_headline(self, "Raw data: Choose x- and y-contents", self.layout)
        self.add_xy_comboboxes() 

        # Constructing a combobox for plotting fitted AC data
        make_headline(self, "Fitted: Choose plot type (AC)", self.layout)
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
        self.pre_header = None
        self.num_meas_freqs = 0
        self.num_meas_temps = 0
        self.num_meas_Hs = 0
        self.temp_subsets = []
        self.meas_temps = []
        self.meas_Hs = []
        self.Tmin, self.Tmax = 0,0
        self.fit_result = None
        self.data_type = None
        self.T_legend = None
        self.DC_data_types = []
        self.data_to_export = None
    


    def add_plot_type_combobox(self): 
        """ Adds a plot_type combobox such that it is possible to switch between viewing 
        the raw data, fitted data and the 3D plot in the plotting window on the right. """

        self.plot_type_combo = QComboBox()
        self.plot_type_combo.currentIndexChanged.connect(self.switch_view)
        self.layout.addWidget(self.plot_type_combo)
        
 

    def add_fit_combobox(self): 
        """ Adds a combobox to choose which fitted data plot you want to view. """
        
        self.fit_combo = QComboBox()
        #self.fit_combo.addItems(['Cole-Cole', 'Freq VS Xp', 'Freq VS Xpp'])
        self.fit_combo.currentIndexChanged.connect(self.plot_from_itemlist)
        self.layout.addWidget(self.fit_combo)
        

    def add_data_treatment_buttons(self): 
        """ Adds the first 6 buttons in the option widget on the left side of the gui window. 
        and connects these to the corresponding functions. """

        make_btn(self, "(1) Load datafile", self.load_data, self.layout)
        make_btn(self, "(2) Make diamagnetic correction", self.update_sample_info_and_make_dia_correction, self.layout)
        self.fit_X_or_get_phi_btn = make_btn(self, "", self.fit_Xp_Xpp_standalone, self.layout)
        self.fit_X_or_get_phi_btn.setEnabled(False)
        #self.get_phi_files = make_btn(self, "(3) Git files for fitting in PHI (DC)", self.get_phi_files, self.layout)
        
        #Copy fit to data analysis and save fit to file layout
        self.copy_save_layout = QHBoxLayout()
        self.copy_fit_btn = make_btn(self, "Copy fit to Data analysis (AC)", self.copy_fit_to_analysis, self.copy_save_layout)
        self.copy_fit_btn.setEnabled(False)
        self.save_fit_X_btn = make_btn(self, "Save fit to file (AC)", self.save_fit_to_file, self.copy_save_layout)
        self.save_fit_X_btn.setEnabled(False)

        self.layout.addLayout(self.copy_save_layout)

    def get_phi_files(self): 

        #try: 
        self.find_phi_path()
        if os.path.isfile(r'{}.{}'.format(self.path_to_export,"exp")): #If the file already exist
            qm = QMessageBox() 
            ans = qm.question(self,'', "File already exists. Do you want to overwrite existing file?", qm.Yes | qm.No)
            if ans == qm.No: 
                MagMessage("File not saved", "Your file has not been saved").exec_()
            if ans == qm.Yes: 
                self.save_phi_files() 
        else: #If file does not already exist
            self.save_phi_files() 


        #except AttributeError: #If raw_df is empty, it will be NoneType. This type has no attribute "to_csv" or "to_excel"
        #    MagMessage("Error", "Cannot export file. No data is loaded in the data treatment tab.").exec_() 

    def get_df_to_save(self, subsets, x_name, y_name, HM_or_XT):
        if len(subsets) > 1: 
            df = pd.DataFrame(subsets[0][x_name])
            for subset in subsets[1:]: 
                df = df.append(pd.DataFrame(subset[x_name]), ignore_index=True)
            len_prev_dfs = 0 
            i = 1
                
            moment_df = pd.DataFrame() 

            for subset in subsets: 

                moment_df = moment_df.iloc[0:0]

                if len_prev_dfs != 0: 
                    moment_df = moment_df.append(pd.DataFrame("#", index = np.arange(len_prev_dfs), columns = [y_name,]), ignore_index= True) 
                    moment_df = moment_df.append(pd.DataFrame(subset[y_name]), ignore_index=True)
                      
                else: 
                    moment_df = moment_df.append(pd.DataFrame(subset[y_name]),ignore_index=True)

                len_prev_dfs += subset[y_name].shape[0]
                idx = df.shape[0] - moment_df.shape[0]
                df_to_append = pd.DataFrame("#", index = np.arange(idx), columns = [y_name,])
                moment_df = moment_df.append(df_to_append, ignore_index=True)
                if HM_or_XT == "HM": 
                    column_name = "M" + str(i)
                if HM_or_XT == "XT": 
                    column_name = "XT" + str(i)
                df.insert(i, column = column_name, value = moment_df) 
                i+=1

        else:
            df = subsets[0][[x_name, y_name]]
        return df
        

    def save_phi_files(self): 
        if "XT" in self.DC_data_types:
            df = self.get_df_to_save(self.XT_subsets, "Temperature (K)", "XT_m (cm3*K/mol)", "XT")
            try: 
                df.to_csv(r'{}.{}'.format(self.path_to_export + "_sus", "exp"), index = False, header = False, sep = " ")
                MagMessage('The data is successfully exported', "The data is saved as an .exp file at: {}.exp".format(self.path_to_export + "_sus")).exec_()  
            except: 
                MagMessage("Error", "Some error occured, you cannot save the H vs. M data").exec_() 

        if "HM" in self.DC_data_types: 
            df = self.get_df_to_save(self.HM_subsets, "Magnetic Field (Oe)", "Moment_m (Na*muB)", "HM")

            try: 
                df.to_csv(r'{}.{}'.format(self.path_to_export + "_mag", "exp"), index = False, header = False, sep = " ")
                MagMessage('The data is successfully exported', "The data is saved as an .exp file at: {}.exp".format(self.path_to_export + "_mag")).exec_()  
            except:
                MagMessage("Error", "Some error occured, you cannot save the H vs. M data").exec_() 




    def find_phi_path(self): 
        #Finds path for datatable csv or excel file to be stored

        splitpath = self.parent.current_file.split('/')
        newpath = ""
        for element in splitpath[:-1]: 
            newpath += element + '/'
        newpath += splitpath[-1].split(".")[0]  
        self.path_to_export = newpath

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
        self.xy_options_cb = make_checkbox(self,self.update_xy_combos,self.xy_options_lo)
        self.layout.addLayout(self.xy_options_lo)

    def add_stackwidget(self): 
        """ Adds the stackwidget which is the plotting window on the right side of the gui window"""

        self.sw = QStackedWidget()
        
        self.raw_plot = PlottingWindow()
        self.fit_plot = PlottingWindow(make_ax="cax")
        self.threeD_plot = PlottingWindow(make_ax = "z") 
        self.HM_plot = PlottingWindow() #make_ax="cax")
        self.XT_plot = PlottingWindow()

        self.sw.addWidget(self.raw_plot)
        self.sw.addWidget(self.fit_plot)
        self.sw.addWidget(self.threeD_plot)
        self.sw.addWidget(self.HM_plot)
        self.sw.addWidget(self.XT_plot)

        self.addWidget(self.sw)

    def add_parameter_controls(self): 
        """ Adds controls for the data plotting where one can pick which datapoints
        and fits to be visualized. """
        
        make_headline(self, "Hide/show raw data and lines", self.layout)
        #self.headline = QLabel( "Hide/show raw data and lines")
        #self.headline.setFont(headline_font)

        self.layout.addWidget(self.headline)
        
        make_line(self, "Double click to edit list", self.layout)
        
        #List of fitted raw data
        self.raw_fit_list = QListWidget()
        self.layout.addWidget(self.raw_fit_list)
        self.raw_fit_list.doubleClicked.connect(self.update_fitted_and_3D_plot)

        make_headline(self, "H vs. M plotting options", self.layout)
        self.T_legend_lo = QHBoxLayout() 
        make_line(self, "Add legend with temperatures to H vs. M plot", self.T_legend_lo)
        self.T_legend_cb = make_checkbox(self, self.add_T_legend, self.T_legend_lo)

        self.layout.addLayout(self.T_legend_lo)

        
    """ Other functions"""
    
    def add_T_legend(self): 
        """ Adds a temperature legend to the H vs. M plot"""
        
        if self.T_legend_cb.isChecked(): 
            self.T_legend = True
        else: 
            self.T_legend = False 
        if self.data_type == "DC": 
            try: 
                self.plot_HM() 
            except: 
                pass 


    def fill_df_data_values(self):
        """ If Xp is in the loaded data, and Mp is not: Mp will be calculated and added to the
        dataframe. If Mp is in the loaded data and Xp is not, Xp will be calculated and added 
        to the dataframe. """
        Htot = self.raw_df['Magnetic Field (Oe)'] + self.raw_df['AC Amplitude (Oe)']
        
        if ('Xp (emu/Oe)' in self.raw_df.columns and not ('Mp (emu)' in self.raw_df.columns)):
            # Susceptibility exists in the data frame, but magnetisation does not
            Mp = self.raw_df['Xp (emu/Oe)']*Htot
            Mpp = self.raw_df['Xpp (emu/Oe)']*Htot
            Htot = self.raw_df['Magnetic Field (Oe)'] + self.raw_df['AC Amplitude (Oe)']
            Mp = self.raw_df['Xp (emu/Oe)']*Htot
            Mpp = self.raw_df['Xpp (emu/Oe)']*Htot
            Xp_idx = self.raw_df.columns.get_loc('Xp (emu/Oe)')
            self.raw_df.insert(Xp_idx, column='Mp (emu)', value=Mp)
            self.raw_df.insert(Xp_idx+1, column='Mpp (emu)', value=Mpp)
            
        elif (not 'Xp (emu/Oe)' in self.raw_df.columns and ('Mp (emu)' in self.raw_df.columns)):
            # Magnetisation exists in the data frame, but susceptibility does not
            Xp = self.raw_df['Mp (emu)']/Htot
            Xpp = self.raw_df['Mpp (emu)']/Htot
            Mp_idx = self.raw_df.columns.get_loc('Mp (emu)')
            self.raw_df.insert(Mp_idx+2, column='Xp (emu/Oe)', value=Xp)
            self.raw_df.insert(Mp_idx+3, column='Xpp (emu/Oe)', value=Xpp)

    def cleanup_loaded_ppms(self):
        """Cleans up the data that is loaded in by removing all empty columns,
        removes Comment column and renames the columns if necessary"""

        # Drop columns where all values are NaN
        self.raw_df.dropna(axis=1, how='all', inplace=True)

        # Drop "Comment" column
        if 'Comment' in self.raw_df.columns:
            self.raw_df.drop(['Comment'], axis='columns', inplace=True)
        
        # Drop all rows where there is still a NaN value
        self.raw_df.dropna(axis=0, inplace=True)
        
        # Make sure that the rows are named correctly after dropping some columns
        old_indices = self.raw_df.index.values
        new_indices = list(range(len(old_indices)))
        self.raw_df.rename(index=dict(zip(old_indices, new_indices)),
                           inplace=True)

    def update_temp_subsets(self):
        """ A function for AC data. Splits datasets into one dataset for each temperature. 
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
        """ Updates the measured temperatures based on the temperatures in self.temp_subsets.
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
        if self.raw_df is None: 
            return #Dont do anything
        self.x_combo.clear()
        self.y_combo.clear()
        if self.data_type == "AC":
            chosenlabels = ["Temperature (K)", "AC Frequency (Hz)", "AC Amplitude (Oe)", "Magnetic Field (Oe)", "Mp (emu)",\
                           "Mpp (emu)", "Xp (emu/Oe)", "Xpp (emu/Oe)"]
            molarlabels = ["Mp_m (emu/mol)", "Mpp_m (emu/mol)", "Xp_m (emu/(Oe*mol))", "Xpp_m (emu/(Oe*mol))"]
        elif self.data_type == "DC": 
            chosenlabels = ["Time Stamp (sec)", "Temperature (K)", "Moment (emu)","Magnetic Field (Oe)"]
            if self.DC_data_types == ["HM"]:
                molarlabels = ["Moment_m (emu/mol)", "Moment_m (Na*muB)"]
            else: 
                molarlabels = ["Moment_m (emu/mol)", "Moment_m (Na*muB)", "X_m (emu/(Oe*mol))", "XT_m (emu*K/(Oe*mol))", "XT_m (cm3*K/mol)"]


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
        """ Determines whether the input data is AC or DC data.
        Tryies to find a column in raw_df named "AC Amplitude (Oe)", 
        if this is present, the data type is assumed to be AC data. 
        Otherwise, the program assumes DC data. """

        try: #Checks if it is AC data
            self.raw_df['AC Amplitude (Oe)']
            self.data_type = "AC"

        except: #Otherwise expects DC data 
            self.data_type = "DC"


    def reset_btns(self):
        """ Blocks the buttons Fit χ' and χ'', Copy fit to Data analysis tab and 
        save fit into file. These buttons cannot be clicked unless enabled """

        self.fit_X_or_get_phi_btn.setEnabled(False)
        #self.fit_X_or_get_phi_btn.clicked.disconnect() 
        self.copy_fit_btn.setEnabled(False)
        self.save_fit_X_btn.setEnabled(False)

    def load_data(self):
        """ This function uses subfunctions to perform the following tasks:  
        Loads the ppms data, clears the dataframe, fills in the dataframe with the loaded raw data, 
        adds Mp or Xp to the dataframe, cleans up the dataframe, updates temperatures subsets, 
        clears the raw data and fitted plots for any old data, plots raw data and updates the data table. """
        self.initialize_attributes()
        succes_in_loading = self.try_load_raw_df()
        if succes_in_loading: ## If the function succeeds at loading the raw_df, this statement will be True

            self.reset_btns() 
            
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
                
                if len(self.HM_subsets) > 0: 
                    self.DC_data_types.append("HM")
                if len(self.XT_subsets) > 0:
                    self.DC_data_types.append("XT")
                
                print("Data types = ", self.DC_data_types)
            #Clearing all the plots
            self.clear_plots() 

            # Clearing axes of "old" drawings and setting front widget to the raw data
            combo_idx = self.plot_type_combo.findText('Raw data')
            self.plot_type_combo.setCurrentIndex(combo_idx)
            
            # Updating analysis combos, which will automatically draw the new data
            self.update_xy_combos()
            
            #Updates "Table of Data" tab with the loaded data
            self.parent.widget_table.update_table()

        QApplication.restoreOverrideCursor()

    def try_load_raw_df(self): 
        """ Tries to load the raw dataframe """
                
        open_file_dialog = QFileDialog()        
        filename_info = open_file_dialog.getOpenFileName(self, 'Open file', self.parent.last_loaded_file)
        filename = filename_info[0]
        #filename = "C:/Users/au592011/OneDrive - Aarhus Universitet/Skrivebord/TestData_MAG/ac-data/ac-data/dy-dbm/20180209DyII_1000.dat"
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            # FileNotFoundError and UnicodeDecodeError will be raised here
            pre_header, potential_header, potential_df = read_ppms_file(filename)
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
            self.parent.last_loaded_file = os.path.split(filename)[0]
            self.parent.current_file = filename
            self.raw_df = potential_df
            self.raw_df_header = potential_header
            self.pre_header = pre_header
            return True
        
        return False


    def clear_plots(self): 
        """ Clears all plots: raw fit, fitted, fitted cax, 3D, HM plot and XT plot """

        self.raw_plot.clear_canvas()
        self.fit_plot.clear_canvas()
        self.fit_plot.cax.clear()
        self.threeD_plot.clear_canvas() 
        self.HM_plot.clear_canvas()
        self.XT_plot.clear_canvas() 
    
    def update_subsets(self): 
        """ Updates the subsets of data for DC data. This function enables you to import a
        datafile that either contains only M VS H data, XT data or both.
        Splits the data based on when the time stamp jumps.
        Assigns HM or XT to each dataset depending on whether the magnetic field is almost constant 
        (within 0.05 Oe) between measurement 1 and 5 or whether the temperature is almost constant
        (within 0.05 K) between measurement 1 and 5"""

        old_time = 0.0
        i = 0
        split_idx = []
        subsets = []
        self.HM_subsets = []
        self.XT_subsets = []
        # Finds indicies to split data and asigns XT or HM to each subset
        for i in range(self.raw_df.shape[0]-5):  
            current_time = self.raw_df["Time Stamp (sec)"].values[i]
            if current_time > old_time + 100: # If time stamp jumps by more than 100 seconds, the data set is split into a new subset.
                split_idx.append(i)
                if isclose(self.raw_df["Magnetic Field (Oe)"].values[i], self.raw_df["Magnetic Field (Oe)"].values[i+5], rel_tol = 0.05): 
                    subsets.append((i,"XT"))
                if isclose(self.raw_df["Temperature (K)"].values[i],self.raw_df["Temperature (K)"].values[i+5],rel_tol = 0.05): #Used to be rel_tol = 0.01 
                    if self.raw_df["Temperature (K)"].values[i] < 250: #Assumption: You would never measure T vs. XT data above 250 K  
                        subsets.append((i,"HM"))
                        
            
            old_time = current_time 
        split_idx.append(self.raw_df.shape[0])

        # Makes subsets of the dataframes 
        for i in range(len(split_idx)): 
            if (split_idx[i], "HM") in subsets: 
                self.HM_subsets.append(self.raw_df.iloc[split_idx[i]:split_idx[i+1]])
            if (split_idx[i], "XT") in subsets: 
                self.XT_subsets.append(self.raw_df.iloc[split_idx[i]:split_idx[i+1]])
        



        self.set_HM_temps() 
        self.set_XT_Hs()




    def set_HM_temps(self): 
        """ Adds each HM subset temperature to the list of measured temperatures
         and finds min and max temperature of all HM subsets"""

        if self.HM_subsets != []: 
            meas_temps = []
            for sub in self.HM_subsets:
                meas_temps.append(sub['Temperature (K)'].mean())
            self.meas_temps = np.array(meas_temps)
            self.num_meas_temps = len(self.meas_temps)
            self.Tmin = self.meas_temps.min()
            self.Tmax = self.meas_temps.max()
    
    def set_XT_Hs(self): 
        """ Adds each HM subset temperature to the list of measured temperatures
         and finds min and max temperature of all HM subsets"""

        if self.XT_subsets != []: 
            meas_Hs = []
            for sub in self.XT_subsets:
                meas_Hs.append(sub['Magnetic Field (Oe)'].mean())
            self.meas_Hs = np.array(meas_Hs)
            self.num_meas_Hs = len(self.meas_Hs)
            self.Hmin = self.meas_Hs.min()
            self.Hmax = self.meas_Hs.max()
        
        self.num_meas_Hs = len(self.meas_Hs)


    def make_diamag_correction_calculation(self):
        """ Calculates the diamagnetic correction if data is loaded, and 
        inserts these into the dataframe in their molar form: Mp_m, Mpp_m, Xp_m, Xpp_m.
        Afterwards, the temperature subsets, xy_combos and the datatable are updated. """

        if self.raw_df is None:
            # Don't do the calculation, if there is nothing to calculate on
            pass
        
        else:
            QApplication.setOverrideCursor(Qt.WaitCursor)

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
            
                self.insert_molar_values(Mp_molar, Mpp_molar, Xp_molar, Xpp_molar) 

                self.update_temp_subsets()
                self.update_xy_combos()

                self.fit_X_or_get_phi_btn.setEnabled(True)
                self.fit_X_or_get_phi_btn.clicked.disconnect()
                self.fit_X_or_get_phi_btn.clicked.connect(self.fit_Xp_Xpp_standalone)
                self.fit_X_or_get_phi_btn.setText("(3) Fit χ' and χ'' (AC)")

            QApplication.restoreOverrideCursor()


    def insert_molar_values(self, Mp_molar, Mpp_molar, Xp_molar, Xpp_molar): 
        """ Updates the molar values of M', M'', χ' and χ'' or adds them to the 
        dataframe (raw_df) if not already present """

        if "Mp_m (emu/mol)" in self.raw_df.columns:
            self.raw_df["Mp_m (emu/mol)"] = Mp_molar
            self.raw_df["Mpp_m (emu/mol)"] = Mpp_molar
            self.raw_df["Xp_m (emu/mol)"] = Xp_molar
            self.raw_df["Xpp_m (emu/mol)"] = Xpp_molar
            self.parent.widget_table.update_table()

        else:
            Mp_idx = self.raw_df.columns.get_loc('Mp (emu)')
            self.raw_df.insert(Mp_idx+1, column="Mp_m (emu/mol)", value=Mp_molar)
            
            Mpp_idx = self.raw_df.columns.get_loc('Mpp (emu)')
            self.raw_df.insert(Mpp_idx+1, column="Mpp_m (emu/mol)", value=Mpp_molar)
            
            Xp_idx = self.raw_df.columns.get_loc('Xp (emu/Oe)')
            self.raw_df.insert(Xp_idx+1, column="Xp_m (emu/(Oe*mol))", value=Xp_molar)
            
            Xpp_idx = self.raw_df.columns.get_loc('Xpp (emu/Oe)')
            self.raw_df.insert(Xpp_idx+1, column="Xpp_m (emu/(Oe*mol))", value=Xpp_molar)
            
            self.parent.widget_table.update_table()


    def update_raw_fit_list(self):   
        """ Updates the raw_fit_list in the ListWidget, where it can be chosen which fits and raw 
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
            if self.Tmin == self.Tmax: 
                t_float = 0.5 
            else: 
                t_float = (T-self.Tmin)/(self.Tmax-self.Tmin) #Makes t_float scaled by temperature from 0 to 1 for colormap

            newitem.setBackground(QColor(to_hex(self.parent.temperature_cmap(t_float))))
            self.raw_fit_list.addItem(newitem)



    def plot_from_itemlist(self):
        """ Plots the fitted plot depending on what is chosen in the fitted combobox 
        (either Freq vs. χ', Freq. vs. χ'' or Cole-Cole) """
        
        if self.data_type == "DC": 
            return
        fit_resolution = 10000

        if self.raw_fit_list.count()==0:
            return
        
        self.fit_plot.ax.clear()

        plot_type, x_name, y_name, fcn_y, x_scale = self.get_plot_settings() 

        for row in range(self.num_meas_temps):
        
            T = self.meas_temps[row]
            rgb = self.parent.temperature_cmap((T-self.Tmin)/(self.Tmax-self.Tmin))

            #Sets x-data to χ' if plot-type is Cole Cole,
            if plot_type == 'Cole-Cole':
                x_data = Xp_(self.temp_subsets[row]['AC Frequency (Hz)'],
                             self.fit_result['ChiS'].iloc[row],
                             self.fit_result['ChiT'].iloc[row],
                             self.fit_result['Tau'].iloc[row],
                             self.fit_result['Alpha'].iloc[row])
            else:
                x_data = self.temp_subsets[row][x_name]
            
            #Plots the data depending on what is chosen in the raw_fit_list
            item = self.raw_fit_list.item(row)
            itemdict = item.data(32)
            if itemdict['raw']:
                self.fit_plot.ax.plot(x_data,
                                            self.temp_subsets[row][y_name],
                                            marker='o',
                                            mec=rgb,
                                            mfc='none',
                                            linestyle='None')
            if itemdict['fit']:
                if plot_type == 'Cole-Cole':
                    freq = self.temp_subsets[row]['AC Frequency (Hz)']
                    freq_for_plotting = np.linspace(min(freq), max(freq), fit_resolution)
                    x_data_for_plotting = Xp_(freq_for_plotting,
                                              self.fit_result['ChiS'].iloc[row],
                                              self.fit_result['ChiT'].iloc[row],
                                              self.fit_result['Tau'].iloc[row],
                                              self.fit_result['Alpha'].iloc[row])
                    y_data_for_plotting = fcn_y(freq_for_plotting,
                                                      self.fit_result['ChiS'].iloc[row],
                                                      self.fit_result['ChiT'].iloc[row],
                                                      self.fit_result['Tau'].iloc[row],
                                                      self.fit_result['Alpha'].iloc[row])                                
                    self.fit_plot.ax.plot(x_data_for_plotting,
                                                y_data_for_plotting,
                                                c=rgb)
                else: #Corresponds to elif plot_type == 'Freq VS Xpp' or plot_type == 'Freq VS Xp':
                    freq_for_plotting = np.logspace(np.log10(min(x_data)), np.log10(max(x_data)), fit_resolution) #Logspace to equally distribute the points on a log-scale
                    self.fit_plot.ax.plot(freq_for_plotting,
                                                fcn_y(freq_for_plotting,
                                                      self.fit_result['ChiS'].iloc[row],
                                                      self.fit_result['ChiT'].iloc[row],
                                                      self.fit_result['Tau'].iloc[row],
                                                      self.fit_result['Alpha'].iloc[row]),
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

        combo_txt = self.plot_type_combo.currentText()
        if combo_txt == "Raw data" or combo_txt == "": 
            self.plot_type_combo.setCurrentText("Fitted (AC)")



    def get_plot_settings(self): 
        """Determines the x-data, y-data, function and x-scale 
        used for plotting in the Fitted window (AC) """

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
         
        return plot_type, x_name, y_name, fcn_y, x_scale
            

    def fit_Xp_Xpp_standalone(self):
        """ Fits χ' and χ'' for all temperature subsets. This gives a tau, tau fit error, alpha, 
        χs, χt and residual to every temperature subset. 
        It also updates the raw_fit_list, plots the fitted data, plots the 3D data
        and updates the xy comboboxes with the new variables introduced in the dataframe """
        
        QApplication.setOverrideCursor(Qt.WaitCursor)


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

            self.get_fit_result() 
            self.update_raw_fit_list()
            self.plot_from_itemlist()
            self.plot3D()
            self.update_xy_combos() 

            #Change plot type combobox text, write into statusbar, change headline for listwidget
            #and enables the copy fit into data analysis and copy fit to file buttons. 
            set_idx = self.plot_type_combo.findText('Fitted')
            self.plot_type_combo.setCurrentIndex(set_idx)
            self.save_fit_X_btn.setEnabled(True)
            self.copy_fit_btn.setEnabled(True)
            self.headline.setText("Fitted/3D plot: Hide/show raw data and fitted lines")
            combo_idx = self.plot_type_combo.findText('Fitted (AC)')
            self.plot_type_combo.setCurrentIndex(combo_idx)

            self.parent.tau_table.update_table()
        
        QApplication.restoreOverrideCursor()

    def get_fit_result(self): 
        """ Fits χ' and χ'' and returns the result of the fit. Fit result contains T, 
        χs, χT, τ, alpha, residual, tau fit error and uncertainty on tau"""

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
        
        self.fit_result = pd.DataFrame(data={'Temp': T,
                                        'Tau': tau,
                                        'dTau': tau_err_RC(tau, tau_fit_err, alpha),
                                        'Alpha': alpha,
                                        'ChiS': Xs,
                                        'ChiT': Xt,
                                        'Residual': resid,
                                        'Tau_Err': tau_fit_err})

        return self.fit_result

    def save_fit_to_file(self):
        """ The fit of χ' and χ'' is saved to a file with columns: 
        Temp, Tau, dTay, Alpha, ChiS, ChiT, Residual and Tau_err. 
        This file can be loaded in the data analysis tab. """
        
        name = QFileDialog.getSaveFileName(self, 'Save File')
        filename = name[0]
        if self.fit_result is None:
            MagMessage("There is nothing to save", "There is probably no fit yet...").exec_()
        elif filename=='':
            pass 
        else:
            df_to_save = self.fit_result.copy()
            df_to_save = df_to_save.reindex(columns=['Temp', 'Tau', 'dTau', 'Alpha','ChiS', 'ChiT', 'Residual', 'Fitting error on τ'])
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
        """Plots 3D plot with temperature vs. frequency vs. X''_molar """

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
        fit_resolution = 1000

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
                y_data = self.temp_subsets[row][self.y_label_3D]
                freq_for_plotting = np.logspace(np.log10(min(y_data)), np.log10(max(y_data)), fit_resolution) #Logspace to equally distribute the points on a log-scale
                self.threeD_plot.ax.plot(np.full(fit_resolution,T), 
                                            np.log10(freq_for_plotting),
                                            Xpp_(freq_for_plotting,
                                                  self.fit_result['ChiS'].iloc[row],
                                                  self.fit_result['ChiT'].iloc[row],
                                                  self.fit_result['Tau'].iloc[row],
                                                  self.fit_result['Alpha'].iloc[row]),
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
                                    linewidth=1)

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

    def update_fitted_and_3D_plot(self): 
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
            if "HM" in self.DC_data_types: 
                self.plot_HM() 
            if "XT" in self.DC_data_types: 
                self.plot_XT()
        self.fit_plot.canvas.draw()
        self.threeD_plot.canvas.draw()
        self.HM_plot.canvas.draw() 


    def copy_fit_to_analysis(self):
        """ Copies the fit of tau obtained after fitting χ' and χ'' 
        to the data analysis tab. """
        QApplication.setOverrideCursor(Qt.WaitCursor)

        try:
            self.parent.data_ana.update_datapoints_wgt() 
            self.parent.data_ana.update_plotting()

        except AttributeError:
            QApplication.restoreOverrideCursor() 
            MagMessage("Fitted data does not exist", "Fitted data does not yet exist in the Data Treatment tab").exec_() 
        except TypeError: 
            QApplication.restoreOverrideCursor() 
            MagMessage("Fitted data does not exist", "Fitted data does not yet exist in the Data Treatment tab").exec_() 
        else: 
            QApplication.restoreOverrideCursor() 
            MagMessage("Succes", "Your fit of tau was succesfully copied to the Data Analysis tab").exec_()

    def make_diamag_correction_calculation_dc(self):
        """ Calculates the diamagnetic correction if data is loaded, and 
        inserts these into the dataframe in their molar form: Mp_m, Mpp_m, Xp_m, Xpp_m.
        Afterwards, the temperature subsets, xy_combos and the datatable are updated. """

        if self.raw_df is None:
            # Don't do the calculation, if there is nothing to calculate on
            pass

        else:
            if "H VS M" not in [self.plot_type_combo.itemText(i) for i in range(self.plot_type_combo.count())]: 
                self.plot_type_combo.addItems(["H VS M", "T VS XT"])
            try:
                m_sample = float(self.m_sample)
                M_sample = float(self.M_sample)
                Xd_sample = float(self.Xd_sample)
                #constant_terms = [float(x) for x in self.constant_terms.split(',')]
                #var_am = [float(x) for x in self.var_am.split(',')]

                #assert len(var_am)%2==0
                #paired_terms = [(var_am[n], var_am[n+1]) for n in range(0,len(var_am),2)]

                if Xd_sample == 0:
                    Xd_sample = -6e-7*M_sample
                
            except (ValueError, AssertionError, AttributeError):
                MagMessage('Error', 'Something wrong in "Sample information"\n').exec_()
            else: 
                H = self.raw_df['Magnetic Field (Oe)']
                H0  = self.raw_df.get('AC Amplitude (Oe)', pd.Series([0]*len(self.raw_df)))
                M = self.raw_df["Moment (emu)"]
                
                # Get molar, corrected values from function in process_ac

                M_molar, X_molar = diamag_correction_dc(H, H0, M, m_sample, M_sample)

                self.insert_molar_values_dc(M_molar, X_molar)

                self.update_subsets()
                self.update_xy_combos()
                self.update_raw_fit_list()
                self.plot_XT_HM() 

                self.fit_X_or_get_phi_btn.setText("(3) Get files for fitting in PHI")
                self.fit_X_or_get_phi_btn.clicked.disconnect() 
                self.fit_X_or_get_phi_btn.clicked.connect(self.get_phi_files)
                self.fit_X_or_get_phi_btn.setEnabled(True)
                MagMessage('Diamagnetic correction',
                            'Diamagnetic correction successful!').exec_()

    def plot_XT_HM(self): 
        """Plotting either XT or HM data or both. Also sets the combo_idx to either one of them
        depending on your data. Also set the headline for the listwidget. """
        
        if "HM" in self.DC_data_types: 
            self.plot_HM() 
            combo_idx = self.plot_type_combo.findText('H VS M')
        else: 
            combo_idx = self.plot_type_combo.findText('T VS XT')    
        if "XT" in self.DC_data_types: 
            self.plot_XT()
          
        self.plot_type_combo.setCurrentIndex(combo_idx)

        self.headline.setText("H VS M: Hide/show raw data and fitted lines") 

    def insert_molar_values_dc(self, M_molar, X_molar): 
        """Insert the molar values of M and χ for DC data into the dataframe. """

        if "Moment_m (emu/mol)" in self.raw_df.columns:
            self.raw_df["Moment_m (emu/mol)"] = M_molar
            self.raw_df["Moment_m (Na*muB)"] = self.raw_df["Moment_m (emu/mol)"].multiply(1.791/(10**4))
        else:            
            M_idx = self.raw_df.columns.get_loc('Moment (emu)')
            self.raw_df.insert(M_idx+1, column = "Moment_m (emu/mol)", value = M_molar)
            self.raw_df.insert(M_idx+2, column = "Moment_m (Na*muB)", value = self.raw_df["Moment_m (emu/mol)"].multiply(1.791/(10**4)))
        
        if "XT" in self.DC_data_types: 
        
            if "X_m (emu/(Oe*mol))" in self.raw_df.columns: 
                self.raw_df["X_m (emu/(Oe*mol))"] = X_molar
                self.raw_df["XT_m (emu*K/(Oe*mol))"] = self.raw_df["X_m (emu/(Oe*mol))"]*self.raw_df["Temperature (K)"]
                self.raw_df["XT_m (cm3*K/mol)"] = self.raw_df["XT_m (emu*K/(Oe*mol))"]
            else:
                self.raw_df.insert(M_idx+3, column = "X_m (emu/(Oe*mol))", value=X_molar)   
                self.raw_df.insert(M_idx+4, column = "XT_m (emu*K/(Oe*mol))", value = self.raw_df["X_m (emu/(Oe*mol))"]*self.raw_df["Temperature (K)"])
                self.raw_df.insert(M_idx+5, column = "XT_m (cm3*K/mol)", value = self.raw_df["XT_m (emu*K/(Oe*mol))"])

        
        self.parent.widget_table.update_table()


    def plot_XT(self): 
        """ Plots XT data """

        self.XT_plot.ax.clear() 

        try: 
            for row in range(self.num_meas_Hs):
                
                H = self.meas_Hs[row]
                if self.Hmin == self.Hmax:  
                    rgb = self.parent.temperature_cmap(0.5)
                else:
                    rgb = self.parent.temperature_cmap((H-self.Hmin)/(self.Hmax-self.Hmin))
    
                x_name = "Temperature (K)"
                y_name = "XT_m (emu*K/(Oe*mol))"
                
                H_label = "{:<6.2f} Oe".format(H)

                self.XT_plot.ax.plot(self.XT_subsets[row][x_name],
                                     self.XT_subsets[row][y_name],
                                     marker='o',
                                     mfc='none',
                                     markersize = 2, 
                                     linewidth = 1, 
                                     label = H_label) 
                
            self.XT_plot.ax.set_xlabel(formatlabel(x_name))
            self.XT_plot.ax.set_ylabel(formatlabel(y_name))

            #if self.XT_legend: 
            #    self.XT_plot.ax.legend() 
            self.XT_plot.canvas.draw()    

        except KeyError: #If errors in sample information  
            pass





    def plot_HM(self): 
        """ Plots H vs. M data. Adds a colorbar if more than one temperature subset 
        of H vs. M data is present. """

        self.HM_plot.ax.clear()

        try: 
            for row in range(self.num_meas_temps):
                
                T = self.meas_temps[row]
                if self.Tmin == self.Tmax:  
                    rgb = self.parent.temperature_cmap(0.5)
                else:
                    rgb = self.parent.temperature_cmap((T-self.Tmin)/(self.Tmax-self.Tmin))
    
                x_name = "Magnetic Field (Oe)"
                y_name = "Moment_m (Na*muB)"
                
                T_label = "{:<6.2f} K".format(T)

                item = self.raw_fit_list.item(row)
                itemdict = item.data(32)
                if itemdict['raw']:
                    self.HM_plot.ax.plot(self.HM_subsets[row][x_name],
                                        self.HM_subsets[row][y_name],
                                        marker='o',
                                        color=rgb,
                                        mfc='none',
                                        markersize = 2, 
                                        linestyle = "none",
                                        label = T_label) 
                if itemdict['fit']:       
                    self.HM_plot.ax.plot(self.HM_subsets[row][x_name],
                                        self.HM_subsets[row][y_name],
                                        linestyle = '-',
                                        color=rgb,
                                        mfc='none',
                                        linewidth = 1,
                                        label = T_label) 
    
                
            self.HM_plot.ax.set_xlabel(formatlabel(x_name))
            self.HM_plot.ax.set_ylabel(formatlabel(y_name))

            #if len(self.HM_subsets) > 1: 
            #    norm = mpl.colors.Normalize(vmin=self.Tmin, vmax=self.Tmax)
            #    self.HM_plot.fig.colorbar(
            #        mpl.cm.ScalarMappable(norm=norm,
            #                              cmap=self.parent.temperature_cmap),
            #                                    orientation='horizontal',
            #        cax=self.HM_plot.cax)
            if self.T_legend: 
                self.HM_plot.ax.legend() 
            self.HM_plot.canvas.draw()    

        except KeyError: #If errors in sample information  
            pass



    def update_sample_info_and_make_dia_correction(self): 
        """ Updates the sample information """
        w = SampleInformation(self.parent, self)
        w.exec_()
        
