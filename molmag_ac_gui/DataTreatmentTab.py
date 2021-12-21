#std packages
import os
import json
from importlib.resources import read_text
from multiprocessing import Pool

#third-party packages 
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, to_hex
import matplotlib.ticker as mticker
from PyQt5.QtWidgets import (QFileDialog, QInputDialog, QListWidgetItem, 
                             QWidget, QVBoxLayout, QPushButton, QLabel, 
                             QHBoxLayout, QComboBox, QStackedWidget, 
                             QCheckBox, QLineEdit, QListWidget, QSplitter)
from PyQt5.QtGui import  QDoubleValidator, QColor


#local imports
from .dialogs import PlottingWindow, MagMessage, FitResultPlotStatus
from .utility import read_ppms_file, update_data_names
from .exceptions import FileFormatError
from .process_ac import (diamag_correction, fit_Xp_Xpp_genDebye, tau_err_RC, 
                         Xp_, Xpp_)
from . import data as pkg_static_data


class DataTreatmentTab(QSplitter): 
    def __init__(self, parent): 
        super(DataTreatmentTab,self).__init__() #Før __init__() peger på QWidget
        self.parent = parent
        self.initUI() 
    
    def initUI(self): 
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

        #Creates dataframe with header and raw_df
        self.raw_df = None
        self.raw_df_header = None
        self.num_meas_freqs = 0
        self.num_meas_temps = 0
        self.temp_subsets = []
        self.meas_temps = []
        self.Tmin, self.Tmax = 0,0
        
        self.raw_data_fit = None
        """Constructing the full layout of tab"""
        ## Making the left column (data loading, fitting and visualization controls)
        self.data_layout = QVBoxLayout()
        self.data_loading_wdgt = QWidget()

        self.load_btn = QPushButton('(1) Load')
        self.load_btn.clicked.connect(self.load_ppms_data)
        self.data_layout.addWidget(self.load_btn)
        
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
        
        ### Constructing data plotting layout
        self.plot_lo = QVBoxLayout()
        
        self.plot_type_header = QLabel('Axis content')
        self.plot_type_header.setFont(self.parent.headline_font)
        self.plot_lo.addWidget(self.plot_type_header)
        
        self.plot_type_combo = QComboBox()
        self.plot_type_combo.addItems(['Raw data', 'Fitted', 'Temp VS Freq VS Xpp'])
        self.plot_type_combo.currentIndexChanged.connect(self.switch_analysis_view)
        self.plot_lo.addWidget(self.plot_type_combo)
        
        self.raw_data_plot_header = QLabel('Raw data plotting')
        self.raw_data_plot_header.setFont(self.parent.headline_font)
        self.plot_lo.addWidget(self.raw_data_plot_header)
        
        ## Constructing the x combobox
        self.x_lo = QHBoxLayout()
        self.x_combo_lbl = QLabel('x')
        self.x_lo.addWidget(self.x_combo_lbl)
        
        self.x_combo = QComboBox()
        self.x_combo.currentIndexChanged.connect(self.plot_from_combo)
        self.x_lo.addWidget(self.x_combo)
        self.plot_lo.addLayout(self.x_lo)
        
        ## Constructing the y combobox
        self.y_lo = QHBoxLayout()
        self.y_combo_lbl = QLabel('y')
        self.y_lo.addWidget(self.y_combo_lbl)
        
        self.y_combo = QComboBox()
        self.y_combo.currentIndexChanged.connect(self.plot_from_combo)
        self.y_lo.addWidget(self.y_combo)
        self.plot_lo.addLayout(self.y_lo)
        
        ## Constructing a combobox for plotting fitted data
        self.fitted_header = QLabel('Fit data plotting')
        self.fitted_header.setFont(self.parent.headline_font)
        self.plot_lo.addWidget(self.fitted_header)
        
        self.fit_combo = QComboBox()
        self.fit_combo.addItems(['ColeCole', 'FreqVSXp', 'FreqVSXpp'])
        self.fit_combo.currentIndexChanged.connect(self.plot_from_itemlist)
        self.plot_lo.addWidget(self.fit_combo)
        
        ## Checkbox for coloring measured points
        self.fit_color_cb_lo = QHBoxLayout()
        self.fit_color_cb_header = QLabel('Use temp. colors')
        self.fit_color_cb_lo.addWidget(self.fit_color_cb_header)
        
        self.fit_data_color_cb = QCheckBox()
        self.fit_data_color_cb.stateChanged.connect(self.plot_from_itemlist)
        self.fit_color_cb_lo.addWidget(self.fit_data_color_cb)
        self.plot_lo.addLayout(self.fit_color_cb_lo)

        ## Finalizing the raw data layout
        self.data_layout.addLayout(self.plot_lo)
        
        ### Finalizing the data loading widget
        self.data_layout.addStretch()
        self.data_loading_wdgt.setLayout(self.data_layout)
        
        ## Making the middle of the tab (data visualization)
        self.sw = QStackedWidget()
        
        self.raw_plot = PlottingWindow()
        self.fit_plot = PlottingWindow(make_ax="cax")
        self.threeD_plot = PlottingWindow(make_ax = "z") 

        self.sw.addWidget(self.raw_plot)
        self.sw.addWidget(self.fit_plot)
        self.sw.addWidget(self.threeD_plot)

        ### Making the right column (parameter controls)
        self.param_wdgt = QWidget()
        self.param_layout = QVBoxLayout()
        
        ## SAMPLE INFO
        self.sample_info_layout = QVBoxLayout()
        
        self.sample_info_header = QLabel('Sample information')
        self.sample_info_header.setFont(self.parent.headline_font)
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
        self.fit_headline = QLabel('Fitted parameters')
        self.fit_headline.setFont(self.parent.headline_font)
        self.param_layout.addWidget(self.fit_headline)
        
        self.raw_fit_list = QListWidget()
        self.param_layout.addWidget(self.raw_fit_list)
        self.raw_fit_list.doubleClicked.connect(self.update_raw_plot)
        
        ## Finalizing layout
        
        self.param_layout.addStretch()
        self.param_wdgt.setLayout(self.param_layout)
        
        self.addWidget(self.data_loading_wdgt)
        self.addWidget(self.sw)
        self.addWidget(self.param_wdgt)

        self.show()

        """End of layout construction"""

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

    def update_meas_temps(self):
        
        meas_temps = []
        for sub in self.temp_subsets:
            meas_temps.append(sub['Temperature (K)'].mean())
        
        self.meas_temps = np.array(meas_temps)
        self.num_meas_temps = len(self.meas_temps)
        self.Tmin = self.meas_temps.min()
        self.Tmax = self.meas_temps.max()

    def update_analysis_combos(self):
    
        self.x_combo.clear()
        self.x_combo.addItems(self.raw_df.columns)
        
        self.y_combo.clear()
        self.y_combo.addItems(self.raw_df.columns)

    def load_ppms_data(self):
        
        #open_file_dialog = QFileDialog()
        #filename_info = open_file_dialog.getOpenFileName(self, 'Open file', self.last_loaded_file)
        #print("filename_info = ", filename_info)
        filename_info =  ('C:/Users/au592011/OneDrive - Aarhus Universitet/Skrivebord/TestData_MAG/ac-data/ac-data/dy-dbm/20180209DyII_1000.dat', 'All Files (*)')
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
            self.parent.current_file = filename
            self.raw_df = potential_df
            self.raw_df_header = potential_header
            
            # Clear old data and set new names
            self.raw_fit_list.clear()
            self.fill_df_data_values()
            
            self.cleanup_loaded_ppms()
            self.num_meas_freqs = len(set(self.raw_df['AC Frequency (Hz)']))
            self.update_temp_subsets()
            self.update_meas_temps()
            
            # Clearing axes of "old" drawings and setting front widget to the raw data
            self.raw_plot.clear_canvas()
            self.fit_plot.clear_canvas()
            self.fit_plot.cax.clear()
            
            combo_idx = self.plot_type_combo.findText('Raw data')
            self.plot_type_combo.setCurrentIndex(combo_idx)
            
            # Updating analysis combos, which will automatically draw the new data
            self.update_analysis_combos()
            
            #Updates "Table of Data" tab with the loaded data
            self.parent.widget_table.updatetable()

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
                self.update_analysis_combos()
                self.parent.widget_table.updatetable()
                MagMessage('Diamagnetic correction',
                               'Diamagnetic correction successful!').exec_()

    def load_sample_data(self):
        #Set like this so I dont have to open it each time: 
        #filename_info = QFileDialog().getOpenFileName(self, 'Open file', self.last_loaded_file)
        filename_info =  ('C:/Users/au592011/OneDrive - Aarhus Universitet/Skrivebord/TestData_MAG/ac-data/ac-data/dy-dbm/dbm_sample_data.dat', 'All Files (*)')
        
        #print("filename_info = ", filename_info)
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
    
    def update_raw_fit_list(self):       
        self.raw_fit_list.clear()
        for i in range(self.num_meas_temps):
            T = self.meas_temps[i]
            newitem = QListWidgetItem()
            newitem.setText('T = {:<6.2f} K, Show raw data: {}, Show fit: {}'.format(round(T,2),True, True)) #Text for Fitted Parameters box
            plotting_dict = {'temp': self.meas_temps[i],
                             'raw': True,
                             'fit': True}
            newitem.setData(32, plotting_dict)
            t_float = (T-self.Tmin)/(self.Tmax-self.Tmin) #Makes t_float scaled by temperature from 0 to 1 for colormap
            newitem.setBackground(QColor(to_hex(self.temperature_cmap(t_float))))
            self.raw_fit_list.addItem(newitem)

    def plot_from_itemlist(self):
        
        if self.raw_fit_list.count()==0:
            return
        
        self.fit_plot.ax.clear()
        plot_type = self.fit_combo.currentText()
        
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
                
            item = self.raw_fit_list.item(row)
            itemdict = item.data(32)
            if itemdict['raw']:
                self.fit_plot.ax.plot(self.temp_subsets[row][x_name],
                                            self.temp_subsets[row][y_name],
                                            marker='o',
                                            mec=markercolor,
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
        self.fit_plot.ax.set_xlabel(x_name)
        self.fit_plot.ax.set_ylabel(y_name)

        norm = mpl.colors.Normalize(vmin=self.Tmin, vmax=self.Tmax)
        self.fit_plot.fig.colorbar(
            mpl.cm.ScalarMappable(norm=norm,
                                  cmap=self.temperature_cmap),
                                        orientation='horizontal',
            cax=self.fit_plot.cax)
        
        self.fit_plot.canvas.draw()

    def fit_Xp_Xpp_standalone(self):
        
        try:
            # Check to see if there has been loaded a data frame
            self.raw_df.columns
            # Check to see if the data to work on is in the data frame
            assert 'Xp_m (emu/(Oe*mol))' in self.raw_df.columns
            
        except AttributeError:
            print("There hasn't been loaded a data frame to work on")
        except AssertionError:
            MagMessage('Error','Calculate diamagnetic correction first\nto make Xp_m and Xpp_m for the algorithm').exec_()
        
        else:
            self.parent.statusBar.showMessage('Running fit...')
            
            # This can't be used currently. Will only work if a separate thread is spawned for fitting.
            #w = QMessageBox()
            #w.setText('Running the fit...\nPlease wait!')
            #w.exec_()
            
            T = [x for x in self.meas_temps]
            Xs, Xt, tau, alpha, resid, tau_fit_err = [],[],[],[],[],[]
            
            v_all = [np.array(self.temp_subsets[t_idx]['AC Frequency (Hz)']) for t_idx in range(self.num_meas_temps)]
            Xp_all = [np.array(self.temp_subsets[t_idx]['Xp_m (emu/(Oe*mol))']) for t_idx in range(self.num_meas_temps)]
            Xpp_all = [np.array(self.temp_subsets[t_idx]['Xpp_m (emu/(Oe*mol))']) for t_idx in range(self.num_meas_temps)]
            
            inputs = tuple(zip(v_all, Xp_all, Xpp_all))
            
            with Pool() as pool:
                res = pool.starmap(fit_Xp_Xpp_genDebye, inputs)
            
            #w.close()
            
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
            set_idx = self.plot_type_combo.findText('Fitted')
            self.plot_type_combo.setCurrentIndex(set_idx)
            
            self.parent.statusBar.showMessage("Fit of X' and X'' complete")

    def save_fit_to_file(self):
        
        name = QFileDialog.getSaveFileName(self, 'Save File')
        filename = name[0]
        if self.raw_data_fit is None:
            MagMessage("There is nothing to save", "There is probably no fit yet...").exec_()
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

    def plot3D(self):

        if self.raw_fit_list.count()==0:
            return

        self.threeD_plot.ax.clear()

        x_label = 'Temperature (K)'
        y_label = "AC Frequency (Hz)"
        z_label = "Xpp_m (emu/(Oe*mol))"
    
                     
        def log_tick_formatter(val, pos=None):
            return f"$10^{{{int(val)}}}$" 
        
        self.threeD_plot.ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
        self.threeD_plot.ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        self.threeD_plot.ax.set_xlabel(x_label)
        self.threeD_plot.ax.set_ylabel(y_label)
        self.threeD_plot.ax.set_zlabel(z_label)            
        
        for row in range(self.num_meas_temps):      
            T = self.meas_temps[row]
            rgb = self.temperature_cmap((T-self.Tmin)/(self.Tmax-self.Tmin))           
                         
            item = self.raw_fit_list.item(row)
            itemdict = item.data(32)
            if itemdict['raw']:
                self.threeD_plot.ax.scatter3D(self.temp_subsets[row][x_label],
                                        np.log10(self.temp_subsets[row][y_label]),
                                        self.temp_subsets[row][z_label], 
                                        color = rgb
                                        )

            if itemdict['fit']:       
                self.threeD_plot.ax.plot(self.temp_subsets[row][x_label], 
                                            np.log10(self.temp_subsets[row][y_label]),
                                            Xpp_(self.temp_subsets[row]['AC Frequency (Hz)'],
                                                  self.raw_data_fit['ChiS'].iloc[row],
                                                  self.raw_data_fit['ChiT'].iloc[row],
                                                  self.raw_data_fit['Tau'].iloc[row],
                                                  self.raw_data_fit['Alpha'].iloc[row]),
                                            c=rgb)  

        norm = mpl.colors.Normalize(vmin=self.Tmin, vmax=self.Tmax)
        self.threeD_plot.fig.colorbar(
            mpl.cm.ScalarMappable(norm=norm,
                                  cmap=self.temperature_cmap),
                                        orientation='horizontal',
            cax=self.threeD_plot.cax)



        self.threeD_plot.canvas.draw()

    def switch_analysis_view(self):

        idx = self.plot_type_combo.currentIndex()
        self.sw.setCurrentIndex(idx)
        if idx == 2: 
            self.plot3D() 

    def plot_from_combo(self):
        
        self.raw_plot.ax.clear()
        
        idx_x = self.x_combo.currentIndex()
        idx_y = self.y_combo.currentIndex()
        x_label = self.raw_df.columns[idx_x]
        y_label = self.raw_df.columns[idx_y]
        
        self.raw_plot.ax.plot(self.raw_df[x_label],
                                    self.raw_df[y_label],
                                    marker='o',
                                    mfc='none',
                                    mec='k',
                                    linestyle='-',
                                    c='k',
                                    linewidth=1,
                                    )
        self.raw_plot.ax.set_xlabel(x_label)
        self.raw_plot.ax.set_ylabel(y_label)
        self.raw_plot.canvas.draw()

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

    def update_itemdict(self, item, itemdict):
        
        item.setData(32, itemdict)
        item.setText('T = {:<6.2f} K, Show raw data: {}, Show fit: {}'.format(round(itemdict['temp'],2),
                                         itemdict['raw'],
                                         itemdict['fit']))

    def update_raw_plot(self):
        
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
            
        self.plot_from_itemlist()
        self.plot3D() 
        self.fit_plot.canvas.draw()

    def copy_fit_to_analysis(self):
    
        try:
            D = np.array(list(zip(self.meas_temps,
                                  self.raw_data_fit['Tau'],
                                  self.raw_data_fit['dTau'])))
            self.parent.data_ana.set_new_t_tau(D)
            self.parent.data_ana.read_indices_for_used_temps()
            self.parent.data_ana.plot_t_tau_on_axes()
            self.parent.data_ana.ana_plot.reset_axes()
        except TypeError:
            print('When the fitted data does not yet exist')
        