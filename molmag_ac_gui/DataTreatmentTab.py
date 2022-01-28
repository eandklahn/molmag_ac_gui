#std packages
import os
from importlib.resources import read_text
from multiprocessing import Pool

#third-party packages 
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import to_hex
import matplotlib.ticker as mticker
from PyQt5.QtWidgets import (QFileDialog, QFormLayout, QInputDialog, QListWidgetItem, 
                             QWidget, QVBoxLayout, QPushButton, QLabel, 
                             QHBoxLayout, QComboBox, QStackedWidget, 
                             QCheckBox, QLineEdit, QListWidget, QSplitter,
                             QGridLayout)
from PyQt5.QtGui import  QDoubleValidator, QColor

from molmag_ac_gui.layout import make_btn, make_headline, make_line


#local imports
from .dialogs import PlottingWindow, MagMessage, FitResultPlotStatus, SampleInformation
from .utility import read_ppms_file, update_data_names, formatlabel
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
        #Creates dataframe with header and raw_df
        self.initialize_attributes()

        """Constructing the full layout of tab"""
        ## Making the left column (data loading, fitting, visualization controls and show/hide option)
        self.layout = QVBoxLayout()
        self.options_wdgt = QWidget()
        
        # Making data import and treatment buttons
        make_headline(self, "Data import and treatment", self.layout)        
        self.add_data_treatment_buttons() 

        # Making Axis content combobox
        make_headline(self, "Axis content", self.layout)
        self.add_plot_type_combobox()

        # Constructing the x and y comboboxes
        make_headline(self, "Raw data plotting", self.layout)
        self.add_xy_comboboxes() 

        # Constructing a combobox for plotting fitted data
        make_headline(self, "Plot type for fitting", self.layout)
        self.add_fit_combobox() 

        # Making the right column (parameter controls)
        self.add_parameter_controls() 
        
        # Finalizing the data loading widget
        self.layout.addStretch()
        self.options_wdgt.setLayout(self.layout)
        self.addWidget(self.options_wdgt)

        # Making data visualization stacking widget
        self.add_stackwidget()

        ## Finalizing layout
        self.setSizes([1,1000])
        self.show()

        """End of layout construction"""

    def initialize_attributes(self): 

        self.raw_df = None
        self.raw_df_header = None
        self.raw_df_reduced = None
        self.num_meas_freqs = 0
        self.num_meas_temps = 0
        self.temp_subsets = []
        self.meas_temps = []
        self.Tmin, self.Tmax = 0,0
        self.raw_data_fit = None

    def add_plot_type_combobox(self): 

        self.plot_type_combo = QComboBox()
        self.plot_type_combo.addItems(['Raw data', 'Fitted', 'Temp VS Freq VS Xpp'])
        self.plot_type_combo.currentIndexChanged.connect(self.switch_analysis_view)
        self.layout.addWidget(self.plot_type_combo)

    def add_fit_combobox(self): 

        self.fit_combo = QComboBox()
        self.fit_combo.addItems(['Cole-Cole', 'Freq VS Xp', 'Freq VS Xpp'])
        self.fit_combo.currentIndexChanged.connect(self.plot_from_itemlist)
        self.layout.addWidget(self.fit_combo)
        
        self.fit_color_cb_lo = QHBoxLayout()

        make_line(self, "Use temperature colors", self.fit_color_cb_lo)
        
        self.fit_data_color_cb = QCheckBox()
        self.fit_data_color_cb.stateChanged.connect(self.plot_from_itemlist)
        self.fit_color_cb_lo.addWidget(self.fit_data_color_cb)
        self.layout.addLayout(self.fit_color_cb_lo)

    def add_data_treatment_buttons(self): 

        make_btn(self, "(1) Load datafile", self.load_ppms_data, self.layout)
        make_btn(self, "(2) Load sample information", self.update_sample_info, self.layout)
        make_btn(self, "(3) Diamagnetic correction", self.make_diamag_correction_calculation, self.layout)
        make_btn(self, "(4) Fit X', X''", self.fit_Xp_Xpp_standalone, self.layout)
        
        #Copy fit to data analysis and save fit to file layout
        self.copy_save_layout = QHBoxLayout()

        make_btn(self, "Copy fit to Data analysis", self.copy_fit_to_analysis, self.copy_save_layout)
        make_btn(self, "Save fit to file", self.save_fit_to_file, self.copy_save_layout)

        self.layout.addLayout(self.copy_save_layout)

    def add_xy_comboboxes(self): 
        self.xy_lo = QFormLayout() 

        self.x_combo = QComboBox()
        self.x_combo.currentIndexChanged.connect(self.plot_from_combo)        
        self.y_combo = QComboBox()
        self.y_combo.currentIndexChanged.connect(self.plot_from_combo)
        
        self.xy_lo.addRow("x: ", self.x_combo) 
        self.xy_lo.addRow("y: ", self.y_combo) 

        self.layout.addLayout(self.xy_lo)

        self.xy_options_lo = QHBoxLayout()

        make_line(self, "Show additional options for x and y ", self.xy_options_lo)
        
        self.xy_options_cb = QCheckBox()
        self.xy_options_cb.stateChanged.connect(self.update_analysis_combos)
        self.xy_options_lo.addWidget(self.xy_options_cb)
        self.layout.addLayout(self.xy_options_lo)

    def add_stackwidget(self): 
        
        self.sw = QStackedWidget()
        
        self.raw_plot = PlottingWindow()
        self.fit_plot = PlottingWindow(make_ax="cax")
        self.threeD_plot = PlottingWindow(make_ax = "z") 

        self.sw.addWidget(self.raw_plot)
        self.sw.addWidget(self.fit_plot)
        self.sw.addWidget(self.threeD_plot)

        self.addWidget(self.sw)


    def add_parameter_controls(self): 
            
        make_headline(self, "Hide/show fitted lines and raw data", self.layout)
        make_line(self, "Double click to edit list", self.layout)
        
        #List of fitted raw data
        self.raw_fit_list = QListWidget()
        self.layout.addWidget(self.raw_fit_list)
        self.raw_fit_list.doubleClicked.connect(self.update_raw_plot)
        

    
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
        self.y_combo.clear()

        chosenlabels = ["Temperature (K)", "AC Frequency (Hz)", "AC Amplitude (Oe)", "Magnetic Field (Oe)", "Mp (emu)",\
                           "Mpp (emu)", "Xp (emu/Oe)", "Xpp (emu/Oe)"]
        
        molarlabels = ["Mp_m (emu/mol)", "Mpp_m (emu/mol)", "Xp_m (emu/(Oe*mol))", "Xpp_m (emu/(Oe*mol))"]
        
        if not self.xy_options_cb.isChecked():
            if all([label in self.raw_df for label in chosenlabels + molarlabels]): #if molar properties have been calculated
                self.x_combo.addItems(chosenlabels + molarlabels)
                self.y_combo.addItems(chosenlabels + molarlabels)
                self.update_reduced_df() 
            elif all([label in self.raw_df for label in chosenlabels]): #If no molar properties has been calculated
                self.update_reduced_df() 
                self.x_combo.addItems(chosenlabels)
                self.y_combo.addItems(chosenlabels)
            else: 
                self.y_combo.addItems(self.raw_df.columns)
                self.x_combo.addItems(self.raw_df.columns)
        else:  #If all properties in chosen labels are not there, show all labels from input file 
            self.y_combo.addItems(self.raw_df.columns)
            self.x_combo.addItems(self.raw_df.columns)
        


    def update_reduced_df(self): 
        chosenlabels = ["Temperature (K)", "AC Frequency (Hz)", "AC Amplitude (Oe)", "Magnetic Field (Oe)", "Mp (emu)",\
                           "Mpp (emu)", "Xp (emu/Oe)", "Xpp (emu/Oe)"]        
        molarlabels = ["Mp_m (emu/mol)", "Mpp_m (emu/mol)", "Xp_m (emu/(Oe*mol))", "Xpp_m (emu/(Oe*mol))"]
        
        if not self.xy_options_cb.isChecked():
            if all([label in self.raw_df for label in chosenlabels + molarlabels]): #if molar properties have been calculated
                self.raw_df_reduced = self.raw_df[chosenlabels + molarlabels]
            elif all([label in self.raw_df for label in chosenlabels]): #If no molar properties has been calculated
                self.raw_df_reduced = self.raw_df[chosenlabels]

    def load_ppms_data(self):
        
        open_file_dialog = QFileDialog()
        filename_info = open_file_dialog.getOpenFileName(self, 'Open file', self.parent.last_loaded_file)
        #filename_info = ('C:/Users/au592011/OneDrive - Aarhus Universitet/Skrivebord/TestData_MAG/ac-data/ac-data/dy-dbm/20180209DyII_1000.dat', 'All Files (*)') 
        filename = filename_info[0]
        try:
            # FileNotFoundError and UnicodeDecodeError will be raised here
            potential_header, potential_df = read_ppms_file(filename)
            if potential_header is None:
                raise FileFormatError(filename)
            summary = update_data_names(potential_df, self.parent.read_options)
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
            print('A data name from self.parent.read_options is showing up more than once in the columns')
            print('OR')
            print('both Mp and Xp are unexpectedly both showing up in the data names')
        
        else:
            # Now that everything has been seen to work,
            # save potential header and potential df as actual header and df
            self.parent.last_loaded_file = os.path.split(filename)[0]
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
                self.update_analysis_combos()
                self.parent.widget_table.updatetable()
                MagMessage('Diamagnetic correction',
                               'Diamagnetic correction successful!').exec_()


    
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
            newitem.setBackground(QColor(to_hex(self.parent.temperature_cmap(t_float))))
            self.raw_fit_list.addItem(newitem)

    def plot_from_itemlist(self):
        
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
            markercolor = 'k'
            if self.fit_data_color_cb.isChecked():
                markercolor = rgb
            
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
        self.fit_plot.ax.set_xlabel(formatlabel(x_name))
        self.fit_plot.ax.set_ylabel(formatlabel(y_name))

        norm = mpl.colors.Normalize(vmin=self.Tmin, vmax=self.Tmax)
        self.fit_plot.fig.colorbar(
            mpl.cm.ScalarMappable(norm=norm,
                                  cmap=self.parent.temperature_cmap),
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
            self.plot3D()
            self.update_analysis_combos() 
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
        self.threeD_plot.ax.set_xlabel(formatlabel(x_label))
        self.threeD_plot.ax.set_ylabel(formatlabel(y_label))
        self.threeD_plot.ax.set_zlabel(formatlabel(z_label))             
        
        for row in range(self.num_meas_temps):      
            T = self.meas_temps[row]
            rgb = self.parent.temperature_cmap((T-self.Tmin)/(self.Tmax-self.Tmin))           
                         
            item = self.raw_fit_list.item(row)
            itemdict = item.data(32)
            if itemdict['raw']:
                self.threeD_plot.ax.scatter3D(self.temp_subsets[row][x_label],
                                        np.log10(self.temp_subsets[row][y_label]),
                                        self.temp_subsets[row][z_label], 
                                        color = rgb, s = 7
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
                                  cmap=self.parent.temperature_cmap),
                                        orientation='horizontal',
            cax=self.threeD_plot.cax)



        self.threeD_plot.canvas.draw()

    def switch_analysis_view(self):

        idx = self.plot_type_combo.currentIndex()
        self.sw.setCurrentIndex(idx)

    def getxylabel(self): 
        idx_x = self.x_combo.currentIndex()
        idx_y = self.y_combo.currentIndex()
        if self.xy_options_cb.isChecked(): #If extra xy_options are chosen with checkbox: Index in raw_df
            x_label = self.raw_df.columns[idx_x]
            y_label = self.raw_df.columns[idx_y]
        else: #If limited xy_options
            if idx_x < len(self.raw_df_reduced.columns): 
                x_label = self.raw_df_reduced.columns[idx_x]
            else: 
                x_label = self.raw_df_reduced.columns[0]
            if idx_y < len(self.raw_df_reduced.columns): 
                y_label = self.raw_df_reduced.columns[idx_y]
            else: 
                y_label = self.raw_df_reduced.columns[0]
        return x_label, y_label 

    def plot_from_combo(self):
        
        self.raw_plot.ax.clear()
        
        x_label, y_label = self.getxylabel()
        
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
            self.parent.data_ana.plot_wdgt.reset_axes()
        except TypeError:
            MagMessage("Fitted data does not exist", "Fitted data does not yet exist in the Data Treatment tab").exec_() 

    def update_sample_info(self): 
        w = SampleInformation(self.parent)
        w.exec_()
