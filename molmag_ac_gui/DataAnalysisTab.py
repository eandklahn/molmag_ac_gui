#std packages 
import os
from re import X
import sys
from collections import deque
import datetime

#third-party packages
import numpy as np
import pandas as pd
import names
from matplotlib.colors import to_hex
from matplotlib._color_data import TABLEAU_COLORS
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QDoubleSpinBox, QFileDialog, QListWidgetItem, 
                             QMessageBox, QWidget, QVBoxLayout, 
                             QLabel, QHBoxLayout, QApplication,
                             QListWidget, QSplitter, QComboBox, QTextEdit, QDialog)


#local imports
from .exceptions import NoGuessExistsError
from .process_ac import (getParameterGuesses, fit_relaxation,
                         default_parameters, add_partial_model)
from .dialogs import (GuessDialog, MagMessage, 
                      PlottingWindow, SimulationDialog)
from .layout import make_checkbox, make_headline, make_btn, make_line, headline_font
from lmfit import fit_report


class DataAnalysisTab(QSplitter): 
    def __init__(self, parent): 
        super(DataAnalysisTab,self).__init__() 
        self.parent = parent
        self.initUI() 
    
    def initUI(self): 

        # Initializes simulations colors and other attributes
        self.initialize_attributes() 

        # The options widget to the left 
        self.options_wdgt = QWidget()
        self.options_layout = QVBoxLayout()
        
        # Adding data loading options
        make_headline(self, "Data loading options", self.options_layout)

        self.add_data_load_options() 

        #Add options to remove outliers
        make_headline(self, "Remove outliers from plotting/fitting", self.options_layout) 
        self.add_outliers_controls() 
        
        # Adding fit controls with checkboxes
        make_headline(self, "Fitting options", self.options_layout)
        self.add_fit_type_checkboxes()

        # Adding temperature controls
        self.add_temp_controls() 

        # Adding a button to run a fit
        make_btn(self, "Run fit!", self.make_the_fit, self.options_layout)

        # Adding a list to hold information about simulations
        make_headline(self, "Plot/simulate fit", self.options_layout)
        self.add_simulations_list()

        # Adding buttons to control simulation list
        self.sim_btn_layout = QHBoxLayout()
        make_btn(self, "New", self.add_sim, self.sim_btn_layout)
        make_btn(self, "Delete", self.delete_sim, self.sim_btn_layout)
        self.options_layout.addLayout(self.sim_btn_layout)


        #Adding view fitted parameters button 
        make_headline(self, "View information about fits", self.options_layout)        
        self.add_fit_parameters_view()     
          

        #Setting the layout of the options widget
        self.options_wdgt.setLayout(self.options_layout)
        self.options_layout.addStretch() 
        self.addWidget(self.options_wdgt)
        
        #Adding plotting widget
        self.add_plot_wdgt() 

        # Finalizing layout of the data analysis tab
        self.setSizes([1,1200])
        self.show() 


    def add_data_load_options(self): 
        data_load_layout = QHBoxLayout() 
        make_btn(self, "Import current fit", self.copy_data_treat_fit, data_load_layout)
        make_btn(self, 'Load fit from file', self.load_t_tau_data, data_load_layout)
        self.options_layout.addLayout(data_load_layout)


    def copy_data_treat_fit(self): 
        """ Copies the fit of tau from the data treatment tab into the data analysis tab. 
        First it clears all previos fits and lists. Then it makes sure the data is not DC data, if not it imports the fit 
        from the Data treatment tab """
        self.clear_fits() 
        try:
            if self.parent.data_treat.data_type == "DC": 
                MagMessage("Error", "The data loaded in data treatment is DC data, not AC data. \nThe functions in this tab requires AC data.").exec_() 
                return
        except:
            pass
        

        self.parent.data_treat.copy_fit_to_analysis()
        


    def clear_fits(self): 
        """ Clears all fits, unchecks checkboxes, lists etc."""

        self.plot_wdgt.ax.clear() 
        self.list_of_simulations.clear()
        self.fit_history.clear() 

        #Now that the fit_history is empty, the fit combobox should be updated (emptied)
        self.update_fit_combo() 

        #Sets all checkboxes to false
        self.orbach_cb.setChecked(False)
        self.raman_cb.setChecked(False)
        self.qt_cb.setChecked(False)
        #Set temperature values to 0
        self.temp_line[1].setValue(0)
        self.temp_line[3].setValue(0)
        
        # I should find out what these below does - are they important?
        self.fit_stat_txt = ""
        #self.fit_summary.setText("")
        #self.fit_title.setText("")

    def add_fit_parameters_view(self): 
        """ Adds fit parameters view, where the user can view information about the
        fits that have been made. Also adds an option to save the fit in a file. """

        self.choose_fit_layout = QHBoxLayout() 
        self.choose_fit_combo = QComboBox()
        self.choose_fit_combo.currentIndexChanged.connect(
                                                  self.show_fit_stat)

        self.choose_fit_line = [QLabel('Choose fit from list: '), self.choose_fit_combo]
        for w in self.choose_fit_line: 
            self.choose_fit_layout.addWidget(w)
        self.options_layout.addLayout(self.choose_fit_layout)

        self.fit_title = QLabel()
        self.fit_title.setFont(headline_font)
        self.options_layout.addWidget(self.fit_title)

        self.fit_summary = QTextEdit()
        self.fit_summary.setReadOnly(True)
        self.options_layout.addWidget(self.fit_summary)

        make_btn(self, "Save fit statistics to file", self.save_fit_statistics, self.options_layout)
    
    def save_fit_statistics(self):
        """ Saves the fit statistics into a file """ 

        name = QFileDialog.getSaveFileName(self, 'Save File')
        filename = name[0]

        additional_notes = "Notes on the printed fit statistics:\n Parameters that are fixed, are not used in the fitting. \n \
Their contributions to tau are basically multiplied by 0 (using the useQT, useR etc. parameters) \n \
t0 has been scaled by a factor of 10^7, such that t0_scaled = t0*10^7.\n \
Ueff is scaled by kB, such that Ueff_by_kB = Ueff/kB.\n \
This is ensure more stable fitting by having similar magnitudes of parameters.\n\n"

        name, ext = os.path.splitext(filename)
        if ext == '':
            ext = '.txt'
        
        try: 
            self.set_fit_stat_txt(all_significant_digits=True) 
            with open(filename + ext, "w") as f:
                f.write(additional_notes)
                f.write(self.fit_stat_txt)
            MagMessage("Succes", "File succesfully written").exec_() 
       
        except IndexError: #If fit history is empty 
            MagMessage("Error", "No fit made yet").exec_() 
       
        except ValueError: #If fit_result format is different, it will save in original format
            fit_idx = self.choose_fit_combo.currentIndex() 
            name, res, time, Tmin, Tmax = self.fit_history[fit_idx]
            self.fit_stat_txt = fit_report(res)
            with open(filename + ext, "w") as f:
                f.write(additional_notes)
                f.write(self.fit_stat_txt)
            MagMessage("Succes", "File succesfully written").exec_() 



    def set_fit_stat_txt(self, all_significant_digits = False): 
        """ Formatting the result from the fit_report and saving the formatted 
        text in self.fit_stat_txt. Displays them with only few siginificant digits """

        def printstring(line, nb_before, nb_float, nb_after,):
            if float(line[nb_float]) > 0.1 and float(line[nb_float]) < 100: 
                string  = "     " + "{}"*nb_before + "{}" + "{}"*nb_after + "\n"
                self.fit_stat_txt += string.format(*line[:nb_float],float(line[nb_float]), *line[nb_float+1:])
            elif all_significant_digits: 
                string = "     " + "{} "*nb_before + "{:e} " + "{} "*nb_after + "\n"
                self.fit_stat_txt += string.format(*line[:nb_float],float(line[nb_float]), *line[nb_float+1:])
            else: 
                string = "     " + "{} "*nb_before + "{:.2e} " + "{} "*nb_after + "\n"
                self.fit_stat_txt += string.format(*line[:nb_float],float(line[nb_float]), *line[nb_float+1:])


        self.fit_stat_txt = "" 

        fit_idx = self.choose_fit_combo.currentIndex() 
        _, res, _, _, _ = self.fit_history[fit_idx]
        joined_result = "".join(fit_report(res))
        lines = [line.split() for line in joined_result.split("\n")]
        
        
        for line in lines: 
            if len(line)>=1 and line[0][:2] == "[[":
                line = [word.replace("[","").replace("]","") for word in line]
                self.fit_stat_txt += " ".join(line) + "\n"
            elif len(line) == 3:
                try: 
                    printstring(line, 2, 2, 0)
                except ValueError: 
                    self.fit_stat_txt += "     "+" ".join(line) + "\n"
            elif len(line) == 4:
                try:
                    printstring(line, 3,3,0)
                except ValueError: 
                    self.fit_stat_txt += "     "+" ".join(line) + "\n"
            elif len(line) == 5:
                try: 
                    printstring(line,4,4,0)
                except ValueError: 
                    self.fit_stat_txt += "     "+" ".join(line) + "\n"
            elif len(line) == 8: 
                try:
                    if all_significant_digits: 
                        self.fit_stat_txt += "     {} {} {} {} {} {} {} {} ) \n".format(line[0], float(line[1]), line[2], float(line[3]), line[4], line[5], line[6], float(line[7][:-1]))
                    elif not all_significant_digits: 
                        self.fit_stat_txt += "     {} {:.2e} {} {:.2e} {} {} {} {:.2e} ) \n".format(line[0], float(line[1]), line[2], float(line[3]), line[4], line[5], line[6], float(line[7][:-1]))
                except ValueError: 
                    self.fit_stat_txt += "     "+" ".join(line) + "\n"
            else: 
                self.fit_stat_txt += "     "+" ".join(line) + "\n"


    def update_fit_combo(self):
        """ Updates the fit combobox with a time and name for each fit in the fit
        history, """

        self.choose_fit_combo.clear() 
        for fit in self.fit_history: 
            name, _, time, _, _ = fit
            self.choose_fit_combo.addItem(f'{time}: {name}')

    def show_fit_stat(self):
        """ Shows and sets the fit statistics text in the parameter view"""

        if self.fit_history == []: 
            return 
        
        fit_idx = self.choose_fit_combo.currentIndex()
        name, res, time, _, _ = self.fit_history[fit_idx]
        title = f'{time}: {name}'
        try: 
            self.set_fit_stat_txt() 
        except (ValueError, IndexError): 
            self.fit_stat_txt = fit_report(res)
        self.fit_summary.setText(self.fit_stat_txt)

        self.fit_title.setText(title)

    def initialize_attributes(self):  
        """ Initializes attributes such as simulation colors, temperature, tau etc."""

        self.simulation_colors = [x for x in TABLEAU_COLORS]
        self.simulation_colors.remove('tab:gray')  
        self.simulation_colors = deque(self.simulation_colors) 
        self.sim_dialog = None

        self.fit_history = list()

        #I am not sure its helpful to use all these, now that we have the dataframe fit_result containing all info
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
        
        self.checked_datapoints = None

        self.fitted_params_dialog = None

    def add_temp_controls(self): 
        """Makes a QHBoxLayout with a two spinboxes where the temperature range is chosen 
        and adds it to the options layout."""

        self.temp_horizontal_layout = QHBoxLayout()

        self.temp_line = [QLabel('Temperature range: '), QDoubleSpinBox(), QLabel('K to'),
                          QDoubleSpinBox(), QLabel('K')]
        
        self.temp_line[1].setRange(0,1000)
        self.temp_line[1].setSingleStep(0.5)
        self.temp_line[3].setRange(0,1000)
        self.temp_line[3].setSingleStep(0.5)

        self.temp_line[1].editingFinished.connect(self.set_new_temp_ranges)
        self.temp_line[3].editingFinished.connect(self.set_new_temp_ranges)

        #Adds all four widgets from the self.temp_lines
        for wgt in self.temp_line:
            self.temp_horizontal_layout.addWidget(wgt)

        self.options_layout.addLayout(self.temp_horizontal_layout)

    

    def delete_items(self):
        QApplication.setOverrideCursor(Qt.WaitCursor)
        df = self.parent.data_treat.fit_result
        if self.checked_datapoints is not None and self.checked_datapoints != []: 
            for data_point in self.checked_datapoints: 
                row_index = df[df.eq(data_point).all(1)].index[0]
                self.datapoints_wgt.takeItem(row_index)
                df = df.drop(row_index)
                df = df.reset_index(drop=True)

            self.parent.data_treat.fit_result = df
            self.checked_datapoints = None
            self.update_datapoints_wgt() 
            self.update_plotting() 
            QApplication.restoreOverrideCursor()
        else: 
            QApplication.restoreOverrideCursor()
            MagMessage("Error", "No datapoints to delete").exec_() 
        
        

    def save_items(self):
        """ Used by the add_outliers_controls function.
        Saves the selected data points in self.checked_datapoints"""

        self.checked_datapoints = []
        for i in range(self.datapoints_wgt.count()):
            item = self.datapoints_wgt.item(i)
            if item.checkState() == Qt.Checked: 
                row_index =  self.datapoints_wgt.row(item)
                self.checked_datapoints.append(self.parent.data_treat.fit_result.values[row_index])

        self.update_plotting()
        
        
    def update_datapoint_colors(self): 
        self.plot_t_tau_on_axes() 
        self.plot_wdgt.canvas.draw()


    def add_outliers_controls(self): 
        """ Adds controls for the data plotting where one can pick which datapoints
        and fits to be visualized. """
        
        make_line(self, "Deletion of data points will also remove them from the datatable tab.", self.options_layout)

        #Creates a QlistWidget showing all datapoints
        self.datapoints_wgt = QListWidget()
        self.options_layout.addWidget(self.datapoints_wgt)
        self.datapoints_wgt.itemChanged.connect(self.save_items)

        #Makes a button where all selected items can be deleted
        make_btn(self, "Delete checked data points", self.delete_items, self.options_layout)
    
    
    def update_datapoints_wgt(self):
        """ Clears the widget containing all the data points. 
        Adds all the data points according to tge current the dataframe of relaxation times """

        self.datapoints_wgt.clear()
        df = self.parent.data_treat.fit_result
        
        for i in range(len(df)): 
            datapoint = df.loc[i]
            T = datapoint['Temp']
            tau = datapoint['Tau']

            item = QListWidgetItem()
            item.setText('T = {:<6.2f} K, tau = {:<.5f} ms'.format(round(T,2),tau*10**3)) #Text for Fitted Parameters box
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable ) #Makes the item checkable
            item.setCheckState(Qt.Unchecked) #Unchecks item

            self.datapoints_wgt.addItem(item)
        
        

    def add_simulations_list(self): 
        """ Makes a QListWidget with all simulations and adds it to the options layout 
        Inspiration (Emil): https://stackoverflow.com/questions/41353653/how-do-i-get-the-checked-items-in-a-qlistview """

        self.list_of_simulations = QListWidget()
        self.list_of_simulations.setDragDropMode(self.list_of_simulations.InternalMove)
        self.list_of_simulations.doubleClicked.connect(self.edit_simulation)
        self.list_of_simulations.itemChanged.connect(self.update_plotting)
        self.options_layout.addWidget(self.list_of_simulations)


    def add_fit_type_checkboxes(self):
        """Creates a QHBoxLayout with three checkboxes and adds these to the options layout"""

        self.fit_type_layout = QHBoxLayout() 
        make_line(self, "Fit types: ", self.fit_type_layout)
        self.orbach_cb = make_checkbox(self, self.read_fit_type_cbs, self.fit_type_layout, "Orbach") 
        self.raman_cb = make_checkbox(self, self.read_fit_type_cbs, self.fit_type_layout, "Raman") 
        self.qt_cb = make_checkbox(self, self.read_fit_type_cbs, self.fit_type_layout, "QT")
        self.direct_cb = make_checkbox(self, self.read_fit_type_cbs, self.fit_type_layout, "Direct")
        self.options_layout.addLayout(self.fit_type_layout)


    def add_plot_wdgt(self): 
        """ Adds a plotting widget for 1/T vs ln(tau) plot """

        self.plot_wdgt = PlottingWindow()
        self.plot_wdgt.ax.set_xlabel('1/T ($K^{-1}$)')
        self.plot_wdgt.ax.set_ylabel(r'$\ln{\tau}/s$') #Maybe this should be ln(τ/s)
        self.addWidget(self.plot_wdgt)



    def reset_analysis_containers(self):
        """ Resets all data analysis containers"""

        self.data_T = None
        self.data_tau = None
        self.data_dtau = None
        
        self.used_T = None
        self.not_used_T = None
        self.used_tau = None
        self.not_used_tau = None
        self.used_dtau = None
        self.not_used_dtau = None
        
        self.checked_datapoints = None
        
        self.used_indices = None


    def read_and_sort_t_tau(self):
        """
        Uses the array D to set new values for T, tau, and alpha

        """
        D = self.parent.data_treat.fit_result.values

        T = D[:,0]
        tau = D[:,1]
        dtau = D[:,2]
        
        sort_indices = T.argsort()

        #Sorts the data according to the temperature (lowest T is the first data point)
        self.data_T = T[sort_indices]
        self.data_tau = tau[sort_indices]
        self.data_dtau = dtau[sort_indices]


    def read_indices_for_used_temps(self):
        """Reads the values for the chosen temperature range, and uses their indicies to 
        update self.used_T, self.used_tau, self.not_used_T and self.not_used_tau  """
        
        min_t = self.temp_line[1].value()
        max_t = self.temp_line[3].value()
        
        try:
            self.used_indices = [list(self.data_T).index(t) for t in self.data_T if t>=min_t and t<=max_t]
            
            self.used_T = self.data_T[self.used_indices]
            self.used_tau = self.data_tau[self.used_indices]
            
            self.not_used_T = np.delete(self.data_T, self.used_indices) #Deletes used indicies from data_T to save in not used T
            self.not_used_tau = np.delete(self.data_tau, self.used_indices) #Deletes used indicies from data_T to save in not used T
            
            
            self.used_dtau = self.data_dtau[self.used_indices]
            self.not_used_dtau = np.delete(self.data_dtau, self.used_indices)
            
        except (AttributeError, TypeError):
            MagMessage("Error", 'No data have been selected yet!').exec_() 

    def plot_t_tau_on_axes(self):
        """ Plots 1/T vs. τ in the plotting window. Also adds errorbars if errors on τ are present """

        self.plot_wdgt.clear_canvas()
        self.plotted_data_pointers = []
        
        #Data points with T within the range, are plotted in blue, the others are plotted in black
        err_used_point, caplines1, barlinecols1 = self.plot_wdgt.ax.errorbar(1/self.used_T,
                                                                            np.log(self.used_tau),
                                                                            yerr=self.used_dtau,
                                                                            fmt='bo',
                                                                            ecolor='b',
                                                                            label='Data',
                                                                            zorder=0.1)
        err_not_used_point, caplines2, barlinecols2 = self.plot_wdgt.ax.errorbar(1/self.not_used_T,
                                                                                np.log(self.not_used_tau),
                                                                                yerr=self.not_used_dtau,
                                                                                fmt='ko',
                                                                                ecolor='k',
                                                                                label='Data',
                                                                                zorder=0.1)

        
        if self.checked_datapoints is not None and self.checked_datapoints != []: 
            for datapoint in self.checked_datapoints: 
                T, tau, dtau = datapoint[0:3]
                self.plot_wdgt.ax.errorbar(1/T, np.log(tau), yerr = dtau, fmt = 'ro', ecolor = 'r', zorder = 0.1)


        self.plotted_data_pointers.append(err_used_point)
        self.plotted_data_pointers.append(err_not_used_point)
        for e in [caplines1, caplines2, barlinecols1, barlinecols2]:
            for line in e:
                self.plotted_data_pointers.append(line)



    def add_T_axis(self): 
        """ Adds Temperature as a second x-axis on top of the 1/T vs. ln(τ) plot """
    
        try: 
            self.plot_wdgt.ax2
        except AttributeError:
            self.plot_wdgt.ax2 = self.plot_wdgt.ax.twiny() 

        self.plot_wdgt.ax2.set_xlabel("Temperature (K)")
        


    def update_T_axis(self):
        """ Updates the second temperature axis (T axis on top of the plot) """
        self.plot_wdgt.ax2.set_xticklabels([])
        
        labels = self.plot_wdgt.ax.get_xticks()
        T_tick_loc = np.array(labels[1:-1])

        def tick_f(X):
            V = []
            for x in X:
                if x != 0:
                    V.append(1 / x)
                else:
                    V.append(0)  # Or any other appropriate value for handling division by zero
            return ["%.2f" % z for z in V]

        self.plot_wdgt.ax2.set_xlim(self.plot_wdgt.ax.get_xlim())
        self.plot_wdgt.ax2.set_xticks(T_tick_loc)
        self.plot_wdgt.ax2.set_xticklabels(tick_f(T_tick_loc))
        self.plot_wdgt.canvas.draw() 








    def load_t_tau_data(self):
        """Load 1/T vs. τ data from file generated in data treatment """

        filename_info = QFileDialog().getOpenFileName(self, 'Open file', self.parent.last_loaded_file)
        #print("filename = ", filename_info[0])
        filename = filename_info[0]
        #filename = "C:/Users/au592011/OneDrive - Aarhus Universitet/Skrivebord/TestData_MAG/ac-data/ac-data/dy-dbm/fakefit.dat"
          
        self.last_loaded_file = os.path.split(filename)[0]
        
        if filename == '':
            pass
        else:
            self.reset_analysis_containers()
        
        try:
            #Need the column names for the dataframe for the tab with relaxation time data. 
            with open(filename, 'r') as f:
                column_names = f.readline().strip().split(';')
            D = np.loadtxt(filename,
                           skiprows=1,
                           delimiter=';')
            
        except (ValueError, OSError) as error_type:
            sys.stdout.flush()
            if error_type == 'ValueError':
                msg = MagMessage("ValueError", 'File format not as expected')
                msg.setIcon(QMessageBox.Warning)
                msg.exec_()
            elif error_type == 'OSError':
                pass
        else:

            self.clear_fits() 
            self.parent.data_treat.fit_result = pd.DataFrame(D, columns=column_names)            
            self.update_datapoints_wgt() 
            self.update_plotting()






    def read_fit_type_cbs(self):
        """ Read the fit type (i.e. which of Quantum Tunneling, Raman and Orbach are 
        included in the fitting). Depends on what is chosen in the comboboxes """
    
        list_of_checked = []
        if self.direct_cb.isChecked(): list_of_checked.append('D')
        if self.qt_cb.isChecked(): list_of_checked.append('QT')
        if self.raman_cb.isChecked(): list_of_checked.append('R')
        if self.orbach_cb.isChecked(): list_of_checked.append('O')

        fitToMake = ''.join(list_of_checked)
        
        return fitToMake

    def set_new_temp_ranges(self):
        """ Sets a new temperature range to fit in, when the spin boxes are changed"""

        self.read_indices_for_used_temps()
        #Changes the color of the data points within the range
        if self.data_T is not None:
            self.update_datapoint_colors() 
        

    def make_the_fit(self):
        """ Makes the fit of 1/T vs. τ"""
        
        #These are used for making a Magmessage pop-up after fitting.
        window_title = 'Fit aborted'
        msg_text = ''
            
        try:
            
            Tmin = self.temp_line[1].value()
            Tmax = self.temp_line[3].value()
            assert Tmin != Tmax
            assert Tmin < Tmax 
            
            #Reads the checked fit types and checks if any has been checked.
            fittypes = self.read_fit_type_cbs()
            assert fittypes != ''

            #Makes initial guesses of the parameters and opens the dialog, where user input can be used as guesses
            guess = getParameterGuesses(self.used_T, self.used_tau, fittypes)
            guess_dialog = GuessDialog(self,
                                       guess,
                                       fittypes)
            accepted = guess_dialog.exec_()
            if not accepted: raise NoGuessExistsError
            

            params = guess_dialog.current_guess
            QApplication.setOverrideCursor(Qt.WaitCursor)

            minimize_res = fit_relaxation(self.used_T, self.used_tau, params)

        except (AssertionError, IndexError):
            msg_text = "Bad temperature or fit settings \nPossible errors: \n \
            - Minimum and maximum temperatures are the same \n \
            - Minimum temperature is larger than maximum temperature \n \
            - No fit options have been selected \n \
            - Can't fit only one data point "
        except RuntimeError:
            msg_text = 'This fit cannot be made within the set temperatures'
        except ValueError:
            msg_text = 'No file has been loaded or there might be some other problem. \nTry to choose another temperature setting or change the fit type, \nif you have already loaded a file.'
        except TypeError:
            msg_text = 'No data has been selected'
        except NoGuessExistsError:
            msg_text = 'Made no guess for initial parameters'
        
        else:
            window_title = 'Fit successful!'
            msg_text = 'Congratulations'
            
            self.add_to_history(minimize_res, fittypes, Tmin, Tmax)
        
            self.update_fit_combo() 
        finally:
            QApplication.restoreOverrideCursor()
            msg = MagMessage(window_title, msg_text)
            msg.setIcon(QMessageBox.Information)
            msg.exec_()  
            if msg_text == 'Congratulations': 
                reply = QMessageBox.question(self, "New simulation obtained!", "Do you want to add the simulation to the plot?",
                                     QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Yes)
                if reply == QMessageBox.Yes:
                    self.add_sim()
        
    
    def add_to_history(self, res, fittype, Tmin, Tmax):
        """ Adds a fit to the top of the fit history"""

        time = datetime.datetime.now()
        time = time.strftime("%d.%m %H:%M:%S")
        self.fit_history.insert(0, (fittype, res, time, Tmin, Tmax))




    def redraw_simulation_lines(self):
        """ Redraw simulation lines """

        for idx in range(self.list_of_simulations.count()):
            item = self.list_of_simulations.item(idx)
            data = item.data(32)
            if item.checkState() == Qt.Checked:
                data['line'].set_visible(True)
            elif item.checkState() == Qt.Unchecked:
                data['line'].set_visible(False)

        for idx in range(self.list_of_simulations.count()):
            item = self.list_of_simulations.item(idx)
            data = item.data(32)
            if data['line'].get_visible():
                self.plot_wdgt.ax.plot(data['line'].get_xdata(), data['line'].get_ydata(), color = data['color'])

        self.plot_wdgt.canvas.draw()

    def get_item_txt(self, T_vals, params):
        """ The string displayed to show information on each fit """
        
        fs = [p for p in params if 'use' in p]
        used = [bool(params[p].value) for p in fs]     

        functions_text = ""
        if used[0]: 
            functions_text += "Quantum Tunneling, "
        if used[1]: 
            functions_text += "Raman, "
        if used[2]: 
            functions_text += "Orbach, "
        
        functions_text = functions_text[:-2]
        functions_text += "."


        text = ['Processes used for fitting: {}\n'.format(functions_text),
                'Temperature range: {} K and {} K\n'.format(*T_vals)]
        return ''.join(text)

    def add_fit_to_list_of_simulations(self, params, T_vals): 
        """Adds a simulation to the list_of_simulations"""

        sim_item = QListWidgetItem()
        sim_item.setFlags( sim_item.flags() | Qt.ItemIsUserCheckable )
        sim_item.setCheckState(Qt.Checked)
            
        self.list_of_simulations.addItem(sim_item)
        color = self.simulation_colors.pop()
        label = names.get_first_name()
        
        line = add_partial_model(self.plot_wdgt.fig,
                                 T_vals[0], 
                                 T_vals[1],
                                 params,
                                 c=color,
                                 label=label)
        
        list_item_data = {'params': params,
                          'T_vals': T_vals,
                          'line': line,
                          'color': color}
        
        new_item_text = self.get_item_txt(T_vals, params)
        
        sim_item.setData(32, list_item_data)
        sim_item.setText(new_item_text)
        sim_item.setBackground(QColor(to_hex(color)))

    def load_new_item(self): 
        """ If new is chosen for the simulation list, the following values are set """
        
        params = default_parameters()
        T_vals = [1,3]
        line = False
        label = None
        color = None

        return params, T_vals, line, color, label 


    def add_sim(self):

        params, T_vals, line, _, _ = self.load_new_item() 
        
        self.sim_dialog = SimulationDialog(parent=self,
                                      fit_history=self.fit_history,
                                      params = params,
                                      min_max_T=T_vals)

        if self.sim_dialog.exec_() == QDialog.Accepted:
            self.add_fit_to_list_of_simulations(params, T_vals) 
            self.update_fit_combo()


    def edit_simulation(self): 
        return 




    def load_chosen_item(self): 
        """ Loads the parameters, T_values, line, color and label of the chosen 
        simulation item """

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

            return sim_item, params, T_vals, line, color, label 


    def set_default_values(self): 
        """ If new is chosen for the simulation list, and no fit is chosen in the 
        fit_history combobox, these will be the default values. """
        
        params = default_parameters()
        T_vals = [1,3]
        line = False
        label = None
        color = None

        return params, T_vals, line, color, label 
    
    def delete_sim(self):
        """ Deleting a simulation/fit from the list"""

        try:
            sim_item = self.list_of_simulations.selectedItems()[0]
        except ValueError as e: 
            print("Value error in delete sim occoured")
            print(e)
        except IndexError: 
            MagMessage("Error", "No simulation has been selected").exec_() 
            return 
        else:
            line_pointer = sim_item.data(32)['line']
            line_color = line_pointer._color
            
            #Deletes the item from the lsit of simulations
            item_row = self.list_of_simulations.row(sim_item)
            sim_item = self.list_of_simulations.takeItem(item_row)
            
            #Adds the color back to the list of simulation colors
            self.simulation_colors.append(line_color)
            
            del sim_item

            self.update_plotting() 
            self.update_datapoints_wgt()





    def update_plotting(self): 
            
            self.plot_wdgt.clear_canvas()


            self.read_and_sort_t_tau()
            self.read_indices_for_used_temps()

            self.plot_t_tau_on_axes()
            self.redraw_simulation_lines()

            self.parent.tau_table.update_table()

            self.add_T_axis()         
            self.update_T_axis()
            self.plot_wdgt.ax.set_xlabel('1/T ($K^{-1}$)')
            self.plot_wdgt.ax.set_ylabel(r'$\ln{\tau}/s$')
        
            self.plot_wdgt.canvas.draw()
