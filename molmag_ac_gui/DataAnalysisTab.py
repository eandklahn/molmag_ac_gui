#std packages 
import os
from re import X
import sys
from collections import deque
import datetime

#third-party packages
import numpy as np
import names
from matplotlib.colors import to_hex
from matplotlib._color_data import TABLEAU_COLORS
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QDoubleSpinBox, QFileDialog, QListWidgetItem, 
                             QMessageBox, QWidget, QVBoxLayout, 
                             QLabel, QHBoxLayout, 
                             QListWidget, QSplitter, QComboBox, QTextEdit)
from scipy.optimize.minpack import curve_fit
import scipy.constants as sc

#local imports
from .exceptions import NoGuessExistsError
from .process_ac import (getParameterGuesses, tau_err_RC, fit_relaxation,
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
        self.startUp = True #Some old stuff about reading a file directly from the terminal. Does not work fully. 

        # Initializes simulations colors and other attributes
        self.initialize_attributes() 

        # The options widget to the left 
        self.options_wdgt = QWidget()
        self.options_layout = QVBoxLayout()
        
        # Adding data loading options
        make_headline(self, "Data loading options", self.options_layout)
        make_btn(self, "Import current fit from Data Treatment", self.copy_data_treat_fit, self.options_layout)
        make_btn(self, 'Load fitting file generated in Data Treatment', self.load_t_tau_data, self.options_layout)

        # Adding fit controls with checkboxes
        make_headline(self, "Fitting options", self.options_layout)
        make_line(self, "Choose fit types to include: ", self.options_layout)
        self.add_fit_type_checkboxes()

        # Adding temperature controls
        make_line(self, "Choose temperature range to fit in: ", self.options_layout)
        self.add_temp_controls() 

        # Adding a button to run a fit
        make_btn(self, "Run fit!", self.make_the_fit, self.options_layout)

         #Adding view fitted parameters button 
        make_headline(self, "View information about fits", self.options_layout)        
        self.add_fit_parameters_view()       

        # Adding a list to hold information about simulations
        make_headline(self, "Plot/simulate fit", self.options_layout)
        self.add_simulations_list()

        # Adding buttons to control simulation list
        self.sim_btn_layout = QHBoxLayout()
        make_btn(self, "New", self.edit_simulation_from_list, self.sim_btn_layout)
        make_btn(self, "Delete", self.delete_sim, self.sim_btn_layout)
        self.options_layout.addLayout(self.sim_btn_layout)


        #Setting the layout of the options widget
        self.options_wdgt.setLayout(self.options_layout)
        self.options_layout.addStretch() 
        self.addWidget(self.options_wdgt)
        
        #Adding plotting widget
        self.add_plot_wdgt() 

        # Finalizing layout of the data analysis tab
        self.setSizes([1,1200])
        self.show() 

    def copy_data_treat_fit(self): 
        self.fit_history.clear() 
        self.update_fit_combo() 
        self.reset_checkboxes_and_T() 
        self.fit_stat_txt = ""
        self.fit_summary.setText("")
        self.fit_title.setText("")
        self.parent.data_treat.copy_fit_to_analysis() 

    def reset_checkboxes_and_T(self): 
        self.orbach_cb.setChecked(False)
        self.raman_cb.setChecked(False)
        self.qt_cb.setChecked(False)
        self.temp_line[1].setValue(0)
        self.temp_line[3].setValue(0)

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
        name, ext = os.path.splitext(filename)
        if ext == '':
            ext = '.txt'
        
        try: 
            self.set_fit_stat_txt(all_significant_digits=True) 
            with open(filename + ext, "w") as f:
                f.write(self.fit_stat_txt)
            MagMessage("Succes", "File succesfully written").exec_() 
       
        except IndexError: #If fit history is empty 
            MagMessage("Error", "No fit made yet").exec_() 
       
        except ValueError: #If fit_result format is different, it will save in original format
            fit_idx = self.choose_fit_combo.currentIndex() 
            name, res, time = self.fit_history[fit_idx]
            self.fit_stat_txt = fit_report(res)
            with open(filename + ext, "w") as f:
                f.write(self.fit_stat_txt)
            MagMessage("Succes", "File succesfully written").exec_() 

    def set_fit_stat_txt(self, all_significant_digits = False): 
        """ Formatting the result from the fit_report and saving the formatted 
        text in self.fit_stat_txt. Displays them with only few siginificant digits """
        
        self.fit_stat_txt = "" 

        fit_idx = self.choose_fit_combo.currentIndex() 
        name, res, time = self.fit_history[fit_idx]
        joined_result = "".join(fit_report(res))
        lines = [line.split() for line in joined_result.split("\n")]
        
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
            name, res, time = fit
            self.choose_fit_combo.addItem(f'{time}: {name}')

    def show_fit_stat(self):
        """ Shows and sets the fit statistics text in the parameter view"""

        if self.fit_history == []: 
            return 
        
        fit_idx = self.choose_fit_combo.currentIndex()
        name, res, time = self.fit_history[fit_idx]
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
        
        self.fitted_params_dialog = None

    def add_temp_controls(self): 
        """Makes a QHBoxLayout with a two spinboxes where the temperature range is chosen 
        and adds it to the options layout."""

        self.temp_horizontal_layout = QHBoxLayout()
        self.temp_line = [QLabel('Temperature range in K: ('), QDoubleSpinBox(), QLabel(','),
                          QDoubleSpinBox(), QLabel(')')]
        
        self.temp_line[1].setRange(0,1000)
        self.temp_line[1].setSingleStep(0.1)
        self.temp_line[3].setRange(0,1000)
        self.temp_line[3].setSingleStep(0.1)

        self.temp_line[1].editingFinished.connect(self.set_new_temp_ranges)
        self.temp_line[3].editingFinished.connect(self.set_new_temp_ranges)
        for w in self.temp_line:
            self.temp_horizontal_layout.addWidget(w)

        self.temp_horizontal_layout.setAlignment(Qt.AlignCenter)
        self.options_layout.addLayout(self.temp_horizontal_layout)



    
    def add_simulations_list(self): 
        """ Makes a QListWidget with all simulations and adds it to the options layout """

        self.list_of_simulations = QListWidget()
        self.list_of_simulations.setDragDropMode(self.list_of_simulations.InternalMove)
        
        self.list_of_simulations.doubleClicked.connect(self.edit_simulation_from_list)
        """https://stackoverflow.com/questions/41353653/how-do-i-get-the-checked-items-in-a-qlistview"""
        self.list_of_simulations.itemChanged.connect(self.redraw_simulation_lines)
        
        self.options_layout.addWidget(self.list_of_simulations)


    def add_fit_type_checkboxes(self):
        """Creates a QHBoxLayout with three checkboxes and adds these to the options layout"""

        self.fit_type_layout = QHBoxLayout() 
        
        self.orbach_cb = make_checkbox(self, self.read_fit_type_cbs, self.fit_type_layout, "Orbach") 
        self.raman_cb = make_checkbox(self, self.read_fit_type_cbs, self.fit_type_layout, "Raman") 
        self.qt_cb = make_checkbox(self, self.read_fit_type_cbs, self.fit_type_layout, "QT")
        
        self.options_layout.addLayout(self.fit_type_layout)


    def add_plot_wdgt(self): 
        """ Adds a plotting widget for 1/T vs ln(tau) plot """

        self.plot_wdgt = PlottingWindow()
        self.plot_wdgt.ax.set_xlabel('1/T ($K^{-1}$)')
        self.plot_wdgt.ax.set_ylabel(r'$\ln{\tau}$ ($\ln{s}$)')
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
        
        self.used_indices = None


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


        if D.shape[1]==3 or D.shape[1] == 8: 
            # Three columns in the array loaded or 8 columns loaded from file
            # Assume the third columns is error
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

    def read_indices_for_used_temps(self):
        """Read the values for the chosen temperature range, and uses their indicies to 
        update self.used_T, self.used_tau, self.not_used_T and self.not_used_tau  """
        
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
            MagMessage("Error", 'No data have been selected yet!').exec_() 

    def plot_t_tau_on_axes(self):
        """ Plots 1/T vs. τ in the plotting window. Also adds errorbars if errors on τ are present """

        if self.plotted_data_pointers is not None:
            for line in self.plotted_data_pointers:
                line.remove()
        self.plotted_data_pointers = []
        
        if self.data_dtau is None:
            used, = self.plot_wdgt.ax.plot(1/self.used_T, np.log(self.used_tau), 'bo', zorder=0.1)
            not_used, = self.plot_wdgt.ax.plot(1/self.not_used_T, np.log(self.not_used_tau), 'ro', zorder=0.1)
            self.plotted_data_pointers.append(used)
            self.plotted_data_pointers.append(not_used)
        else:

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
                                                                                    fmt='ro',
                                                                                    ecolor='r',
                                                                                    label='Data',
                                                                                    zorder=0.1)

                

            self.plotted_data_pointers.append(err_used_point)
            self.plotted_data_pointers.append(err_not_used_point)
            for e in [caplines1, caplines2, barlinecols1, barlinecols2]:
                for line in e:
                    self.plotted_data_pointers.append(line)
        
            self.plot_wdgt.ax.relim() 
            self.plot_wdgt.ax.autoscale_view() 
        self.plot_wdgt.canvas.draw()

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
            V = 1/X 
            return ["%.2f" % z for z in V]

        self.plot_wdgt.ax2.set_xlim(self.plot_wdgt.ax.get_xlim())
        self.plot_wdgt.ax2.set_xticks(T_tick_loc)
        self.plot_wdgt.ax2.set_xticklabels(tick_f(T_tick_loc))
        self.plot_wdgt.canvas.draw() 

    def load_t_tau_data(self):
        """Load 1/T vs. τ data from file generated in data treatment """

        if self.startUp:
            try:
                filename = sys.argv[1]
            except IndexError:
                pass
            finally:
                self.startUp = False
                return 0
        else:
            filename_info = QFileDialog().getOpenFileName(self, 'Open file', self.parent.last_loaded_file)
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
                msg = MagMessage("ValueError", 'File format not as expected')
                msg.setIcon(QMessageBox.Warning)
                msg.exec_()
            elif error_type == 'OSError':
                pass
        else:
            #self.plot_wdgt.ax.clear() #Clearing data in plotting wdgt 
            #self.fit_history.clear() 

            self.set_new_t_tau(D)
            self.read_indices_for_used_temps()
            self.plot_t_tau_on_axes()
            self.add_T_axis()         
            self.plot_wdgt.reset_axes()
            self.update_T_axis()
            
            #self.update_fit_combo()



    def read_fit_type_cbs(self):
        """ Read the fit type (i.e. which of Quantum Tunneling, Raman and Orbach are 
        included in the fitting). Depends on what is chosen in the comboboxes """
    
        list_of_checked = []
        if self.qt_cb.isChecked(): list_of_checked.append('QT')
        if self.raman_cb.isChecked(): list_of_checked.append('R')
        if self.orbach_cb.isChecked(): list_of_checked.append('O')
        fitToMake = ''.join(list_of_checked)
        
        return fitToMake

    def set_new_temp_ranges(self):
        """ Sets a new temperature range to fit in, when the spin boxes are changed"""

        self.read_indices_for_used_temps()
        if self.data_T is not None:
            self.plot_t_tau_on_axes()
        

    def make_the_fit(self):
        """ Makes the fit of 1/T vs. τ"""

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
            assert Tmin < Tmax 
            
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
            msg_text = "Bad temperature or fit settings \nPossible errors: \n \
            - Minimum and maximum temperatures are the same \n \
            - Minimum temperature is larger than maximum temperature \n \
            - No fit options have been selected \n \
            - Can't fit only one data point "
        except RuntimeError:
            msg_text = 'This fit cannot be made within the set temperatures'
        #except ValueError as e:
        #    msg_text = 'No file has been loaded'
        except TypeError as e:
            msg_text = 'No data has been selected'
            print(e)
        except NoGuessExistsError:
            msg_text = 'Made no guess for initial parameters'
        
        else:
            window_title = 'Fit successful!'
            msg_text = 'Congratulations'
            
            self.add_to_history(minimize_res, fitwith)
            self.update_fit_combo() 
        finally:
            msg = MagMessage(window_title, msg_text)
            msg.setIcon(QMessageBox.Information)
            msg.exec_()  
      
    
    def add_to_history(self, p_fit, perform_this_fit):
        """ Adds a fit to the top of the fit history"""

        now = datetime.datetime.now()
        now = now.strftime("%d.%m %H:%M:%S")

        if len(self.fit_history)>9:
            self.fit_history.pop()
        self.fit_history.insert(0, (perform_this_fit, p_fit, now))



    def redraw_simulation_lines(self):
        """ Redraw simulation lines """

        for idx in range(self.list_of_simulations.count()):
            item = self.list_of_simulations.item(idx)
            data = item.data(32)
            
            if item.checkState() == Qt.Checked:
                data['line']._visible = True
            elif item.checkState() == Qt.Unchecked:
                data['line']._visible = False
        
        self.plot_wdgt.canvas.draw()

    def edit_simulation_from_list(self):
        """ Edit simulations/fits from the list"""

        action = self.get_action() 

        if action == 'Edit':
            params, T_vals, line, color, label  = self.load_chosen_item() 

        elif action == 'New':
            params, T_vals, line, color, label = self.load_new_item() 

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
                self.plot_wdgt.ax.lines.remove(line)
            else:
                # In this case, there was no old line and therefore also no sim_item
                """https://stackoverflow.com/questions/55145390/pyqt5-qlistwidget-with-checkboxes-and-drag-and-drop"""
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
            
            new_item_text = self.represent_simulation(T_vals, params)
            
            sim_item.setData(32, list_item_data)
            sim_item.setText(new_item_text)
            sim_item.setBackground(QColor(to_hex(color)))
            
            self.redraw_simulation_lines()

    def get_action(self): 
        """ Get the action (edit or new) depending on the sender chosen. 
        Used for edit_simulations_from_list"""
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
        return action

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

            return params, T_vals, line, color, label 


    def load_new_item(self): 
        """ If new is chosen for the simulation list, the following values are set """
        
        params = default_parameters()
        T_vals = [1,3]
        line = False
        label = None
        color = None
        if len(self.simulation_colors)<1:
            self.statusBar.showMessage("ERROR: can't make any more simulations")
            return

        return params, T_vals, line, color, label 
    
    def delete_sim(self):
        """ Deleting a simulation/fit from the list"""

        try:
            sim_item = self.list_of_simulations.selectedItems()[0]
        except IndexError:
            pass
        else:
            line_pointer = sim_item.data(32)['line']
            line_color = line_pointer._color
            
            self.plot_wdgt.ax.lines.remove(line_pointer)
            self.plot_wdgt.canvas.draw()
            
            item_row = self.list_of_simulations.row(sim_item)
            sim_item = self.list_of_simulations.takeItem(item_row)
            
            self.simulation_colors.append(line_color)
            
            del sim_item

    def represent_simulation(self, T_vals, params):
        """ The string displayed to show information on each fit """
        
        fs = [p for p in params if 'use' in p]
        used = [bool(params[p].value) for p in fs]
        text = ['Using QT: {}, Raman: {}, Orbach: {}\n'.format(*used),
                'to plot between {} K and {} K\n'.format(*T_vals)]
        return ''.join(text)