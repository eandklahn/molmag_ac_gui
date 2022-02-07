#std packages 
import os
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
                             QLabel, QHBoxLayout, QCheckBox, 
                             QListWidget, QSplitter)
from scipy.optimize.minpack import curve_fit
import scipy.constants as sc

#local imports
from .exceptions import NoGuessExistsError
from .process_ac import (getParameterGuesses, tau_err_RC, fit_relaxation,
                         default_parameters, add_partial_model)
from .dialogs import (GuessDialog, MagMessage, ParamDialog, 
                      PlottingWindow, SimulationDialog)
from .layout import make_headline, make_btn, make_line

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
        make_btn(self, "Import current fit from Data Treatment", self.parent.data_treat.copy_fit_to_analysis, self.options_layout)
        make_btn(self, 'Load file generated in Data Treatment', self.load_t_tau_data, self.options_layout)

        # Adding fit controls with checkboxes
        make_headline(self, "Fitting options", self.options_layout)
        make_line(self, "Choose fit types to include: ", self.options_layout)
        self.add_fit_type_checkboxes()

        # Adding temperature controls
        make_line(self, "Choose temperature range to fit in: ", self.options_layout)
        self.add_temp_controls() 

        # Adding a button to run a fit
        make_btn(self, "Run fit!", self.make_the_fit, self.options_layout)
        

        # Adding a list to hold information about simulations
        make_headline(self, "Simulations", self.options_layout)
        self.add_simulations_list()

        # Adding buttons to control simulation list
        self.sim_btn_layout = QHBoxLayout()
        make_btn(self, "New", self.edit_simulation_from_list, self.sim_btn_layout)
        make_btn(self, "Delete", self.delete_sim, self.sim_btn_layout)
        make_btn(self, "Edit", self.edit_simulation_from_list, self.sim_btn_layout)
        self.options_layout.addLayout(self.sim_btn_layout)
        
        #Adding view fitted parameters button 
        make_headline(self, "View fitted parameters", self.options_layout)        
        make_btn(self, "Fitted params", self.show_fitted_params, self.options_layout)

        #Setting the layout of the options widget
        self.options_wdgt.setLayout(self.options_layout)
        self.options_layout.addStretch() 
        self.addWidget(self.options_wdgt)
        
        #Adding plotting widget
        self.add_plot_wdgt() 

        # Finalizing layout of the data analysis tab
        self.setSizes([1,1200])
        self.show() 


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
        
        self.temp_line[1].setRange(0,self.temp_line[3].value())
        self.temp_line[1].setSingleStep(0.1)
        self.temp_line[3].setRange(self.temp_line[1].value(),1000)
        self.temp_line[3].setSingleStep(0.1)
        
        self.temp_line[1].editingFinished.connect(self.set_new_temp_ranges)
        self.temp_line[3].editingFinished.connect(self.set_new_temp_ranges)
        for w in self.temp_line:
            self.temp_horizontal_layout.addWidget(w)

        self.temp_horizontal_layout.setAlignment(Qt.AlignCenter)
        self.options_layout.addLayout(self.temp_horizontal_layout)

    def add_simulations_list(self): 
        """Makes a QListWidget with all simulations and adds it to the options layout"""

        self.list_of_simulations = QListWidget()
        self.list_of_simulations.setDragDropMode(self.list_of_simulations.InternalMove)
        
        self.list_of_simulations.doubleClicked.connect(self.edit_simulation_from_list)
        """https://stackoverflow.com/questions/41353653/how-do-i-get-the-checked-items-in-a-qlistview"""
        self.list_of_simulations.itemChanged.connect(self.redraw_simulation_lines)
        
        self.options_layout.addWidget(self.list_of_simulations)


    def add_fit_type_checkboxes(self):
        """Creates a QHBoxLayout with three checkboxes and adds these to the options layout"""

        self.fit_type_layout = QHBoxLayout() 
        
        self.orbach_cb = QCheckBox('Orbach')
        self.orbach_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_type_layout.addWidget(self.orbach_cb)
        
        self.raman_cb = QCheckBox('Raman')
        self.raman_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_type_layout.addWidget(self.raman_cb)
        
        self.qt_cb = QCheckBox('QT')
        self.qt_cb.stateChanged.connect(self.read_fit_type_cbs)
        self.fit_type_layout.addWidget(self.qt_cb)
        
        self.options_layout.addLayout(self.fit_type_layout)


    def add_plot_wdgt(self): 
        """Adds a plotting widget for 1/T vs ln(tau) plot """

        self.plot_wdgt = PlottingWindow()
        self.plot_wdgt.ax.set_xlabel('1/T ($K^{-1}$)')
        self.plot_wdgt.ax.set_ylabel(r'$\ln{\tau}$ ($\ln{s}$)')
        self.addWidget(self.plot_wdgt)


    def show_fitted_params(self):
        
        try:
            fit = self.fit_history[0]
        except (IndexError, TypeError):
            w = MagMessage('Fit history error', 'There is no fit history yet!')
        else:
            w = ParamDialog(self, self.fit_history)
        finally:
            w.exec_()

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

    def plot_t_tau_on_axes(self):
        
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
        
        self.plot_wdgt.canvas.draw()

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
            self.set_new_t_tau(D)
            self.read_indices_for_used_temps()
            self.plot_t_tau_on_axes()
            self.plot_wdgt.reset_axes()


    def read_fit_type_cbs(self):
    
        list_of_checked = []
        if self.qt_cb.isChecked(): list_of_checked.append('QT')
        if self.raman_cb.isChecked(): list_of_checked.append('R')
        if self.orbach_cb.isChecked(): list_of_checked.append('O')
        fitToMake = ''.join(list_of_checked)
        
        return fitToMake

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
            msg.setIcon(QMessageBox.Information)
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

    def redraw_simulation_lines(self):
        
        for idx in range(self.list_of_simulations.count()):
            item = self.list_of_simulations.item(idx)
            data = item.data(32)
            
            if item.checkState() == Qt.Checked:
                data['line']._visible = True
            elif item.checkState() == Qt.Unchecked:
                data['line']._visible = False
        
        self.plot_wdgt.canvas.draw()

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

    def delete_sim(self):
        
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
        
        fs = [p for p in params if 'use' in p]
        used = [bool(params[p].value) for p in fs]
        text = ['Using QT: {}, Raman: {}, Orbach: {}\n'.format(*used),
                'to plot between {} K and {} K\n'.format(*T_vals)]
        return ''.join(text)