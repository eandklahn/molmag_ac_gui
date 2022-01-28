#std packages 
import os
import sys
from collections import deque

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
from .process_ac import (addPartialModel, getFittingFunction, 
                         getParameterGuesses, getStartParams, readPopt,
                         tau_err_RC, _QT, _R, _O)
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
        
        #Adding view fitted parameters button 
        make_headline(self, "View fitted parameters", self.options_layout)        
        make_btn(self, "Fitted params", self.show_fitted_params, self.options_layout)

        # Adding a list to hold information about simulations
        make_headline(self, "Simulations", self.options_layout)
        self.add_simulations_list()

        # Adding buttons to control simulation list
        self.sim_btn_layout = QHBoxLayout()
        make_btn(self, "New", self.edit_simulation_from_list, self.sim_btn_layout)
        make_btn(self, "Delete", self.delete_sim, self.sim_btn_layout)
        make_btn(self, "Edit", self.edit_simulation_from_list, self.sim_btn_layout)
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



    def initialize_attributes(self):  

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
            assert self.fitted_params_dialog is None
            
        except IndexError:
            pass
        except AssertionError:
            self.fitted_params_dialog.show()
            self.fitted_params_dialog.activateWindow()
        else:
            self.fitted_params_dialog = ParamDialog(fit_history=self.fit_history)
            finished = self.fitted_params_dialog.exec_()

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

    def add_to_history(self, p_fit, perform_this_fit):
    
        if len(self.fit_history)>9:
            self.fit_history.pop()
        self.fit_history.insert(0, (perform_this_fit, p_fit))

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

    def fit_relaxation_lmfit(self, guess, perform_this_fit):

        fn = self.make_fitting_function(perform_this_fit)
        
        def objective(params, T, tau_obs, err=None):
            
            tQT = params['tQT']
            Cr = params['Cr']
            n = params['n']
            t0 = params['t0']
            Ueff = params['Ueff']
            
            tau_calc = fn(T, tQT, Cr, n, t0, Ueff)

            residuals = tau_calc - tau_obs
            
            if err is not None:
                residuals = residuals/err
                
            return residuals

    def make_the_fit(self):
        window_title = 'Fit aborted'
        msg_text = ''
        msg_details = ''
            
        try:
            
            # This will raise TypeError and IndexError first
            # to warn that no data was loaded
            guess = getParameterGuesses(self.used_T, self.used_tau)
            
            Tmin = self.temp_line[1].value()
            Tmax = self.temp_line[3].value()
            perform_this_fit = self.read_fit_type_cbs()
            
            assert Tmin != Tmax
            assert perform_this_fit != ''
            
            guess_dialog = GuessDialog(guess=guess,
                                       fit_history=self.fit_history)
            accepted = guess_dialog.exec_()
            if not accepted: raise NoGuessExistsError
            
            # If both fit and temperature setting are good,
            # and the GuessDialog was accepted, get the
            # guess and perform fitting
            guess = guess_dialog.return_guess
            p_fit = self.fit_relaxation(guess, perform_this_fit)
            
        except (AssertionError, IndexError):
            msg_text = 'Bad temperature or fit settings'
            msg_details = """Possible errors:
            - min and max temperatures are the same
            - no fit options have been selected
            - can't fit only one data point"""
        except RuntimeError:
            msg_text = 'This fit cannot be made within the set temperatures'
        except ValueError as e:
            print(e)
            msg_text = 'No file has been loaded'
        except TypeError:    
            msg_text = 'No data has been selected'
        except NoGuessExistsError:
            msg_text = 'Made no guess for initial parameters'
        
        else:
            window_title = 'Fit successful!'
            msg_text = 'Congratulations'
            
            self.add_to_history(p_fit, perform_this_fit)
            
        finally:
            msg = MagMessage(window_title,msg_text)
            msg.setIcon(QMessageBox.Warning)
            msg.setDetailedText(msg_details)
            msg.exec_()        

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
                print('Did not find any selected line')
                return
            else:

                # Reading off information from the selected item
                old_data = sim_item.data(32)
                old_plot_type_list = old_data['plot_type']
                old_p_fit = old_data['p_fit']
                old_T_vals = old_data['T_vals']
                old_line = old_data['line']
                old_color = old_line._color
                old_label = old_line._label
        
        elif action == 'New':
            
            old_plot_type_list = []
            old_p_fit = {'tQT': 0.1, 'Cr': 0.1, 'n': 0.1, 't0': 0.1, 'Ueff': 0.1}
            old_T_vals = [0,0]
            old_line = False
            old_label = None
            old_color = None
            if len(self.simulation_colors)<1:
                self.parent.statusBar.showMessage("ERROR: can't make any more simulations")
                return
            
        sim_dialog = SimulationDialog(fit_history=self.fit_history,
                                      plot_type_list=old_plot_type_list,
                                      plot_parameters=old_p_fit,
                                      min_and_max_temps=old_T_vals)
        finished_value = sim_dialog.exec_()
        
        try:
            assert finished_value
            
            new_plot_type = sim_dialog.plot_type_list
            assert len(new_plot_type)>0
            
        except AssertionError:
            pass
            
        else:
            
            new_p_fit = sim_dialog.plot_parameters
            new_T_vals = sim_dialog.min_and_max_temps
            plot_to_make = ''.join(new_plot_type)
            
            if old_line:
                self.plot_wdgt.ax.lines.remove(old_line)
            else:
                # In this case, there was no old line and therefore also no sim_item
                
                """https://stackoverflow.com/questions/55145390/pyqt5-qlistwidget-with-checkboxes-and-drag-and-drop"""
                sim_item = QListWidgetItem()
                sim_item.setFlags( sim_item.flags() | Qt.ItemIsUserCheckable )
                sim_item.setCheckState(Qt.Checked)
                
                self.list_of_simulations.addItem(sim_item)
                old_color = self.simulation_colors.pop()
                old_label = names.get_first_name()
            
            new_line = addPartialModel(self.plot_wdgt.fig,
                                       new_T_vals[0],
                                       new_T_vals[1],
                                       self.prepare_sim_dict_for_plotting(new_p_fit),
                                       plotType=plot_to_make,
                                       c=old_color,
                                       label=old_label)
            
            list_item_data = {'plot_type': new_plot_type,
                              'p_fit': new_p_fit,
                              'T_vals': new_T_vals,
                              'line': new_line,
                              'color': old_color}
            
            new_item_text = f"{old_label},\n({new_T_vals[0]:.1f},{new_T_vals[1]:.1f}),\n"
            new_item_text += f"tQT: {new_p_fit['tQT']:.2e}, Cr: {new_p_fit['Cr']:.2e}, "
            new_item_text += f"n: {new_p_fit['n']:.2e}, t0: {new_p_fit['t0']:.2e}, "
            new_item_text += f"Ueff: {new_p_fit['Ueff']:.2e}"
            
            sim_item.setData(32, list_item_data)
            sim_item.setText(new_item_text)
            sim_item.setBackground(QColor(to_hex(old_color)))
            
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

