from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QLabel, QHBoxLayout, QComboBox, QStackedWidget, QCheckBox, QLineEdit, QListWidget, QSplitter
from .dialogs import (PlottingWindow)

from PyQt5.QtGui import  QDoubleValidator

class DataTreatmentTab(QSplitter): 
    def __init__(self, parent): 
        super(DataTreatmentTab,self).__init__() #Før __init__() peger på QWidget
        self.parent = parent
        self.initUI() 
    
    def initUI(self): 
        self.layout = QVBoxLayout()
        #self = QSplitter()
        self.data_layout = QVBoxLayout()
        ### Making the left column (data loading and visualization controls)
        self.data_loading_wdgt = QWidget()
        self.data_layout = QVBoxLayout()
        
        self.raw_data_load_btn = QPushButton('(1) Load')
        #self.raw_data_load_btn.clicked.connect(self.load_ppms_data)
        self.data_layout.addWidget(self.raw_data_load_btn)
        
        self.diamag_correction_btn = QPushButton("(2) Diamagnetic correction")
        #self.diamag_correction_btn.clicked.connect(self.make_diamag_correction_calculation)
        self.data_layout.addWidget(self.diamag_correction_btn)
        #
        self.fit_Xp_Xpp_btn = QPushButton("(3) Fit X', X''")
        #self.fit_Xp_Xpp_btn.clicked.connect(self.fit_Xp_Xpp_standalone)
        self.data_layout.addWidget(self.fit_Xp_Xpp_btn)
        #
        self.copy_fit_to_ana_btn = QPushButton('(4) Copy for analysis')
        #self.copy_fit_to_ana_btn.clicked.connect(self.copy_fit_to_analysis)
        self.data_layout.addWidget(self.copy_fit_to_ana_btn)
        #
        self.save_fit_to_file_btn = QPushButton('(5) Save fit to file')
        #self.save_fit_to_file_btn.clicked.connect(self.save_fit_to_file)
        self.data_layout.addWidget(self.save_fit_to_file_btn)
        
        ### Constructing data plotting layout
        self.raw_data_plot_lo = QVBoxLayout()
        
        self.plot_type_header = QLabel('Axis content')
        self.plot_type_header.setFont(self.parent.headline_font)
        self.raw_data_plot_lo.addWidget(self.plot_type_header)
        
        self.plot_type_combo = QComboBox()
        self.plot_type_combo.addItems(['Raw data', 'Fitted', 'Temp VS Freq VS Xpp'])
        #self.plot_type_combo.currentIndexChanged.connect(self.switch_analysis_view)
        self.raw_data_plot_lo.addWidget(self.plot_type_combo)
        #
        self.raw_data_plot_header = QLabel('Raw data plotting')
        self.raw_data_plot_header.setFont(self.parent.headline_font)
        self.raw_data_plot_lo.addWidget(self.raw_data_plot_header)
        #
        ## Constructing the x combobox
        self.x_lo = QHBoxLayout()
        self.x_combo_lbl = QLabel('x')
        self.x_lo.addWidget(self.x_combo_lbl)
        #
        self.x_combo = QComboBox()
        #self.x_combo.currentIndexChanged.connect(self.plot_from_combo)
        self.x_lo.addWidget(self.x_combo)
        self.raw_data_plot_lo.addLayout(self.x_lo)
        #
        ## Constructing the y combobox
        self.y_lo = QHBoxLayout()
        self.y_combo_lbl = QLabel('y')
        self.y_lo.addWidget(self.y_combo_lbl)
        
        self.y_combo = QComboBox()
        #self.y_combo.currentIndexChanged.connect(self.plot_from_combo)
        self.y_lo.addWidget(self.y_combo)
        self.raw_data_plot_lo.addLayout(self.y_lo)
        #
        ## Constructing a combobox for plotting fitted data
        self.fitted_data_plot_header = QLabel('Fit data plotting')
        self.fitted_data_plot_header.setFont(self.parent.headline_font)
        self.raw_data_plot_lo.addWidget(self.fitted_data_plot_header)
        #
        self.fit_data_plot_combo = QComboBox()
        self.fit_data_plot_combo.addItems(['ColeCole', 'FreqVSXp', 'FreqVSXpp'])
        #self.fit_data_plot_combo.currentIndexChanged.connect(self.plot_from_itemlist)
        self.raw_data_plot_lo.addWidget(self.fit_data_plot_combo)
#
#
        #
        ## Checkbox for coloring measured points
        self.fit_data_color_cb_lo = QHBoxLayout()
        self.fit_data_color_cb_header = QLabel('Use temp. colors')
        self.fit_data_color_cb_lo.addWidget(self.fit_data_color_cb_header)
        #
        self.fit_data_color_cb = QCheckBox()
        #self.fit_data_color_cb.stateChanged.connect(self.plot_from_itemlist)
        self.fit_data_color_cb_lo.addWidget(self.fit_data_color_cb)
        #
        #self.raw_data_plot_lo.addLayout(self.fit_data_color_cb_lo)
    #
        ## Finalizing the raw data layout
        self.data_layout.addLayout(self.raw_data_plot_lo)
        #
        ### Finalizing the data loading widget
        self.data_layout.addStretch()
        self.data_loading_wdgt.setLayout(self.data_layout)
        #
        ## Making the middle of the tab (data visualization)
        self.sw = QStackedWidget()
        #
        self.raw_plot = PlottingWindow()
        self.fit_plot = PlottingWindow(make_ax="cax")
        self.threeDplot = PlottingWindow(make_ax = "z") 
#
        self.sw.addWidget(self.raw_plot)
        self.sw.addWidget(self.fit_plot)
        self.sw.addWidget(self.threeDplot)
#
        ### Making the right column (parameter controls)
        self.param_wdgt = QWidget()
        self.param_layout = QVBoxLayout()
        #
        ## SAMPLE INFO
        self.sample_info_layout = QVBoxLayout()
        
        self.sample_info_header = QLabel('Sample information')
        self.sample_info_header.setFont(self.parent.headline_font)
        self.sample_info_layout.addWidget(self.sample_info_header)
        #
        ## Sample mass edit
        self.sample_mass_layout = QHBoxLayout()
        self.sample_info_layout.addLayout(self.sample_mass_layout)
        #
        self.sample_mass_lbl = QLabel('m (sample) [mg]')
        self.sample_mass_layout.addWidget(self.sample_mass_lbl)
        #
        self.sample_mass_inp = QLineEdit()
        self.sample_mass_inp.setValidator(QDoubleValidator())
        self.sample_mass_layout.addWidget(self.sample_mass_inp)
#
        ## Sample molar mass edit
        self.molar_mass_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.molar_mass_lo)
        #
        self.molar_mass_lbl = QLabel('M (sample) [g/mol]')
        self.molar_mass_lo.addWidget(self.molar_mass_lbl)
        #
        self.molar_mass_inp = QLineEdit()
        self.molar_mass_inp.setValidator(QDoubleValidator())
        self.molar_mass_lo.addWidget(self.molar_mass_inp)
#
        ## Sample Xd edit
        self.sample_xd_lo = QHBoxLayout()
        self.sample_info_layout.addLayout(self.sample_xd_lo)
        

        self.sample_xd_lbl = QLabel(u"<a href={}>X\u1D05</a>".format(
                                    self.parent.tooltips_dict['English']['Xd_link'])
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
        #self.load_sample_data_btn.clicked.connect(self.load_sample_data)
        self.sample_data_lo.addWidget(self.load_sample_data_btn)
        
        self.save_sample_data_btn = QPushButton('Save sample data')
        #self.save_sample_data_btn.clicked.connect(self.save_sample_data)
        self.sample_data_lo.addWidget(self.save_sample_data_btn)
        
        self.sample_data_lo.addStretch()
        
        self.param_layout.addLayout(self.sample_info_layout)
        
        # List of fitted raw data
        self.raw_fit_headline = QLabel('Fitted parameters')
        self.raw_fit_headline.setFont(self.parent.headline_font)
        self.param_layout.addWidget(self.raw_fit_headline)
        
        self.raw_fit_list = QListWidget()
        self.param_layout.addWidget(self.raw_fit_list)
        #self.raw_fit_list.doubleClicked.connect(self.update_raw_plot)
        
        ## Finalizing layout
        
        self.param_layout.addStretch()
        self.param_wdgt.setLayout(self.param_layout)
        
        self.addWidget(self.data_loading_wdgt)
        self.addWidget(self.sw)
        self.addWidget(self.param_wdgt)


        self.show()