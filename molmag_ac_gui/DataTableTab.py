#std packages
import os

#third-party packages
from PyQt5.QtWidgets import (QHBoxLayout, QLabel, QPushButton, 
                             QVBoxLayout, QWidget, QTableWidget,
                             QTableWidgetItem, QMessageBox)

from molmag_ac_gui.layout import make_btn, make_line 
from .dialogs import MagMessage

class DataTableTab(QWidget): #First tab er Qwidget, så nu er det self.layout fx i stedet for self.widget.layout

    def __init__(self, parent): 

        super(DataTableTab,self).__init__() #Før __init__() peger på QWidget
        self.parent = parent
        self.initUI() 
    
    def initUI(self):

        self.path_to_export = ""
        self.export_filetype = ""
    
        self.layout = QVBoxLayout() 

        self.line_edit = make_line(self,"When a file is loaded in the \"Data Treatment\" tab, a table of the loaded data will be shown here", self.layout )
        make_line(self,"To load a new file, go to the \"Data Tretment\" tab and click \"(1) Load datafile\" ", self.layout )
        self.tableWidget = QTableWidget() 
        self.layout.addWidget(self.tableWidget)

        make_btn(self, "Update Table of data", self.update_table_from_tab, self.layout)
        self.layoutH = QHBoxLayout()  
        make_btn(self, "Export table to .xlsx file", self.export_table_excel, self.layoutH)
        make_btn(self, "Export table to .csv file", self.export_table_csv, self.layoutH)
        make_btn(self, "Export table to .dat file", self.export_table_dat, self.layoutH)
        self.layout.addLayout(self.layoutH)
        
        self.setLayout(self.layout) 

        self.show()
    
    def export_table_dat(self): 

        self.export_filetype = "dat"
        self.export_table() 

    def update_table_from_tab(self): 

        self.update_table(ask_question=False)

    def update_table(self, ask_question = True):
        """ Updates the table """

        try: 
            df = self.parent.data_treat.raw_df
        except: 
            return
        if ask_question: 
            max_cols = 8000
            if df.shape[0] > max_cols: 
                qm = QMessageBox() 
                ans = qm.question(self,'Trying to load large dataset into Data table tab', "The dataframe is very large (more than {} rows). Updating the Table of data, might take some time. Do you want to update the Table of data anyway? \
                    The table can be updated manually from the Table of data tab".format(max_cols), qm.Yes | qm.No, defaultButton=qm.No)

                if ans == qm.No: 
                    return 
        
        self.tableWidget.setRowCount(len(df))
        self.tableWidget.setColumnCount(len(df.columns))
            
        for row in range(len(df)):
            for col in range(len(df.columns)): 
                self.tableWidget.setItem(row,col,QTableWidgetItem("{}".format(df.iloc[row,col])))

        #Makes a list of column names from df: 
        Colnames = []
        for i in range(len(df.columns)): 
            Colnames.append(df.columns.values[i])

        #Labels all collumns 
        self.tableWidget.setHorizontalHeaderLabels(Colnames)

        #Makes table read-only
        self.tableWidget.setEditTriggers(self.tableWidget.NoEditTriggers)

        #Resizes column width with respect to contests
        self.tableWidget.resizeColumnsToContents()
        
        self.line_edit.setText("The data is loaded from the file located at: {}".format(self.parent.current_file))


    def findnewpath(self): 
        #Finds path for datatable csv or excel file to be stored

        splitpath = self.parent.current_file.split('/')
        newpath = ""
        for element in splitpath[:-1]: 
            newpath += element + '/'
        newpath += splitpath[-1].split(".")[0] + "_datatable"    
        self.path_to_export = newpath
    
    def export_table_excel(self):
        """ Sets export_filetype to xlsx and exports the table """

        self.export_filetype = "xlsx"
        self.export_table() 
        
    def export_table_csv(self):
        """ Sets export_filetype to csv and exports the table """

        self.export_filetype = "csv"
        self.export_table() 

    def savetofile(self): 
        """Saves the fit to file depending on export_filetype """

        if self.export_filetype == "csv": 
            self.parent.data_treat.raw_df.to_csv(r'{}.{}'.format(self.path_to_export, self.export_filetype), index = False)
        if self.export_filetype == "xlsx": 
            self.parent.data_treat.raw_df.to_excel(r'{}.{}'.format(self.path_to_export, self.export_filetype), index = False)
        if self.export_filetype == "dat": 
            self.parent.data_treat.raw_df.to_csv(r'{}.{}'.format(self.path_to_export, self.export_filetype), index = False, header = self.parent.data_treat.raw_df_header)
            with open(r'{}.{}'.format(self.path_to_export, self.export_filetype), "r") as f: 
                csv_lines = f.readlines()  
            with open(r'{}.{}'.format(self.path_to_export, self.export_filetype), "w") as f: 
                for line in self.parent.data_treat.pre_header: 
                    f.write(line)
                for line in csv_lines: 
                    f.write(line)
        MagMessage('The data is successfully exported', "The data is saved as a .{} file at: {}.{}".format(self.export_filetype, self.path_to_export, self.export_filetype)).exec_() 


    def export_table(self):
        """ Saves file with table, checking whether the file already exists and asks the user 
        whether or not the file should be overwritten if present. """

        filetype = self.export_filetype
        try: 
            self.findnewpath()
            if os.path.isfile(r'{}.{}'.format(self.path_to_export,filetype)): #If the file already exist
                qm = QMessageBox() 
                ans = qm.question(self,'', "File already exists. Do you want to overwrite existing file?", qm.Yes | qm.No)
                if ans == qm.No: 
                    MagMessage("File not saved", "Your file has not been saved").exec_()
                if ans == qm.Yes: 
                    self.savetofile() 
            else: #If file does not already exist
                self.savetofile() 
                       
        except AttributeError: #If raw_df is empty, it will be NoneType. This type has no attribute "to_csv" or "to_excel"
            MagMessage("Error", "Cannot export file. No data is loaded in the data treatment tab.").exec_() 