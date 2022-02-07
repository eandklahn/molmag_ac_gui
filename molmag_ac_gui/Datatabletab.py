#std packages
import os

#third-party packages
import pandas as pd
from PyQt5.QtWidgets import (QHBoxLayout, QLabel, QPushButton, 
                             QVBoxLayout, QWidget, QTableWidget,
                             QTableWidgetItem, QMessageBox) 
from .dialogs import MagMessage

class DataTableTab(QWidget): #First tab er Qwidget, så nu er det self.layout fx i stedet for self.widget.layout

    def __init__(self, parent): 
        super(DataTableTab,self).__init__() #Før __init__() peger på QWidget
        self.parent = parent
        self.initUI() 
    
    def initUI(self): 
        self.path_to_export = ""
        self.export_filetype = ""
    
        self.layout = QVBoxLayout() #Vertical box layout. Første ting øverst, så nedenunder osv. 
        self.layoutH = QHBoxLayout() #Vertical box layout. Første ting øverst, så nedenunder osv. 

        self.line_edit = QLabel("When a file is loaded in the \"Data Treatment\" tab, a table of the loaded data will be shown here")

        self.tableWidget = QTableWidget() 
        
        self.export_table_excel_btn = QPushButton('Export table to Excel (.xlsx file)')
        self.export_table_excel_btn.clicked.connect(self.export_table_excel)

        self.export_table_csv_btn = QPushButton('Export table to .csv file')
        self.export_table_csv_btn.clicked.connect(self.export_table_csv)

        self.layout.addWidget(self.line_edit)     
        self.layout.addWidget(self.tableWidget)
        self.layoutH.addWidget(self.export_table_excel_btn)
        self.layoutH.addWidget(self.export_table_csv_btn)

        self.setLayout(self.layout) 
        self.layout.addLayout(self.layoutH)

        
        self.show()
    
    def updatetable(self):
        
        self.line_edit.setText("The data is loaded from the file located at: {}".format(self.parent.current_file))

        df = self.parent.data_treat.raw_df
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


    def findnewpath(self): 
        #Finds path for datatable csv or excel file to be stored
        splitpath = self.parent.current_file.split('/')
        newpath = ""
        for element in splitpath[:-1]: 
            newpath += element + '/'
        newpath += splitpath[-1] + "_datatable"    
        self.path_to_export = newpath
    
    def export_table_excel(self):
        self.export_filetype = "xlsx"
        self.export_table() 
        
    def export_table_csv(self):
        self.export_filetype = "csv"
        self.export_table() 

    def savetofile(self): 
        if self.export_filetype == "csv": 
            self.parent.data_treat.raw_df.to_csv(r'{}.{}'.format(self.path_to_export, self.export_filetype), index = False)
        if self.export_filetype == "xlsx": 
            self.parent.data_treat.raw_df.to_excel(r'{}.{}'.format(self.path_to_export, self.export_filetype), index = False)
        MagMessage('The data is successfully exported', "The data is saved as a .{} file at: {}.{}".format(self.export_filetype, self.path_to_export, self.export_filetype)).exec_() 

    def export_table(self):
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