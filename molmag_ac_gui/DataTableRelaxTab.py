#std packages
import os

#third-party packages
from PyQt5.QtWidgets import (QHBoxLayout, QLabel, QPushButton, 
                             QVBoxLayout, QWidget, QTableWidget,
                             QTableWidgetItem, QMessageBox, QSizePolicy)

from molmag_ac_gui.layout import make_btn, make_line 
from .dialogs import MagMessage

class DataTableRelaxTab(QWidget): #First tab er Qwidget, så nu er det self.layout fx i stedet for self.widget.layout

    def __init__(self, parent): 

        super(DataTableRelaxTab,self).__init__() #Før __init__() peger på QWidget
        self.parent = parent
        self.initUI() 
    
    def initUI(self):

        self.path_to_export = ""
        self.export_filetype = ""
    
        self.layout = QVBoxLayout() 

        make_line(self, "When data has been imported into the Data Analysis tab, it will be shown here.", self.layout)
        make_line(self, "If you delete datapoints (outliers) in the Data Analysis tab, the datapoints will also be deleted from the table here. You can then save the data in a new .dat file by using the button below", self.layout)

        self.tableWidget = QTableWidget() 
        self.layout.addWidget(self.tableWidget)
        
        make_btn(self,  "Save fit to file (AC)", self.parent.data_treat.save_fit_to_file, self.layout)

        self.setLayout(self.layout) 
        self.show()
    
    def export_table_dat(self): 

        self.export_filetype = "dat"
        self.export_table() 

    def update_table_from_tab(self): 

        self.update_table()

    def update_table(self):
        """ Updates the table """
        
        try: 
            df = self.parent.data_treat.fit_result
        except: 
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
        


