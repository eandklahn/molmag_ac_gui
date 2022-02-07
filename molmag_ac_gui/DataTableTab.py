from PyQt5.QtWidgets import QLabel, QLineEdit, QPushButton, QVBoxLayout, QWidget, QTableWidget,QTableWidgetItem
import pandas as pd


class DataTableTab(QWidget): #First tab er Qwidget, så nu er det self.layout fx i stedet for self.widget.layout

    def __init__(self, parent): 
        super(DataTableTab,self).__init__() #Før __init__() peger på QWidget
        self.parent = parent
        self.initUI() 
    
    def initUI(self): 
        self.layout = QVBoxLayout() #Vertical box layout. Første ting øverst, så nedenunder osv. 
        
        self.line_edit = QLabel("When a file is loaded in the \"Data Treatment\" tab, a table of the loaded data will be shown here")

        self.tableWidget = QTableWidget() 
        self.layout.addWidget(self.line_edit)     
        self.layout.addWidget(self.tableWidget)
        self.setLayout(self.layout) 
        self.show()
    
    def updatetable(self):
        
        self.line_edit.setText("The data is loaded from the file located at: {}".format(self.parent.current_file))

        df = self.parent.raw_df
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


        