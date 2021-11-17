from PyQt5.QtWidgets import QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget, QTableWidget,QTableWidgetItem
import pandas as pd


class Datatabletab(QWidget): #First tab er Qwidget, så nu er det self.layout fx i stedet for self.widget.layout

    def __init__(self, parent): 
        super(Datatabletab,self).__init__() #Før __init__() peger på QWidget
        self.parent = parent
        self.initUI() 
    
    def initUI(self): 
        self.path_to_export = ""
    
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


    def findnewpath(self): 
        #Finds path for datatable csv or excel file to be stored
        splitpath = self.parent.current_file.split('/')
        newpath = ""
        for element in splitpath[:-1]: 
            newpath += element + '/'
        newpath += splitpath[-1] + "_datatable"    
        self.path_to_export = newpath
    
    def export_table_excel(self):
        self.findnewpath()
        self.parent.raw_df.to_excel(r'{}.xlsx'.format(self.path_to_export), index=False)


    def export_table_csv(self):
        self.findnewpath()
        self.parent.raw_df.to_csv(r'{}.csv'.format(self.path_to_export), index = False)

      