from PyQt5.QtWidgets import (QCheckBox, QPushButton, QLabel)
from PyQt5.QtGui import QFont

headline_font = QFont() #Defines headline font
headline_font.setBold(True) #Makes it bold
        
def make_headline(self, headline_string, layout): 
    self.headline = QLabel(headline_string)
    self.headline.setFont(headline_font)
    layout.addWidget(self.headline)

def make_btn(self, btn_string, function, layout): 
    self.btn = QPushButton(btn_string)
    self.btn.clicked.connect(function)
    layout.addWidget(self.btn)
    return self.btn

def make_line(self, line_string, layout): 
    self.line = QLabel(line_string)
    layout.addWidget(self.line)

def make_checkbox(self, function, layout, cb_string = ""): 
    self.checkbox = QCheckBox(cb_string) 
    self.checkbox.stateChanged.connect(function)
    layout.addWidget(self.checkbox)
    return self.checkbox
