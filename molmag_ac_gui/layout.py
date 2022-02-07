from PyQt5.QtWidgets import (QPushButton, QLabel)
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

def make_line(self, line_string, layout): 
    self.line = QLabel(line_string)
    layout.addWidget(self.line)

