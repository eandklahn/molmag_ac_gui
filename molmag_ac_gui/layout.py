from PyQt5.QtWidgets import (QDoubleSpinBox, QFileDialog, QListWidgetItem, 
                             QMessageBox, QWidget, QVBoxLayout, 
                             QPushButton, QLabel, QHBoxLayout, QCheckBox, 
                             QListWidget, QSplitter, QSizePolicy, QGridLayout)
from sqlalchemy import func


def make_headline(self, headline_string, layout): 
    self.headline = QLabel(headline_string)
    self.headline.setFont(self.parent.headline_font)
    layout.addWidget(self.headline)

def make_btn(self, btn_string, function, layout): 
    self.btn = QPushButton(btn_string)
    self.btn.clicked.connect(function)
    layout.addWidget(self.btn)

def make_line(self, line_string, layout): 
    self.line = QLabel(line_string)
    layout.addWidget(self.line)

