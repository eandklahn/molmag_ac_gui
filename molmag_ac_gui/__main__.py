#std packages
import ctypes
import sys

#third-party packages
from PyQt5.QtWidgets import QApplication

#local imports
from molmag_ac_gui.ac_gui import ACGui

if __name__ == '__main__':
    
    myappid = 'AC Processing v1.0'
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    
    app = QApplication(sys.argv)
    w = ACGui()
    sys.exit(app.exec_())