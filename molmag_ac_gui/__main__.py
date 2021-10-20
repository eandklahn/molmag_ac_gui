#std packages
import ctypes
import sys

#third-party packages
from PyQt5.QtWidgets import QApplication

#local imports
from . import ac_gui
from . import process_ac

if __name__ == '__main__':
    
    myappid = 'AC Processing v1.0'
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    
    app = QApplication(sys.argv)
    w = ac_gui.ACGui()
    sys.exit(app.exec_())