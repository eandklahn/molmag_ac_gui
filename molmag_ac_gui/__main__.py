"""
Molmag AC GUI
"""

#std packages
import ctypes
import platform
import sys

#third-party packages
from PyQt5.QtWidgets import QApplication

#local imports
from . import ac_gui
from . import process_ac
from .__init__ import __version__

if __name__ == '__main__':
    
    if platform.system() == 'Windows':
        myappid = 'Molmag AC Gui {}'.format(__version__)
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    
    app = QApplication(sys.argv)
    w = ac_gui.ACGui()
    sys.exit(app.exec_())