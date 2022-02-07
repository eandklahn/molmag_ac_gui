# v. 0.3.0
Published on 07/02/2022
Changed layout of Data Treatment and Data Analysis tab and adds Data Table tab. 
Makes internal changes to fitting procedures.

# v. 0.2.9
Published on 26/10/2021
Makes Windows-specific code only run on Windows

# v. 0.2.8
Published on 26/10/2021
Updated read_options so SQUID-files can be read

# v. 0.2.7
Published on 25/10/2021

# v. 0.2.6
Published on 20/10/2021
Removed requirements from setup.py.
Will instead rely on third-party packages being imported in __init__.py

# v. 0.2.5
Published on 12/5/2021
Added option to color the measured data points

# v. 0.2.4:
Published on 11/3/2021
Quick fix to show Ueff in kelvin in the fitted parameters dialog

# v. 0.2.3:
Published on 11/3/2021
Added a fit history so the last 10 fits can be retrived to see the parameters
or use the parameters as starting guesses for new fits or to use the parameters
to plot the simulation

# v. 0.2.2:
Published on 4/3/2021
Adding a data-folder to hold static variables. This will also  
make it easier in the future to add a function to change which  
header values can be read by the interface.

# v. 0.2.1:
Published on 3/3/2021
Changed fitting procedures so  
 - function to fit Xp and Xpp together is now in a separate function in process_ac
 - fitting of Xp and Xpp in the GUI is now done by multiprocessing

# v. 0.2.0:
Published on 1/3/2021
Fixed issues relating to  
 - loading of PPMS-data (issue #8)
 - loading of T,tau-data (issue #5, issue #9)
 - calculation of diamagnetic corrections (issue #6)
 - now using matplotlib to handle color mapping (issue #11)

# v. 0.1.10:
Fixed issue with loading ACMS data files (issue #4 on Github)

# v. 0.1.9:
Fixed import issue in importing process_ac

# v. 0.1.8:
Removed a dependency on my own library of scientific constants and replaced with scipy.constants
