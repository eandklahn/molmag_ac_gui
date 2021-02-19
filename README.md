# Molmag AC GUI

Import `process_ac` to use in scripts with your own data.

Run the GUI with
```
python -m molmag_ac_gui
```
to use the GUI for fitting relaxation times. Your own file with data of (T, tau) can be loaded and fitted in the "Analysis"-tab
when the file is formatted as

```
Temp;Tau #header is mandatory, but content is optional
T1;tau1(;dtau1)
T2;tau2(;dtau2)
T3;tau3(;dtau3)

...
```
